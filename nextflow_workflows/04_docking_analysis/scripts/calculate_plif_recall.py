"""
Calculate PLIF recall for docked poses against their reference crystal structures.

For each pose in a docked SDF, the script:
1. Looks up the reference crystal structure in the fixed fragalysis cache
2. Runs PLIP on both (protein + crystal ligand) and (protein + docked ligand) in memory
3. Computes Tversky recall: |docked ∩ crystal| / |crystal|

Output CSV is designed to be joined back to the combined docking parquet on
(compound_name, ReferenceStructureName, Pose_ID).

Usage:
    python calculate_plif_recall.py \\
        --docked-sdf docking_results.sdf \\
        --cache-dir /path/to/mpro_fragalysis-04-01-24_curated_cache_fixed \\
        --output-csv plif_recall.csv
"""

import tempfile
import time
from pathlib import Path

import click
import pandas as pd
from drugforge.data.backend.openeye import (
    oechem,
    combine_protein_ligand,
    load_openeye_pdb,
    save_openeye_pdb,
)
from drugforge.data.readers.molfile import MolFileFactory
from harbor.analysis.utils import FileLogger

from harbor.pli.plip_analysis_schema import (
    PLIntReport,
    FingerprintLevel,
    calculate_fingerprint,
    calculate_tversky,
)



def build_cache_lookup(cache_dir: Path) -> dict[str, Path]:
    """Map structure name (e.g. 'Mpro-x11548_0A') to its cache directory."""
    lookup = {}
    for d in cache_dir.iterdir():
        if not d.is_dir():
            continue
        # Directory name format: {structure_name}-{hash}+{inchikey}
        structure_name = d.name.split("-")[0] + "-" + d.name.split("-")[1]
        lookup[structure_name] = d
    return lookup


def find_crystal_ligand_sdf(cache_entry: Path) -> Path | None:
    """Return the single SDF file in a cache directory (the crystal ligand)."""
    sdfs = list(cache_entry.glob("*.sdf"))
    if len(sdfs) == 1:
        return sdfs[0]
    if len(sdfs) > 1:
        logger.warning(f"Multiple SDFs in {cache_entry}, using first: {sdfs[0]}")
        return sdfs[0]
    return None


def plip_report_from_protein_and_ligand(
    protein_pdb_path: Path,
    ligand_mol: oechem.OEMol,
    ligand_id: str = "UNK",
) -> PLIntReport:
    """Run PLIP on a protein+ligand pair combined in memory via a temp file."""
    protein_mol = load_openeye_pdb(protein_pdb_path)
    combined = combine_protein_ligand(protein_mol, ligand_mol, lig_name=ligand_id)
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
        tmp_path = Path(f.name)
    try:
        save_openeye_pdb(combined, tmp_path)
        report = PLIntReport.from_complex_path(tmp_path, ligand_id=ligand_id)
    finally:
        tmp_path.unlink(missing_ok=True)
    return report


def compute_plif_recall(
    reference_report: PLIntReport,
    query_report: PLIntReport,
    level: FingerprintLevel,
) -> dict:
    """
    Compute Tversky recall (alpha=1, beta=0) of query vs reference at the given
    fingerprint level. Returns NaN scores if the reference has no interactions.
    """
    ref_fp = calculate_fingerprint(reference_report, level)
    query_fp = calculate_fingerprint(query_report, level)

    n_ref = sum(ref_fp.values())
    n_query = sum(query_fp.values())

    if n_ref == 0:
        return {
            "plif_tversky_recall": float("nan"),
            "plif_tanimoto": float("nan"),
            "n_reference_interactions": 0,
            "n_query_interactions": n_query,
            "n_matched_interactions": 0,
        }

    recall = calculate_tversky(ref_fp, query_fp, alpha=1, beta=0)
    tanimoto = calculate_tversky(ref_fp, query_fp, alpha=1, beta=1)

    return {
        "plif_tversky_recall": recall.score,
        "plif_tanimoto": tanimoto.score,
        "n_reference_interactions": recall.number_of_interactions_in_reference,
        "n_query_interactions": recall.number_of_interactions_in_query,
        "n_matched_interactions": recall.number_of_interactions_in_intersection,
    }


@click.command()
@click.option(
    "--docked-sdf",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Path to docking_results.sdf from a single docking job",
)
@click.option(
    "--cache-dir",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Path to the fixed fragalysis cache directory",
)
@click.option(
    "--output-csv",
    required=True,
    type=click.Path(path_type=Path),
    help="Path to output CSV file",
)
@click.option(
    "--fingerprint-level",
    type=click.Choice([lvl.name for lvl in FingerprintLevel]),
    default="ByInteractionTypeAndResidueTypeAndNumber",
    show_default=True,
    help="Fingerprint granularity for PLIF comparison",
)
@click.option(
    "--ligand-id",
    default="UNK",
    show_default=True,
    help="Residue name to assign the ligand in the combined PDB",
)
def main(docked_sdf, cache_dir, output_csv, fingerprint_level, ligand_id):
    """Calculate PLIF Tversky recall for each docked pose vs its crystal reference."""
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    logger = FileLogger(
        logname="calculate_plif_recall",
        path=output_csv.parent,
        logfile="calculate_plif_recall.log",
    ).getLogger()

    level = FingerprintLevel[fingerprint_level]

    cache_lookup = build_cache_lookup(cache_dir)
    logger.info(f"Cache lookup built: {len(cache_lookup)} entries from {cache_dir}")

    docked_poses = MolFileFactory(filename=docked_sdf).load()
    logger.info(f"Loaded {len(docked_poses)} docked poses from {docked_sdf}")

    results = []
    plip_times = []
    for posed_lig in docked_poses:
        compound_name = posed_lig.compound_name
        ref_name = posed_lig.tags["ReferenceStructureName"]
        pose_id = posed_lig.tags["Pose_ID"]

        cache_entry = cache_lookup.get(ref_name)
        if cache_entry is None:
            logger.warning(f"No cache entry for {ref_name}, skipping {compound_name}")
            continue

        protein_pdb = cache_entry / f"{ref_name}.pdb"
        if not protein_pdb.exists():
            logger.warning(f"Protein PDB not found: {protein_pdb}, skipping")
            continue

        crystal_sdf = find_crystal_ligand_sdf(cache_entry)
        if crystal_sdf is None:
            logger.warning(f"No crystal ligand SDF in {cache_entry}, skipping")
            continue

        crystal_ligs = MolFileFactory(filename=crystal_sdf).load()
        if not crystal_ligs:
            logger.warning(f"Could not read crystal ligand from {crystal_sdf}, skipping")
            continue

        try:
            t_pose = time.perf_counter()
            ref_report = plip_report_from_protein_and_ligand(protein_pdb, crystal_ligs[0].to_oemol(), ligand_id)
            query_report = plip_report_from_protein_and_ligand(protein_pdb, posed_lig.to_oemol(), ligand_id)
            scores = compute_plif_recall(ref_report, query_report, level)
            plip_times.append(time.perf_counter() - t_pose)
        except Exception as e:
            logger.warning(f"PLIP failed for {compound_name} / {ref_name}: {e}")
            scores = {
                "plif_tversky_recall": float("nan"),
                "plif_tanimoto": float("nan"),
                "n_reference_interactions": -1,
                "n_query_interactions": -1,
                "n_matched_interactions": -1,
            }

        results.append({
            "compound_name": compound_name,
            "Pose_ID": pose_id,
            "ReferenceStructureName": ref_name,
            "fingerprint_level": fingerprint_level,
            **scores,
        })

    if plip_times:
        avg = sum(plip_times) / len(plip_times)
        logger.info(f"PLIP timing: {len(plip_times)} poses, avg {avg:.3f}s/pose, total {sum(plip_times):.1f}s")

    pd.DataFrame(results).to_csv(output_csv, index=False)
    logger.info(f"Wrote {len(results)} records to {output_csv}")


if __name__ == "__main__":
    main()
