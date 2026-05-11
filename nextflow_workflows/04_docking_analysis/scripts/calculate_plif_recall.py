"""
Calculate PLIF recall for docked poses against their original crystal structures.

For each pose in a docked SDF, the script:
1. Looks up the compound's original crystal structure(s) via structure_to_cmpd_dict.json
   (the ReferenceStructureName tag in the SDF is the *docking target*, not the crystal origin)
2. Finds which of those original structures are present in the fixed fragalysis cache
3. Runs PLIP on (original crystal protein + crystal ligand) and
              (original crystal protein + docked ligand) for each cached structure
4. Takes the best (max) Tversky recall across all valid original structures
5. Computes Tversky recall: |docked ∩ crystal| / |crystal|

Output CSV is designed to be joined back to the combined docking parquet on
(compound_name, ReferenceStructureName, Pose_ID).

Usage:
    python calculate_plif_recall.py \\
        --docked-sdf docking_results.sdf \\
        --cache-dir /path/to/mpro_fragalysis-04-01-24_curated_cache_fixed \\
        --cmpd-dict /path/to/cmpd_date_dict/structure_to_cmpd_dict.json \\
        --output-csv plif_recall.csv
"""

import json
import tempfile
import time
from collections import defaultdict
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
    """Map full structure name with suffix (e.g. 'Mpro-x11548_0A') to cache dir."""
    lookup = {}
    for d in cache_dir.iterdir():
        if not d.is_dir():
            continue
        # Directory: {structure_name}-{hash}+{inchikey}
        # structure_name includes suffix: Mpro-x11548_0A
        structure_name = d.name.split("-")[0] + "-" + d.name.split("-")[1]
        lookup[structure_name] = d
    return lookup


def build_base_to_cache_keys(cache_lookup: dict[str, Path]) -> dict[str, list[str]]:
    """
    Build a mapping from base structure name (no chain/alt suffix) to the list of
    cache keys that match it.

    cache_lookup keys look like 'Mpro-x0464_0A' or 'Mpro-P2214_0B'.
    structure_to_cmpd_dict keys look like 'Mpro-x0464' (no suffix).

    We strip the trailing _XX suffix to get the base name.
    """
    base_to_keys = defaultdict(list)
    for key in cache_lookup:
        # 'Mpro-x0464_0A' -> base 'Mpro-x0464'
        base = key.rsplit("_", 1)[0]
        base_to_keys[base].append(key)
    return dict(base_to_keys)


def build_cmpd_to_cache_entries(
    cmpd_dict_path: Path,
    cache_lookup: dict[str, Path],
    base_to_cache_keys: dict[str, list[str]],
) -> dict[str, list[tuple[str, Path]]]:
    """
    Build a mapping from compound_name to the list of (structure_name, cache_dir)
    tuples for that compound's original crystal structures that exist in the cache.
    """
    struct_to_cmpd = json.loads(cmpd_dict_path.read_text())
    cmpd_to_structs = defaultdict(list)
    for struct_base, cmpd in struct_to_cmpd.items():
        cmpd_to_structs[cmpd].append(struct_base)

    cmpd_to_entries = {}
    for cmpd, struct_bases in cmpd_to_structs.items():
        entries = []
        for base in struct_bases:
            for cache_key in base_to_cache_keys.get(base, []):
                cache_dir = cache_lookup[cache_key]
                entries.append((cache_key, cache_dir))
        cmpd_to_entries[cmpd] = entries
    return cmpd_to_entries


def find_crystal_ligand_sdf(cache_entry: Path, logger) -> Path | None:
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
    Compute Tversky recall (alpha=1, beta=0) of query vs reference.
    Returns NaN scores if the reference has no interactions.
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


def best_plif_recall_across_structures(
    posed_lig_mol: oechem.OEMol,
    crystal_structures: list[tuple[str, Path]],
    level: FingerprintLevel,
    ligand_id: str,
    logger,
) -> tuple[dict, str]:
    """
    Run PLIP against all cached original crystal structures for a compound and
    return the scores from the best-matching structure (max Tversky recall) plus
    the winning structure name.
    """
    best_scores = None
    best_struct = None

    for struct_name, cache_dir in crystal_structures:
        protein_pdb = cache_dir / f"{struct_name}.pdb"
        if not protein_pdb.exists():
            logger.warning(f"Protein PDB not found: {protein_pdb}, skipping")
            continue

        crystal_sdf = find_crystal_ligand_sdf(cache_dir, logger)
        if crystal_sdf is None:
            logger.warning(f"No crystal ligand SDF in {cache_dir}, skipping")
            continue

        crystal_ligs = MolFileFactory(filename=crystal_sdf).load()
        if not crystal_ligs:
            logger.warning(f"Could not read crystal ligand from {crystal_sdf}, skipping")
            continue

        try:
            ref_report = plip_report_from_protein_and_ligand(
                protein_pdb, crystal_ligs[0].to_oemol(), ligand_id
            )
            query_report = plip_report_from_protein_and_ligand(
                protein_pdb, posed_lig_mol, ligand_id
            )
            scores = compute_plif_recall(ref_report, query_report, level)
        except Exception as e:
            logger.warning(f"PLIP failed for {struct_name}: {e}")
            continue

        recall = scores["plif_tversky_recall"]
        if best_scores is None or (
            recall is not None
            and not pd.isna(recall)
            and (best_scores["plif_tversky_recall"] is None
                 or pd.isna(best_scores["plif_tversky_recall"])
                 or recall > best_scores["plif_tversky_recall"])
        ):
            best_scores = scores
            best_struct = struct_name

    return best_scores, best_struct


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
    "--cmpd-dict",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Path to structure_to_cmpd_dict.json (maps structure -> compound_name)",
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
)
@click.option("--ligand-id", default="UNK", show_default=True)
def main(docked_sdf, cache_dir, cmpd_dict, output_csv, fingerprint_level, ligand_id):
    """Calculate PLIF Tversky recall for each docked pose vs its original crystal structure."""
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    logger = FileLogger(
        logname="calculate_plif_recall",
        path=output_csv.parent,
        logfile="calculate_plif_recall.log",
    ).getLogger()

    level = FingerprintLevel[fingerprint_level]

    cache_lookup = build_cache_lookup(cache_dir)
    base_to_cache_keys = build_base_to_cache_keys(cache_lookup)
    cmpd_to_entries = build_cmpd_to_cache_entries(cmpd_dict, cache_lookup, base_to_cache_keys)

    logger.info(f"Cache lookup built: {len(cache_lookup)} entries from {cache_dir}")
    logger.info(f"Compound dict loaded: {len(cmpd_to_entries)} compounds")

    docked_poses = MolFileFactory(filename=docked_sdf).load()
    logger.info(f"Loaded {len(docked_poses)} docked poses from {docked_sdf}")

    results = []
    plip_times = []
    n_no_crystal = 0

    for posed_lig in docked_poses:
        compound_name = posed_lig.compound_name
        docking_ref = posed_lig.tags["ReferenceStructureName"]
        pose_id = posed_lig.tags["Pose_ID"]

        crystal_structures = cmpd_to_entries.get(compound_name, [])
        if not crystal_structures:
            logger.warning(f"No original crystal structures found for {compound_name}, skipping")
            n_no_crystal += 1
            continue

        t0 = time.perf_counter()
        scores, best_struct = best_plif_recall_across_structures(
            posed_lig.to_oemol(), crystal_structures, level, ligand_id, logger
        )
        plip_times.append(time.perf_counter() - t0)

        if scores is None:
            logger.warning(f"All PLIP attempts failed for {compound_name}")
            scores = {
                "plif_tversky_recall": float("nan"),
                "plif_tanimoto": float("nan"),
                "n_reference_interactions": -1,
                "n_query_interactions": -1,
                "n_matched_interactions": -1,
            }
            best_struct = None

        results.append({
            "compound_name": compound_name,
            "Pose_ID": pose_id,
            "ReferenceStructureName": docking_ref,
            "OriginalCrystalStructure": best_struct,
            "n_crystal_structures_checked": len(crystal_structures),
            "fingerprint_level": fingerprint_level,
            **scores,
        })

    if plip_times:
        avg = sum(plip_times) / len(plip_times)
        logger.info(
            f"PLIP timing: {len(plip_times)} poses, avg {avg:.3f}s/pose, total {sum(plip_times):.1f}s"
        )
    if n_no_crystal:
        logger.warning(f"{n_no_crystal} poses had no original crystal structure in dict")

    pd.DataFrame(results).to_csv(output_csv, index=False)
    logger.info(f"Wrote {len(results)} records to {output_csv}")


if __name__ == "__main__":
    main()
