from pathlib import Path
import warnings
import click
import MDAnalysis as mda
from MDAnalysis.analysis.rms import rmsd

from protein_rmsd_schema import ProteinRMSD

# Binding site residue numbers for Mpro (used when --binding_site is set)
BINDING_SITE_RESIDUES = [
    # P1
    142,
    141,
    140,
    172,
    163,
    143,
    144,
    # P1'
    25,
    26,
    27,
    # P2
    41,
    49,
    54,
    # P3-4-5
    189,
    190,
    191,
    192,
    168,
    167,
    166,
    165,
]


def calculate_rmsd(
    ref_pdb: Path,
    mobile_pdb: Path,
    chain: str = "A",
    binding_site_residues: list[int] | None = None,
) -> float:
    """
    Calculate the Cα RMSD between two protein structures after superposition.

    Uses MDAnalysis to load both structures, select the desired atoms on the
    specified chain, superpose the mobile onto the reference, and return the RMSD.

    Parameters
    ----------
    ref_pdb : Path
        Path to the reference structure PDB file.
    mobile_pdb : Path
        Path to the mobile structure PDB file.
    chain : str
        Chain ID to use for alignment (default: "A").
    binding_site_residues : list[int] or None
        If provided, restrict alignment and RMSD to Cα atoms of these residue
        numbers only. If None, all Cα atoms on the chain are used.

    Returns
    -------
    float
        Cα RMSD in Ångströms after superposition.
    """
    ref_u = mda.Universe(str(ref_pdb))
    mobile_u = mda.Universe(str(mobile_pdb))

    if binding_site_residues is not None:
        resid_sel = " or ".join(f"resid {r}" for r in binding_site_residues)
        sel = f"protein and chainID {chain} and name CA and ({resid_sel})"
    else:
        sel = f"protein and chainID {chain} and name CA"

    ref_atoms = ref_u.select_atoms(sel)
    mobile_atoms = mobile_u.select_atoms(sel)

    if len(ref_atoms) == 0:
        warnings.warn(
            f"No atoms selected in reference {ref_pdb.stem} with selection '{sel}'. "
            f"Trying without chain filter."
        )
        fallback_sel = (
            "protein and name CA"
            if binding_site_residues is None
            else (f"protein and name CA and ({resid_sel})")
        )
        ref_atoms = ref_u.select_atoms(fallback_sel)
        mobile_atoms = mobile_u.select_atoms(fallback_sel)

    if len(ref_atoms) != len(mobile_atoms):
        warnings.warn(
            f"Atom count mismatch between {ref_pdb.stem} ({len(ref_atoms)}) "
            f"and {mobile_pdb.stem} ({len(mobile_atoms)}). RMSD may be unreliable."
        )
        # trim to the smaller set by matching residue IDs present in both
        ref_resids = set(ref_atoms.resids)
        mob_resids = set(mobile_atoms.resids)
        common = sorted(ref_resids & mob_resids)
        resid_filter = " or ".join(f"resid {r}" for r in common)
        ref_atoms = ref_u.select_atoms(f"({sel}) and ({resid_filter})")
        mobile_atoms = mobile_u.select_atoms(f"({sel}) and ({resid_filter})")

    return float(rmsd(ref_atoms.positions, mobile_atoms.positions, superposition=True))


def _find_pdb(subdir: Path) -> Path | None:
    """Return the single PDB file in a cache subdirectory, or None."""
    pdbs = list(subdir.glob("*.pdb"))
    if len(pdbs) == 1:
        return pdbs[0]
    warnings.warn(f"Expected 1 PDB in {subdir}, found {len(pdbs)} — skipping.")
    return None


@click.command("calculate_protein_RMSD")
@click.option(
    "--ref_dir",
    type=click.Path(exists=True, file_okay=False),
    required=True,
    help="Path to the reference structure subdirectory (contains one PDB and one SDF)",
)
@click.option(
    "--mobile_dir",
    type=click.Path(exists=True, file_okay=False),
    default=None,
    help="Path to a single mobile structure subdirectory (use instead of --cache_dir)",
)
@click.option(
    "--cache_dir",
    type=click.Path(exists=True, file_okay=False),
    default=None,
    help="Root cache directory containing structure subdirectories to compare against the reference",
)
@click.option(
    "--output_csv",
    type=click.Path(),
    required=True,
    help="Path to output CSV file (columns: Reference_Structure, Query_Structure, RMSD, Binding_Site_Only)",
)
@click.option(
    "--binding_site",
    is_flag=True,
    default=False,
    help="Restrict alignment and RMSD to binding-site Cα atoms only",
)
@click.option(
    "--chain",
    default="A",
    show_default=True,
    help="Chain ID to use for alignment",
)
def main(ref_dir, mobile_dir, cache_dir, output_csv, binding_site, chain):
    """Calculate pairwise protein Cα RMSD using MDAnalysis.

    Provide either --mobile_dir for a single pairwise calculation, or --cache_dir to
    compare the reference against all structure subdirectories in the cache.
    Results are written to --output_csv with columns: Reference_Structure, Query_Structure,
    RMSD, Binding_Site_Only.
    """
    if mobile_dir is None and cache_dir is None:
        raise click.UsageError("Provide either --mobile_dir or --cache_dir.")
    if mobile_dir is not None and cache_dir is not None:
        raise click.UsageError("Provide either --mobile_dir or --cache_dir, not both.")

    bs_residues = BINDING_SITE_RESIDUES if binding_site else None

    # Resolve the reference PDB and ID from the ref subdir
    ref_pdb = _find_pdb(Path(ref_dir))
    if ref_pdb is None:
        raise click.UsageError(f"Could not find a PDB file in --ref_dir {ref_dir}")
    ref_id = ref_pdb.stem

    # Collect (mobile_id, mobile_pdb) pairs
    if mobile_dir is not None:
        mobile_pdb = _find_pdb(Path(mobile_dir))
        if mobile_pdb is None:
            raise click.UsageError(
                f"Could not find a PDB file in --mobile_dir {mobile_dir}"
            )
        mobile_pairs = [(mobile_pdb.stem, mobile_pdb)]
    else:
        mobile_pairs = []
        for subdir in sorted(Path(cache_dir).iterdir()):
            if not subdir.is_dir():
                continue
            pdb = _find_pdb(subdir)
            if pdb is None or pdb.stem == ref_id:
                continue
            mobile_pairs.append((pdb.stem, pdb))

    rows = []
    for mobile_id, mobile_pdb in mobile_pairs:
        r = calculate_rmsd(
            ref_pdb=ref_pdb,
            mobile_pdb=mobile_pdb,
            chain=chain,
            binding_site_residues=bs_residues,
        )
        rows.append(
            ProteinRMSD.from_superposition(
                ref_id=ref_id,
                mobile_id=mobile_id,
                rmsd=r,
                binding_site_only=binding_site,
            )
        )

    output_path = Path(output_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ProteinRMSD.construct_dataframe(rows).to_csv(output_path, index=False)


if __name__ == "__main__":
    main()
