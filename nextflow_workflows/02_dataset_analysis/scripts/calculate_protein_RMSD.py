from pathlib import Path
import warnings
import click
import MDAnalysis as mda
from MDAnalysis.analysis.rms import rmsd

from protein_rmsd_schema import AtomSelection, ProteinRMSD

# MDAnalysis selection string for each AtomSelection option
ATOM_SELECTION_MAP = {
    AtomSelection.all_atom: "protein",
    AtomSelection.heavy_atom: "protein and not name H*",
    AtomSelection.c_alpha: "protein and name CA",
}

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
    atom_selection: AtomSelection = AtomSelection.heavy_atom,
) -> tuple[float, int]:
    """
    Calculate the protein RMSD between two structures after superposition.

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
        If provided, restrict alignment and RMSD to atoms of these residue
        numbers only. If None, all atoms matching atom_selection are used.
    atom_selection : AtomSelection
        Which atoms to use: all_atom, heavy_atom (default), or c_alpha.

    Returns
    -------
    float
        RMSD in Ångströms after superposition.
    """
    ref_u = mda.Universe(str(ref_pdb))
    mobile_u = mda.Universe(str(mobile_pdb))

    base_sel = ATOM_SELECTION_MAP[atom_selection]
    resid_sel = (
        " or ".join(f"resid {r}" for r in binding_site_residues)
        if binding_site_residues is not None
        else None
    )

    if resid_sel is not None:
        sel = f"({base_sel}) and chainID {chain} and ({resid_sel})"
    else:
        sel = f"({base_sel}) and chainID {chain}"

    ref_atoms = ref_u.select_atoms(sel)
    mobile_atoms = mobile_u.select_atoms(sel)

    # Ensure we only compare atoms present in both structures (e.g. missing residues).
    # Match by (resid, resname, atom name) — atom indices are universe-local and
    # cannot be compared across different Universe objects.
    ref_key_to_id = {(a.resid, a.resname, a.name): str(a.id) for a in ref_atoms}
    mobile_key_to_id = {(a.resid, a.resname, a.name): str(a.id) for a in mobile_atoms}
    common_keys = set(ref_key_to_id) & set(mobile_key_to_id)

    if not common_keys:
        raise ValueError(
            f"No common atoms found between {ref_pdb} and {mobile_pdb} "
            f"with selection '{sel}'. Check chain ID and atom selection."
        )

    ref_ids = [ref_key_to_id[k] for k in sorted(common_keys)]
    mobile_ids = [mobile_key_to_id[k] for k in sorted(common_keys)]

    ref_atoms = ref_u.select_atoms(f"id {' '.join(ref_ids)}")
    mobile_atoms = mobile_u.select_atoms(f"id {' '.join(mobile_ids)}")

    n_atoms = len(ref_ids)
    return float(rmsd(ref_atoms.positions, mobile_atoms.positions, superposition=True)), n_atoms


def _find_pdb(subdir: Path) -> Path | None:
    """Return the single PDB file in a cache subdirectory, or None."""
    pdbs = list(subdir.glob("*.pdb"))
    if len(pdbs) == 1:
        return pdbs[0]
    if len(pdbs) != 0:
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
    help="Restrict alignment and RMSD to binding-site residues only (atom type controlled by --atom_selection)",
)
@click.option(
    "--chain",
    default="A",
    show_default=True,
    help="Chain ID to use for alignment",
)
@click.option(
    "--atom_selection",
    type=click.Choice([a.value for a in AtomSelection], case_sensitive=False),
    default=AtomSelection.heavy_atom.value,
    show_default=True,
    help="Atom selection for alignment and RMSD: all_atom, heavy_atom, or c_alpha",
)
def main(
    ref_dir, mobile_dir, cache_dir, output_csv, binding_site, chain, atom_selection
):
    """Calculate pairwise protein RMSD using MDAnalysis.

    Provide either --mobile_dir for a single pairwise calculation, or --cache_dir to
    compare the reference against all structure subdirectories in the cache.
    Results are written to --output_csv with columns: Reference_Structure, Query_Structure,
    RMSD, Binding_Site_Only, Atom_Selection.
    """
    if mobile_dir is None and cache_dir is None:
        raise click.UsageError("Provide either --mobile_dir or --cache_dir.")
    if mobile_dir is not None and cache_dir is not None:
        raise click.UsageError("Provide either --mobile_dir or --cache_dir, not both.")

    bs_residues = BINDING_SITE_RESIDUES if binding_site else None
    atom_sel = AtomSelection(atom_selection)

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
        r, n_atoms = calculate_rmsd(
            ref_pdb=ref_pdb,
            mobile_pdb=mobile_pdb,
            chain=chain,
            binding_site_residues=bs_residues,
            atom_selection=atom_sel,
        )
        rows.append(
            ProteinRMSD.from_superposition(
                ref_id=ref_id,
                mobile_id=mobile_id,
                rmsd=r,
                n_atoms=n_atoms,
                binding_site_only=binding_site,
                atom_selection=atom_sel,
            )
        )

    output_path = Path(output_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ProteinRMSD.construct_dataframe(rows).to_csv(output_path, index=False)


if __name__ == "__main__":
    main()
