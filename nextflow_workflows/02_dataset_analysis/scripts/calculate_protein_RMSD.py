import MDAnalysis as mda
import numpy as np
from pydantic import BaseModel, Field
from drugforge.modeling.modeling import superpose_molecule
from drugforge.data.schema.target import Target
from drugforge.modeling.schema import PreppedComplex
from drugforge.modeling.modeling import find_component_chains
from drugforge.data.backend.openeye import oechem, oespruce
import warnings
import click

# Define pocket residues
pockets = {
    "P1": [142, 141, 140, 172, 163, 143, 144],
    "P1_prime": [25,26,27],
    "P2": [41,49,54],
    "P3_4_5": [189, 190, 191, 192, 168, 167, 166, 165]
}

class BindingSite(BaseModel):
    """
    BindingSite
    """
    residues: list[int] = Field(..., description="List of residue numbers to use for alignment")

class AlignParams(BaseModel):
    """
    AlignParams
    """
    align: bool = Field(False, description="Whether to align the structures before calculating RMSD")
    ref_chain: str = Field("A", description="Chain to use for reference structure")
    mobile_chain: str = Field("A", description="Chain to use for mobile structure")


def superpose_molecule(ref_mol, mobile_mol, ref_chain="A", mobile_chain="A", binding_site: BindingSite = None):
    """
    Superpose `mobile_mol` onto `ref_mol`.

    Parameters
    ----------
    ref_mol : oechem.OEGraphMol
        Reference molecule to align to.
    mobile_mol : oechem.OEGraphMol
        Molecule to align.
    ref_chain : Reference chain to align to
    mobile_chain : Mobile chain to use for alignment (the whole molecule will move as well though)
    binding_site : BindingSite to use for alignment.

    Returns
    -------
    oechem.OEGraphMol
        New aligned molecule.
    float
        RMSD between `ref_mol` and `mobile_mol` after alignment.
    """
    chains_in_ref = find_component_chains(ref_mol, "protein", sort_by="size")
    if ref_chain not in chains_in_ref or ref_chain is None:
        warnings.warn(
            f"Chain {ref_chain} not found in reference molecule: chains {chains_in_ref}, using largest chain as reference {chains_in_ref[0]}"
        )
        ref_chain = chains_in_ref[0]

    chains_in_mobile = find_component_chains(mobile_mol, "protein", sort_by="size")
    if mobile_chain not in chains_in_mobile or mobile_chain is None:
        warnings.warn(
            f"Chain {mobile_chain} not found in mobile molecule: chains {chains_in_mobile}, using largest chain {chains_in_mobile[0]}"
        )
        mobile_chain = chains_in_mobile[0]

    if ref_chain != mobile_chain:
        warnings.warn(
            f"Chains {ref_chain} and {mobile_chain} are not the same, this may not be what you want"
        )
    ref_pred = oechem.OEHasChainID(ref_chain)
    mobile_pred = oechem.OEHasChainID(mobile_chain)

    if binding_site is not None:
        # Build an OEOrAtom predicate that matches any of the binding site residue numbers
        residues = binding_site.residues
        res_pred = oechem.OEHasResidueNumber(residues[0])
        for res_num in residues[1:]:
            res_pred = oechem.OEOrAtom(res_pred, oechem.OEHasResidueNumber(res_num))
        # AND the residue predicate with the chain predicate so we only use
        # binding-site atoms on the correct chain for each structure
        ref_pred = oechem.OEAndAtom(ref_pred, res_pred)
        mobile_pred = oechem.OEAndAtom(mobile_pred, res_pred)

    # Create object to store results
    aln_res = oespruce.OESuperposeResults()

    # Set up superposing object and set reference molecule
    superpos = oespruce.OESuperpose()
    superpos.SetupRef(ref_mol, ref_pred)

    # Perform superposing
    superpos.Superpose(aln_res, mobile_mol, mobile_pred)
    # print(f"RMSD: {aln_res.GetRMSD()}")

    # Create copy of molecule and transform it to the aligned position
    mobile_mol_aligned = mobile_mol.CreateCopy()
    aln_res.Transform(mobile_mol_aligned)
    return mobile_mol_aligned, aln_res.GetRMSD()

@click.command("calculate_protein_RMSD")
@click.option("--ref_json", type=click.Path(exists=True), required=True, help="Path to reference json file")
@click.option("--mobile_json", type=click.Path(exists=True), required=True, help="Path to mobile json file")
@click.option("--binding_site", is_flag=True, default=False, help="Use binding site residues for alignment")
def main(ref_json, mobile_json, binding_site):
    """Calculate RMSD between two protein structures, optionally using only binding site residues for alignment."""
    ref = PreppedComplex.from_json_file(ref_json)
    mobile = PreppedComplex.from_json_file(mobile_json)

    if binding_site:
        binding_site_residues = pockets["P1"] + pockets["P1_prime"] + pockets["P2"] + pockets["P3_4_5"]
        binding_site = BindingSite(residues=binding_site_residues)
    else:
        binding_site = None

    _, rmsd = superpose_molecule(ref.target.to_oemol(), mobile.target.to_oemol(), binding_site=binding_site)
    _, rmsd = superpose_molecule(ref.target.to_oemol(), mobile.target.to_oemol(), binding_site=binding_site)
    print(f"RMSD: {rmsd:.2f} Å")

if __name__ == "__main__":
    main()