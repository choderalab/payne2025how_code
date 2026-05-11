"""
Script to calculate the Maximum Common Substructure between reference and query ligands.
"""

from openeye import oechem
from asapdiscovery.data.readers.molfile import MolFileFactory
from asapdiscovery.data.schema.ligand import Ligand
import pandas as pd
import argparse
from pathlib import Path
from asapdiscovery.data.util.logging import FileLogger
import numpy as np
import multiprocessing as mp
from chemical_similarity_schema import MCSSimilarity


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate maximum common substructure tanimoto between ligands"
    )
    parser.add_argument(
        "--ref-ligand-sdf",
        type=Path,
        required=True,
        help="Path to directory containing prepped reference ligand sdf.",
    )
    parser.add_argument(
        "--output-dir", required=True, type=Path, help="Path to output directory"
    )
    parser.add_argument(
        "--query-ligand-sdf",
        type=Path,
        required=False,
        help="Path to directory containing prepped query ligand sdf. If false, ref-ligand-sdf will be used.",
    )
    parser.add_argument(
        "--ncpus",
        type=int,
        default=1,
        help="Number of CPUs to use for parallelization.",
    )
    return parser.parse_args()


def one_to_many_mcs(refmol: oechem.OEMol, querymols: list[oechem.OEMol]):
    """
    Get the number of atoms in the maximum common substructure and union between each pair of molecules.
    :param refmol: Reference molecule
    :param querymols: List of query molecules to compare against
    :return: Arrays of MCS atom counts and union atom counts
    """
    atomexpr = (
        oechem.OEExprOpts_Aromaticity
        | oechem.OEExprOpts_AtomicNumber
        | oechem.OEExprOpts_FormalCharge
    )
    bondexpr = oechem.OEExprOpts_Aromaticity | oechem.OEExprOpts_BondOrder

    mcs_num_atoms = np.zeros(len(querymols), dtype=int)
    union_num_atoms = np.zeros(len(querymols), dtype=int)
    pattern_query = oechem.OEQMol(refmol)
    pattern_query.BuildExpressions(atomexpr, bondexpr)
    mcss = oechem.OEMCSSearch(pattern_query)
    mcss.SetMCSFunc(oechem.OEMCSMaxAtomsCompleteCycles())

    for j, querymol in enumerate(querymols):
        try:
            mcs = next(iter(mcss.Match(querymol, True)))
            mcs_num_atoms[j] = mcs.NumAtoms()
        except StopIteration:
            mcs_num_atoms[j] = 0
        # TODO: asapdiscovery Ligand.from_oemol calls OEAddExplicitHydrogens before storing
        # the SDF string, so to_oemol() returns a mol with explicit H. NumAtoms() with no
        # predicate counts all atoms including H. Fix by calling
        # oechem.OESuppressHydrogens(mol) on refmol/querymol before this loop, or use
        # mol.NumAtoms(oechem.OEIsHeavy()) here and in the union calculation below.
        # Union = Total atoms - Overlap
        union_num_atoms[j] = refmol.NumAtoms() + querymol.NumAtoms() - mcs_num_atoms[j]

    return mcs_num_atoms, union_num_atoms


def parallelize(ref: Ligand, query_ligands: list[Ligand], logger):
    """
    Calculate the MCS between a reference ligand and a list of query ligands.
    :param ref: Reference ligand
    :param query_ligands: List of query ligands
    :return: Dataframe with MCS results
    """
    logger.info(f"Calculating MCS for {ref.compound_name}...")
    refmol = ref.to_oemol()
    query_mols = [query.to_oemol() for query in query_ligands]
    num_atoms_mcs_array, num_atoms_union_array = one_to_many_mcs(refmol, query_mols)
    tanimoto_array = num_atoms_mcs_array / num_atoms_union_array
    return MCSSimilarity.construct_dataframe(
        [
            MCSSimilarity(
                Reference_Ligand=ref.compound_name,
                Query_Ligand=query.compound_name,
                Tanimoto=tanimoto,
                N_Atoms_in_MCS=mcs,
                N_Atoms_in_Union=union,
            )
            for (query, tanimoto, mcs, union) in zip(
                query_ligands,
                tanimoto_array,
                num_atoms_mcs_array,
                num_atoms_union_array,
            )
        ]
    )


def main():
    args = parse_args()
    output_dir = args.output_dir
    output_dir.mkdir(exist_ok=True, parents=True)

    logger = FileLogger(
        "calculate_mcs_tanimoto",
        output_dir,
        logfile="calculate_mcs_tanimoto.log",
    ).getLogger()

    logger.info("Loading molecules...")
    references = MolFileFactory(filename=args.ref_ligand_sdf).load()
    queries = (
        MolFileFactory(filename=args.query_ligand_sdf).load()
        if args.query_ligand_sdf
        else references.copy()
    )

    logger.info(f"Loaded {len(references)} reference molecules.")
    logger.info(f"Loaded {len(queries)} query molecules.")

    logger.info("Calculating similarities...")
    # Parallelize the MCS calculation
    cpus = min(args.ncpus, mp.cpu_count())
    logger.info(f"Using {cpus} CPUs for parallelization.")
    with mp.Pool(cpus) as pool:
        results = pool.starmap(
            parallelize, [(ref, queries, logger) for ref in references]
        )
        results = [result for result in results if result is not None]
    # Save results
    logger.info("Saving results...")
    df = pd.concat(results, ignore_index=True)
    output_path = output_dir / "mcs_tanimoto.csv"
    df.to_csv(output_path, index=False)
    logger.info(f"Results saved to {output_path}")


if __name__ == "__main__":
    main()
