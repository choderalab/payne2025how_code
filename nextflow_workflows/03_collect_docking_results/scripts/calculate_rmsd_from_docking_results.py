# Description: Calculate RMSD between ligand poses

from pathlib import Path
from tqdm import tqdm
from drugforge.data.schema.ligand import Ligand
from drugforge.data.readers.molfile import MolFileFactory
from drugforge.data.backend.openeye import oechem
import argparse
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description="Calculate RMSD between ligand poses")
    parser.add_argument(
        "-d",
        "--results_dir",
        type=Path,
        required=True,
        help="Path to directory containing docking results",
    )
    parser.add_argument(
        "-l",
        "--ligands",
        type=Path,
        required=True,
        help="Path to original ligand sdf file.",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        type=Path,
        required=False,
        default="rmsd_results.csv",
        help="Path to output file",
    )
    parser.add_argument(
        "--cutoff", type=float, default=2.0, help="RMSD cutoff for distinct poses."
    )
    return parser.parse_args()


def calculate_ligand_rmsd(ref: Ligand, fit: Ligand, append_rsmd=True) -> Ligand:
    from asapdiscovery.data.backend.openeye import oechem

    fitmol = fit.to_oemol()
    refmol = ref.to_oemol()
    nConfs = fit.num_poses
    vecRmsd = oechem.OEDoubleArray(nConfs)
    success = oechem.OERMSD(refmol, fitmol, vecRmsd)
    if not success:
        print("RMSD calculation failed")

    if append_rsmd:
        fit.set_SD_data({"RMSD": list(vecRmsd)})
    return fit


def calculate_ligand_rmsd_oemol(ref: oechem.OEMol, fit: oechem.OEMol) -> float:
    return oechem.OERMSD(ref, fit)


def get_filtered_poses(posed_ligands: list[Ligand], cutoff):
    """
    Filter out poses with RMSD above cutoff.
    Heavily based on code by Benjamin Kaminow.
    """

    # sort by pose id
    posed_ligands.sort(key=lambda x: x.tags["Pose_ID"])
    all_oemols = [posed_ligand.to_oemol() for posed_ligand in posed_ligands]
    filtered_results_idx = []
    for i, oemol1 in enumerate(all_oemols):

        # Check if this oemol is similar to any already selected
        for idx in filtered_results_idx:
            oemol2 = all_oemols[idx]
            if calculate_ligand_rmsd_oemol(oemol1, oemol2) <= cutoff:
                # Similar to an already selected mol so don't need this one
                break
        else:
            # if not similar to any in filtered_results_idx, add index to list
            filtered_results_idx.append(i)

    # pull Ligand objects from indices
    filtered_results = [posed_ligands[i] for i in filtered_results_idx]
    return filtered_results


def main():
    args = get_args()
    results_dir = args.results_dir

    out_dir = args.output_file.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    mff = MolFileFactory(filename=args.ligands)
    ligs = mff.load()
    lig_dict = {lig.compound_name: lig for lig in ligs}

    # load docked poses
    docked_poses = MolFileFactory(filename=results_dir / "docking_results.sdf").load()
    print(f"Loaded {len(docked_poses)} docked poses")

    # We don't want to filter across targets (at least at first)
    from collections import defaultdict

    pose_dict = defaultdict(list)
    _ = [
        pose_dict[lig.tags["ReferenceStructureName"]].append(lig)
        for lig in docked_poses
    ]
    filtered_results = []
    for target_name, posed_mols in pose_dict.items():
        filtered_results.extend(get_filtered_poses(posed_mols, cutoff=args.cutoff))

    print(f"Calculating RMSD for {len(filtered_results)} poses")
    records = []
    for posed_lig in tqdm(filtered_results):
        ref = lig_dict[posed_lig.compound_name]

        # no need to return anything because the posed_lig is modified directly
        calculate_ligand_rmsd(ref, posed_lig)
        records.append(
            {
                "Query_Ligand": posed_lig.compound_name,
                "Pose_ID": int(posed_lig.tags["Pose_ID"]),
                "RMSD": posed_lig.tags["RMSD"],
                "Reference_Structure": posed_lig.tags["ReferenceStructureName"],
                "Reference_Ligand": posed_lig.tags["ReferenceLigandName"],
                "docking-confidence-POSIT": posed_lig.tags["docking-confidence-POSIT"],
                "POSIT_Method": posed_lig.tags["_POSIT_method"],
                "SMILES": posed_lig.smiles,
            }
        )

    print("Writing output")
    df = pd.DataFrame.from_records(records)
    og_df = results_dir / "docking_scores_raw.csv"
    og_df = pd.read_csv(og_df)
    og_df = og_df[
        [
            "docking-structure-POSIT",
            "pose_id",
            "ligand_id",
            "docking-score-POSIT",
        ]
    ]
    og_df.columns = [
        "Reference_Structure",
        "Pose_ID",
        "Query_Ligand",
        "Docking_Score",
    ]
    # make sure Pose_ID is an int
    og_df["Pose_ID"] = og_df["Pose_ID"].astype(int)
    df["Pose_ID"] = df["Pose_ID"].astype(int)
    df = pd.merge(
        df, og_df, on=["Reference_Structure", "Pose_ID", "Query_Ligand"], how="left"
    )
    df.to_csv(args.output_file, index=False)


if __name__ == "__main__":
    main()
