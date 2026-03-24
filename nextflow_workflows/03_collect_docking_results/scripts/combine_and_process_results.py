"""
This script combines the results csvs and generates an input CSV that will be used to run the full cross docking evaluation.
"""

import pandas as pd
from pathlib import Path
from harbor.analysis.cross_docking import (
    DataFrameModel,
    DataFrameType,
    DockingDataModel,
)
from harbor.analysis.utils import FileLogger
import click
import numpy as np
import json


@click.command()
@click.argument("pose-data", nargs=-1, type=click.Path(exists=True), required=True)
@click.option("--tc-data", type=click.Path(exists=True))
@click.option("--ecfp-data", type=click.Path(exists=True))
@click.option("--mcs-data", type=click.Path(exists=True))
@click.option("--scaffold-data", type=click.Path(exists=True))
@click.option(
    "--protein-rmsd-data",
    type=click.Path(exists=True),
    default=None,
    help="Path to combined_protein_rmsd.csv produced by the PROTEIN_RMSD_ANALYSIS workflow",
)
@click.option(
    "--date-dict",
    type=click.Path(exists=True),
    required=True,
    help="Path to date_dict.json file",
)
@click.option(
    "--structure-cmpd-dict", required=True, help="Path to structure_to_cmpd_dict"
)
@click.option("--deduplicate/--no-deduplicate", default=False)
@click.option("--output-file-prefix", required=True, help="Output file suffix")
@click.option("--add-padding/--no-add-padding", default=True)
def main(
    pose_data,
    tc_data,
    ecfp_data,
    mcs_data,
    scaffold_data,
    protein_rmsd_data,
    date_dict,
    structure_cmpd_dict,
    deduplicate,
    output_file_prefix,
    add_padding,
):
    logger = FileLogger(
        logname="combine_and_process_results",
        path="combine_and_process_results.log",
    ).getLogger()

    dfms = []

    report_dict = {"err_msg": []}

    pose_df = pd.concat([pd.read_csv(csv) for csv in pose_data])

    # Add Date Information
    logger.info("Adding date information")
    with open(date_dict, "r") as f:
        date_dict = json.load(f)
    missing = [
        ref_structure
        for ref_structure in pose_df.Reference_Structure.unique()
        if ref_structure[:-3] not in date_dict.keys()
    ]
    if len(missing) > 0:
        report_dict["err_msg"].append(
            f"The following Reference_Structure were not in date_dict.json:"
        )
        report_dict["missing_reference_structures"] = missing

    with open(structure_cmpd_dict, "r") as f:
        structure_cmpd_dict = json.load(f)

    # combine data
    structure_cmpd_date_df = pd.DataFrame(
        structure_cmpd_dict.items(), columns=["Structure", "Ligand"]
    )
    structure_cmpd_date_df["Date"] = [
        date_dict.get(x[:-3], None) for x in structure_cmpd_date_df["Structure"]
    ]

    # Fix incorrect compound_id pulled by MetaStructureFactory from the metadata.csv
    ref_to_ligand_df = pose_df.groupby("Reference_Structure").head(1)
    incorect_ref_to_ligand_dict = {
        ref: lig
        for ref, lig in zip(
            ref_to_ligand_df.Reference_Structure, ref_to_ligand_df.Reference_Ligand
        )
    }
    incorrect_lig_to_correct_lig_dict = {
        lig: structure_cmpd_dict.get(ref[:-3])
        for ref, lig in incorect_ref_to_ligand_dict.items()
    }
    correct_ref_to_ligand_dict = {
        ref: structure_cmpd_dict.get(ref[:-3])
        for ref in pose_df.Reference_Structure.unique()
    }

    # Build reverse mapping: ligand -> structure, used to add Query_Structure column.
    # Where multiple structures share the same ligand, take the first alphabetically
    # (consistent with how reference structures are chosen elsewhere).
    ligand_to_structure_dict = {}
    for ref, lig in sorted(correct_ref_to_ligand_dict.items()):
        if lig not in ligand_to_structure_dict:
            ligand_to_structure_dict[lig] = ref

    # make an incorrect_to_correct ligand mapping
    pose_df["Reference_Ligand"] = pose_df["Reference_Ligand"].replace(
        incorrect_lig_to_correct_lig_dict
    )
    pose_df["Query_Ligand"] = pose_df["Query_Ligand"].replace(
        incorrect_lig_to_correct_lig_dict
    )

    # Add Query_Structure column derived from Query_Ligand via ligand_to_structure_dict
    pose_df["Query_Structure"] = pose_df["Query_Ligand"].map(ligand_to_structure_dict)

    # drop any query ligands that are not in the reference structures
    pose_df = pose_df[
        pose_df["Query_Ligand"].isin(pose_df["Reference_Ligand"].unique())
    ]

    # add padding to pose_df
    if add_padding:
        logger.info("Padding the data with the missing pairs")
        query_ligs = pose_df["Query_Ligand"].unique()
        ref_structures = pose_df["Reference_Structure"].unique()

        posed_pairs = set(zip(pose_df["Query_Ligand"], pose_df["Reference_Structure"]))

        from itertools import product

        possible_pairs = {
            (query, ref) for query, ref in product(query_ligs, ref_structures)
        }

        missing_pairs = possible_pairs - posed_pairs

        logger.info(f"Found {len(missing_pairs)} missing pairs to pad")

        null_df = pd.DataFrame(
            [
                {
                    "Reference_Structure": ref_struct,
                    "Query_Ligand": query_lig,
                    "Query_Structure": ligand_to_structure_dict.get(query_lig),
                    "Reference_Ligand": correct_ref_to_ligand_dict[ref_struct],
                    "RMSD": np.nan,
                    "Pose_ID": 0,
                    "POSIT_Method": "Failed",
                }
                for query_lig, ref_struct in missing_pairs
            ]
        )

        pose_df = pd.concat([pose_df, null_df])

        refs = pose_df.Reference_Ligand
        queries = pose_df.Query_Ligand
        pairs = {(ref, query) for ref, query in zip(refs, queries)}

        padding_success = len(pairs) == len(possible_pairs)

        report_dict["padding_pairs"] = [
            f"{query} - {ref}" for query, ref in missing_pairs
        ]
        report_dict["padding_success"] = padding_success
        if not padding_success:
            raise ValueError(
                f"Expected {len(possible_pairs)} pairs after padding, got {len(pairs)} pairs"
            )

    refdf = pd.DataFrame(
        {
            "Reference_Structure": list(pose_df.Reference_Structure.unique()),
            "Reference_Ligand": [
                structure_cmpd_dict.get(x[:-3], None)
                for x in pose_df.Reference_Structure.unique()
            ],
            "Date": [
                date_dict.get(x[:-3], None)
                for x in pose_df.Reference_Structure.unique()
            ],
        }
    )

    # add any with missing dates to report
    missing_dates = refdf[refdf["Date"].isnull()]
    if not missing_dates.empty:
        report_dict["err_msg"].append(
            "The following Reference_Structure had no date in date_dict.json:"
        )
        report_dict["missing_dates"] = missing_dates["Structure"].tolist()

    query_lig_set = set(pose_df["Query_Ligand"].unique())
    ref_lig_set = set(pose_df["Reference_Ligand"].unique())

    report_dict["never_docked"] = list(ref_lig_set - query_lig_set)
    report_dict["never_used_as_ref"] = list(query_lig_set - ref_lig_set)

    ref_data = DataFrameModel(
        name="RefData",
        type=DataFrameType.REFERENCE,
        dataframe=refdf,
        key_columns=["Reference_Structure", "Reference_Ligand"],
    )

    pose_dfm = DataFrameModel(
        name="PoseData",
        type=DataFrameType.POSE,
        dataframe=pose_df,
        key_columns=[
            "Reference_Structure",
            "Query_Structure",
            "Query_Ligand",
            "Reference_Ligand",
            "Pose_ID",
        ],
    )

    dfms.extend([pose_dfm])

    common_key_cols = ["Reference_Ligand", "Query_Ligand"]

    def get_dataframe(
        df_path: Path,
        deduplicate: bool,
        param_args: list = [],
    ):
        df = pd.read_csv(df_path)
        df["Query_Ligand"] = df["Query_Ligand"].apply(
            lambda x: incorrect_lig_to_correct_lig_dict.get(x, x)
        )
        df["Reference_Ligand"] = df["Reference_Ligand"].apply(
            lambda x: incorrect_lig_to_correct_lig_dict.get(x, x)
        )
        if deduplicate:
            df = df.groupby(common_key_cols + param_args).head(1)
        return df

    if tc_data:
        dfms.append(
            DataFrameModel(
                name="TanimotoComboData",
                type=DataFrameType.CHEMICAL_SIMILARITY,
                dataframe=get_dataframe(tc_data, deduplicate, ["Aligned"]),
                key_columns=common_key_cols,
                param_columns=["Aligned"],
            )
        )
    if ecfp_data:
        dfms.append(
            DataFrameModel(
                name="ECFPData",
                type=DataFrameType.CHEMICAL_SIMILARITY,
                dataframe=get_dataframe(ecfp_data, deduplicate, ["radius", "bitsize"]),
                key_columns=common_key_cols,
                param_columns=["radius", "bitsize"],
            )
        )
    if mcs_data:
        df = get_dataframe(mcs_data, deduplicate)
        dfms.append(
            DataFrameModel(
                name="MCSData",
                type=DataFrameType.CHEMICAL_SIMILARITY,
                dataframe=df,
                key_columns=common_key_cols,
            )
        )
    if protein_rmsd_data:
        logger.info("Adding protein RMSD data")
        prot_rmsd_df = pd.read_csv(protein_rmsd_data)
        dfms.append(
            DataFrameModel(
                name="ProteinRMSDData",
                type=DataFrameType.STRUCTURE_SIMILARITY,
                dataframe=prot_rmsd_df,
                key_columns=["Reference_Structure", "Query_Structure"],
                param_columns=["Binding_Site_Only", "Atom_Selection"],
            )
        )
        logger.info(f"Added protein RMSD data with {len(prot_rmsd_df)} rows")
    if scaffold_data:
        query_data = pd.read_csv(scaffold_data)
        query_data.columns = [
            "Query_Ligand",
            "Scaffold_ID",
            "Scaffold_Smarts",
            "Scaffold_Type",
        ]
        query_data["Query_Ligand"] = query_data["Query_Ligand"].apply(
            lambda x: incorrect_lig_to_correct_lig_dict.get(x, x)
        )
        dfms.append(
            DataFrameModel(
                name="QueryData",
                type=DataFrameType.QUERY,
                dataframe=(
                    query_data.groupby("Query_Ligand").head(1)
                    if deduplicate
                    else query_data
                ),
                key_columns=["Query_Ligand"],
                param_columns=["Scaffold_Type"],
            )
        )

        ref_scaffold_df = pd.read_csv(scaffold_data)
        ref_scaffold_df.columns = [
            "Reference_Ligand",
            "Scaffold_ID",
            "Scaffold_Smarts",
            "Scaffold_Type",
        ]
        ref_scaffold_df["Reference_Ligand"] = ref_scaffold_df["Reference_Ligand"].apply(
            lambda x: incorrect_lig_to_correct_lig_dict.get(x, x)
        )
        ref_scaffold_df = ref_scaffold_df.merge(
            ref_data.dataframe, on="Reference_Ligand", how="outer"
        )
        ref_data = DataFrameModel(
            name="RefData",
            type=DataFrameType.REFERENCE,
            dataframe=(
                ref_scaffold_df.groupby("Reference_Ligand").head(1)
                if deduplicate
                else ref_scaffold_df
            ),
            key_columns=["Reference_Ligand", "Reference_Structure"],
            param_columns=["Scaffold_Type"],
        )
    dfms.append(ref_data)
    ddm = DockingDataModel.from_models(dfms)
    ddm.serialize(output_file_prefix)


if __name__ == "__main__":
    main()
