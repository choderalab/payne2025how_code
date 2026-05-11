"""
Filter the PLIF-merged docking parquet to a scaffold-specific subset and
re-serialize as a DockingDataModel so downstream steps can deserialize it.

Filter types:
  self_docked        — ref scaffold == query scaffold, both in top-4 (IDs 0-3)
  not_top4_to_top4   — ref scaffold in top-4, query scaffold NOT in top-4
"""

import click
from pathlib import Path
from harbor.analysis.cross_docking import DockingDataModel

TOP4_SCAFFOLD_IDS = [0, 1, 2, 3]

FILTER_FUNCS = {
    "self_docked": lambda df: (
        df["RefData_Scaffold_ID"].isin(TOP4_SCAFFOLD_IDS)
        & (df["RefData_Scaffold_ID"] == df["QueryData_Scaffold_ID"])
    ),
    "not_top4_to_top4": lambda df: (
        df["RefData_Scaffold_ID"].isin(TOP4_SCAFFOLD_IDS)
        & (~df["QueryData_Scaffold_ID"].isin(TOP4_SCAFFOLD_IDS))
    ),
}


@click.command()
@click.option("--input-parquet", required=True, type=Path)
@click.option(
    "--filter-type",
    required=True,
    type=click.Choice(list(FILTER_FUNCS)),
)
@click.option("--output-stem", required=True, type=Path)
def main(input_parquet, filter_type, output_stem):
    model = DockingDataModel.deserialize(input_parquet)
    df = model.dataframe

    mask = FILTER_FUNCS[filter_type](df)
    filtered_df = df[mask].reset_index(drop=True)

    n_refs = filtered_df["Reference_Structure"].nunique()
    n_queries = filtered_df["Query_Structure"].nunique()
    print(
        f"Filter '{filter_type}': {len(filtered_df)} rows, "
        f"{n_refs} unique references, {n_queries} unique queries"
    )

    if len(filtered_df) == 0:
        raise ValueError(f"Filter '{filter_type}' produced an empty dataframe.")

    filtered_model = DockingDataModel(
        dataframe=filtered_df,
        name=model.name,
        type=model.type,
        data_types_dict=model.data_types_dict,
        key_columns_dict=model.key_columns_dict,
        param_columns_dict=model.param_columns_dict,
        value_columns_dict=model.value_columns_dict,
    )
    filtered_model.serialize(output_stem)
    print(f"Written to {output_stem}.parquet + {output_stem}.json")


if __name__ == "__main__":
    main()
