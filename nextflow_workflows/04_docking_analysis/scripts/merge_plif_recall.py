"""
Merge PLIF recall results (from COMBINE_PLIF_RECALL) into the main docking parquet.

Joins on (compound_name, ReferenceStructureName) and adds the plif_tversky_recall
column (and optionally plif_tanimoto, n_reference_interactions) to the parquet.

Loads and re-serializes the DockingDataModel so that downstream steps that call
DockingDataModel.deserialize() find a matching .json alongside the output .parquet.
The input .json must sit next to the input .parquet (same stem, same directory) —
Nextflow satisfies this by staging both files in the work directory.
"""

import click
import pandas as pd
from pathlib import Path
from harbor.analysis.cross_docking import DockingDataModel


@click.command()
@click.option("--input-parquet", required=True, type=Path)
@click.option("--plif-recall-csv", required=True, type=Path)
@click.option("--output-parquet", required=True, type=Path)
@click.option(
    "--plif-recall-cutoff",
    type=float,
    default=0.5,
    help="Informational only — not used in merge, just echoed to logs.",
)
def main(input_parquet, plif_recall_csv, output_parquet, plif_recall_cutoff):
    model = DockingDataModel.deserialize(input_parquet)
    df = model.dataframe
    plif = pd.read_csv(plif_recall_csv)

    print(f"Docking parquet rows: {len(df)}")
    print(f"PLIF recall rows:     {len(plif)}")
    print(f"PLIF recall columns:  {list(plif.columns)}")

    # Keep the best-fingerprint-level row per (compound_name, ReferenceStructureName)
    # calculate_plif_recall.py writes one row per fingerprint level; use the most
    # granular one (ByInteractionTypeAndResidueTypeAndNumber) if present.
    preferred_level = "ByInteractionTypeAndResidueTypeAndNumber"
    if "fingerprint_level" in plif.columns:
        plif_best = plif[plif["fingerprint_level"] == preferred_level].copy()
        if len(plif_best) == 0:
            plif_best = plif.copy()
    else:
        plif_best = plif.copy()

    keep_cols = [
        "compound_name",
        "ReferenceStructureName",
        "plif_tversky_recall",
        "plif_tanimoto",
        "n_reference_interactions",
        "n_query_interactions",
        "n_matched_interactions",
    ]
    keep_cols = [c for c in keep_cols if c in plif_best.columns]
    plif_best = plif_best[keep_cols].drop_duplicates(
        subset=["compound_name", "ReferenceStructureName"]
    )

    # Identify join keys in the docking parquet
    # Typical column names: Query_Ligand / Reference_Structure (or similar)
    ref_col = next(
        (c for c in df.columns if "Reference" in c and "Structure" in c), None
    )
    lig_col = next(
        (c for c in df.columns if "Query" in c and "Ligand" in c), None
    )
    if ref_col is None or lig_col is None:
        # Fall back: try compound_name directly
        available = list(df.columns)
        raise ValueError(
            f"Could not identify Reference/Query columns. Available: {available}"
        )

    print(f"Joining on parquet[{lig_col}] <-> plif[compound_name]")
    print(f"           parquet[{ref_col}] <-> plif[ReferenceStructureName]")

    merged = df.merge(
        plif_best.rename(columns={
            "compound_name": lig_col,
            "ReferenceStructureName": ref_col,
        }),
        on=[lig_col, ref_col],
        how="left",
    )

    n_matched = merged["plif_tversky_recall"].notna().sum()
    print(f"Rows with PLIF recall:  {n_matched} / {len(merged)}")

    # Rename PLIF value columns to PLIFData_ prefix so they are registered as a
    # proper sub-model in the DockingDataModel JSON (same pattern as ProteinRMSDData).
    plif_raw_cols = [c for c in merged.columns if c not in df.columns]
    rename_map = {c: f"PLIFData_{c}" for c in plif_raw_cols}
    merged = merged.rename(columns=rename_map)
    plif_prefixed_cols = [f"PLIFData_{c}" for c in plif_raw_cols]
    print(f"PLIFData columns registered: {plif_prefixed_cols}")

    # PLIFData joins on the same keys as PoseData (one PLIF value per pose row).
    plif_key_cols = [lig_col, ref_col]

    updated = DockingDataModel(
        dataframe=merged,
        name=model.name,
        type=model.type,
        data_types_dict={**model.data_types_dict, "PLIFData": "PoseData"},
        key_columns_dict={**model.key_columns_dict, "PLIFData": plif_key_cols},
        param_columns_dict={**model.param_columns_dict, "PLIFData": []},
        value_columns_dict={**model.value_columns_dict, "PLIFData": plif_prefixed_cols},
    )
    output_prefix = output_parquet.with_suffix("")
    updated.serialize(output_prefix)
    print(f"Written to {output_prefix}.parquet + {output_prefix}.json")


if __name__ == "__main__":
    main()
