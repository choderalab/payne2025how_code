"""
Create evaluators for a scaffold-filtered PLIF parquet.

Generates evaluators for all three PLIF Tversky recall cutoffs (0.5, 0.75, 1.0)
with RMSD and POSIT_Probability as scorers and date + random reference splits.
N_reference_structures values are derived from the actual filtered dataset via
logarithmic scaling.

Intended to be called by the CREATE_EVALUATORS_MODULAR Nextflow process:
    python3 create_evaluators_filtered_plif.py --output <name> --input-parquet <path>
"""

import click
import pandas as pd
from pathlib import Path
from harbor.analysis.cross_docking import DockingDataModel, EvaluatorFactory

PLIF_CUTOFFS = [0.5, 0.75, 1.0]


def make_plif_factory(cutoff: float, name: str) -> EvaluatorFactory:
    evf = EvaluatorFactory(name=name)

    evf.success_rate_evaluator_settings.use = True
    evf.success_rate_evaluator_settings.success_rate_column = "PLIFData_plif_tversky_recall"
    evf.success_rate_evaluator_settings.rmsd_cutoff = cutoff
    evf.success_rate_evaluator_settings.below_cutoff_is_good = False  # above cutoff is success

    evf.scorer_settings.rmsd_scorer_settings.use = True
    evf.scorer_settings.rmsd_scorer_settings.rmsd_column_name = "PoseData_RMSD"

    evf.scorer_settings.posit_scorer_settings.use = True
    evf.scorer_settings.posit_scorer_settings.posit_score_column_name = (
        "PoseData_docking-confidence-POSIT"
    )

    evf.reference_split_settings.use = True
    evf.reference_split_settings.date_split_settings.use = True
    evf.reference_split_settings.date_split_settings.reference_structure_date_column = (
        "RefData_Date"
    )
    evf.reference_split_settings.random_split_settings.use = True
    evf.reference_split_settings.update_reference_settings.use = True
    evf.reference_split_settings.update_reference_settings.use_logarithmic_scaling = True

    return evf


@click.command()
@click.option(
    "-i",
    "--input-parquet",
    required=True,
    type=click.Path(exists=True, path_type=Path),
)
@click.option("-o", "--output", required=True, type=Path)
def main(input_parquet, output):
    data = DockingDataModel.deserialize(input_parquet)
    n_refs = data.dataframe["Reference_Structure"].nunique()
    print(f"Loaded {len(data.dataframe)} rows, {n_refs} unique references")

    all_evaluators = []
    for cutoff in PLIF_CUTOFFS:
        name = f"filtered_plif{cutoff:g}"
        evf = make_plif_factory(cutoff, name)
        evs = evf.create_evaluators(data)
        print(f"  cutoff={cutoff}: {len(evs)} evaluators")
        all_evaluators.extend(evs)

    print(f"Total evaluators: {len(all_evaluators)}")

    output.mkdir(exist_ok=True, parents=True)
    for i, ev in enumerate(all_evaluators):
        ev.to_json_file(output / f"evaluator_{i}.json")

    summary = pd.DataFrame.from_records([ev.get_records() for ev in all_evaluators])
    summary.to_csv(output / "evaluators_summary.csv", index=False)
    print(f"Written {len(all_evaluators)} evaluator JSONs to {output}")


if __name__ == "__main__":
    main()
