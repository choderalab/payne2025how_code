"""
Create EvaluatorFactory settings files using PLIF Tversky recall as the
*scoring method* (higher recall → selected pose) instead of POSIT probability.

Generates one settings YAML per evaluation metric:
  - RMSD < 2 Å
  - PLIF Tversky recall >= 0.5
  - PLIF Tversky recall >= 0.75
  - PLIF Tversky recall >= 1.0

Each YAML covers the date + random reference split comparison only.
The RMSD and POSIT scorers are disabled so the output isolates the PLIF scorer.
"""

import click
from pathlib import Path
from harbor.analysis.cross_docking import EvaluatorFactory


EVAL_METRICS = [
    # (yaml_label, success_column, cutoff, below_cutoff_is_good)
    ("rmsd2",    "PoseData_RMSD",                    2.0,  True),
    ("plif0.5",  "PLIFData_plif_tversky_recall",     0.5,  False),
    ("plif0.75", "PLIFData_plif_tversky_recall",     0.75, False),
    ("plif1.0",  "PLIFData_plif_tversky_recall",     1.0,  False),
]


@click.command()
@click.option(
    "-o", "--output",
    type=Path,
    required=False,
    default="./",
)
def main(output):
    output.mkdir(exist_ok=True, parents=True)

    for label, success_col, cutoff, below_is_good in EVAL_METRICS:
        evf = EvaluatorFactory(name=f"reference_split_comparison_plif_scorer_{label}")

        # Evaluation metric
        evf.success_rate_evaluator_settings.use = True
        evf.success_rate_evaluator_settings.success_rate_column = success_col
        evf.success_rate_evaluator_settings.rmsd_cutoff = cutoff
        evf.success_rate_evaluator_settings.below_cutoff_is_good = below_is_good

        # PLIF scorer only — disable RMSD and POSIT
        evf.scorer_settings.rmsd_scorer_settings.use = False
        evf.scorer_settings.posit_scorer_settings.use = False
        evf.scorer_settings.plif_scorer_settings.use = True

        # Date + random reference split
        evf.reference_split_settings.use = True
        evf.reference_split_settings.date_split_settings.use = True
        evf.reference_split_settings.date_split_settings.reference_structure_date_column = (
            "RefData_Date"
        )
        evf.reference_split_settings.random_split_settings.use = True
        evf.reference_split_settings.update_reference_settings.use = True
        evf.reference_split_settings.update_reference_settings.use_logarithmic_scaling = True

        evf.to_yaml_file(output)
        print(f"Wrote {evf.name}.yaml")


if __name__ == "__main__":
    main()
