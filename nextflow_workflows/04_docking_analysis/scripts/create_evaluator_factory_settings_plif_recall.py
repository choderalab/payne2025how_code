"""
Create EvaluatorFactory settings files using PLIF Tversky recall as the
success metric (recall >= plif_cutoff) instead of RMSD.

Produces the same set of split configs as create_evaluator_factory_settings.py
but with a _plif<cutoff> name suffix to keep settings files distinct.
"""

import click
from pathlib import Path
from harbor.analysis.cross_docking import (
    EvaluatorFactory,
    ScaffoldSplitOptions,
)


@click.command()
@click.option(
    "-o",
    "--output",
    type=Path,
    required=False,
    default="./",
)
@click.option(
    "--plif-cutoff",
    type=float,
    default=0.5,
    help="PLIF Tversky recall threshold to define docking success (>= cutoff). Default 0.5.",
)
def main(output, plif_cutoff):
    output.mkdir(exist_ok=True, parents=True)
    suffix = f"_plif{plif_cutoff:g}"

    default = EvaluatorFactory(name="default")

    default.success_rate_evaluator_settings.use = True
    default.success_rate_evaluator_settings.success_rate_column = "PLIFData_plif_tversky_recall"
    default.success_rate_evaluator_settings.rmsd_cutoff = plif_cutoff
    default.success_rate_evaluator_settings.below_cutoff_is_good = False

    default.scorer_settings.rmsd_scorer_settings.use = True
    default.scorer_settings.rmsd_scorer_settings.rmsd_column_name = "PoseData_RMSD"

    default.scorer_settings.posit_scorer_settings.use = True
    default.scorer_settings.posit_scorer_settings.posit_score_column_name = (
        "PoseData_docking-confidence-POSIT"
    )

    # date / random reference split
    evf = default.__deepcopy__()
    evf.name = "reference_split_comparison" + suffix
    evf.reference_split_settings.use = True
    evf.reference_split_settings.date_split_settings.use = True
    evf.reference_split_settings.date_split_settings.reference_structure_date_column = (
        "RefData_Date"
    )
    evf.reference_split_settings.random_split_settings.use = True
    evf.reference_split_settings.update_reference_settings.use = True
    evf.reference_split_settings.update_reference_settings.use_logarithmic_scaling = True
    evf.to_yaml_file(output)

    # scaffold split helpers
    default_scaffold = default.__deepcopy__()
    default_scaffold.name = "default_scaffold_settings" + suffix
    default_scaffold.pairwise_split_settings.use = True
    default_scaffold.pairwise_split_settings.scaffold_split_settings.use = True
    default_scaffold.pairwise_split_settings.scaffold_split_settings.reference_scaffold_id_column = (
        "RefData_Scaffold_ID"
    )
    default_scaffold.pairwise_split_settings.scaffold_split_settings.query_scaffold_id_column = (
        "QueryData_Scaffold_ID"
    )
    default_scaffold.pairwise_split_settings.scaffold_split_settings.reference_scaffold_min_count = 1
    default_scaffold.pairwise_split_settings.scaffold_split_settings.query_scaffold_min_count = 1

    for option, name_key in [
        (ScaffoldSplitOptions.X_TO_X, "x_to_x_scaffold_split"),
        (ScaffoldSplitOptions.X_TO_Y, "x_to_y_scaffold_split"),
        (ScaffoldSplitOptions.X_TO_NOT_X, "x_to_not_x_scaffold_split"),
        (ScaffoldSplitOptions.NOT_X_TO_X, "not_x_to_x_scaffold_split"),
    ]:
        evf = default_scaffold.__deepcopy__()
        evf.name = name_key + suffix
        evf.pairwise_split_settings.scaffold_split_settings.scaffold_split_option = option
        evf.to_yaml_file(output)

    # x_to_x 5 refs
    evf = default_scaffold.__deepcopy__()
    evf.name = "x_to_x_scaffold_split_5_refs" + suffix
    evf.pairwise_split_settings.scaffold_split_settings.scaffold_split_option = ScaffoldSplitOptions.X_TO_X
    evf.dataset_before_similarity = False
    evf.combine_reference_and_similarity_splits = True
    evf.reference_split_settings.use = True
    evf.reference_split_settings.random_split_settings.use = True
    evf.reference_split_settings.n_reference_structures = [5]
    evf.pairwise_split_settings.scaffold_split_settings.reference_scaffold_min_count = 5
    evf.pairwise_split_settings.scaffold_split_settings.query_scaffold_min_count = 5
    evf.to_yaml_file(output)

    # x_to_y 5 refs
    evf = default_scaffold.__deepcopy__()
    evf.name = "x_to_y_scaffold_split_5_refs" + suffix
    evf.pairwise_split_settings.scaffold_split_settings.scaffold_split_option = ScaffoldSplitOptions.X_TO_Y
    evf.dataset_before_similarity = False
    evf.combine_reference_and_similarity_splits = True
    evf.reference_split_settings.use = True
    evf.reference_split_settings.random_split_settings.use = True
    evf.reference_split_settings.n_reference_structures = [5]
    evf.pairwise_split_settings.scaffold_split_settings.reference_scaffold_min_count = 5
    evf.pairwise_split_settings.scaffold_split_settings.query_scaffold_min_count = 5
    evf.to_yaml_file(output)

    # not_x_to_x 5 refs
    evf = default_scaffold.__deepcopy__()
    evf.name = "not_x_to_x_scaffold_split_5_refs" + suffix
    evf.pairwise_split_settings.scaffold_split_settings.scaffold_split_option = ScaffoldSplitOptions.NOT_X_TO_X
    evf.combine_reference_and_similarity_splits = True
    evf.dataset_before_similarity = False
    evf.reference_split_settings.use = True
    evf.reference_split_settings.random_split_settings.use = True
    evf.reference_split_settings.n_reference_structures = [5]
    evf.pairwise_split_settings.scaffold_split_settings.reference_scaffold_min_count = 5
    evf.pairwise_split_settings.scaffold_split_settings.query_scaffold_min_count = 5
    evf.to_yaml_file(output)

    # x_to_not_x with logarithmic N refs
    evf = default_scaffold.__deepcopy__()
    evf.name = "x_to_not_x_scaffold_split" + suffix
    evf.pairwise_split_settings.scaffold_split_settings.scaffold_split_option = ScaffoldSplitOptions.X_TO_NOT_X
    evf.dataset_before_similarity = True
    evf.combine_reference_and_similarity_splits = True
    evf.reference_split_settings.use = True
    evf.reference_split_settings.random_split_settings.use = True
    evf.reference_split_settings.update_reference_settings.use = True
    evf.reference_split_settings.update_reference_settings.use_logarithmic_scaling = True
    evf.to_yaml_file(output)

    # TC similarity split
    sim_split = default.__deepcopy__()
    sim_split.name = "increasing_similarity_tanimoto_combo_aligned" + suffix
    sim_split.pairwise_split_settings.use = True
    sim_split.pairwise_split_settings.similarity_split_settings.use = True
    sim_split.pairwise_split_settings.similarity_split_settings.similarity_column_name = (
        "TanimotoComboData_Tanimoto"
    )
    sim_split.pairwise_split_settings.similarity_split_settings.include_similar = False
    sim_split.pairwise_split_settings.similarity_split_settings.similarity_groupby_dict = {
        "TanimotoComboData_Type": "TanimotoCombo",
        "TanimotoComboData_Aligned": True,
    }
    sim_split.pairwise_split_settings.similarity_split_settings.update_reference_settings.use = True
    sim_split.pairwise_split_settings.similarity_split_settings.update_reference_settings.use_logarithmic_scaling = True
    sim_split.to_yaml_file(output)

    # MCS similarity split
    sim_split = default.__deepcopy__()
    sim_split.name = "increasing_similarity_mcs" + suffix
    sim_split.pairwise_split_settings.use = True
    sim_split.pairwise_split_settings.similarity_split_settings.use = True
    sim_split.pairwise_split_settings.similarity_split_settings.include_similar = False
    sim_split.pairwise_split_settings.similarity_split_settings.similarity_column_name = "MCSData_Tanimoto"
    sim_split.pairwise_split_settings.similarity_split_settings.similarity_groupby_dict = {
        "MCSData_Type": "MCS"
    }
    sim_split.pairwise_split_settings.similarity_split_settings.update_reference_settings.use = True
    sim_split.pairwise_split_settings.similarity_split_settings.update_reference_settings.use_logarithmic_scaling = True
    sim_split.to_yaml_file(output)

    # ECFP4 similarity split
    sim_split = default.__deepcopy__()
    sim_split.name = "increasing_similarity_ecfp4" + suffix
    sim_split.pairwise_split_settings.use = True
    sim_split.pairwise_split_settings.similarity_split_settings.use = True
    sim_split.pairwise_split_settings.similarity_split_settings.include_similar = False
    sim_split.pairwise_split_settings.similarity_split_settings.similarity_column_name = "ECFPData_Tanimoto"
    sim_split.pairwise_split_settings.similarity_split_settings.similarity_groupby_dict = {
        "ECFPData_fingerprint": "ECFP4_2048"
    }
    sim_split.pairwise_split_settings.similarity_split_settings.update_reference_settings.use = True
    sim_split.pairwise_split_settings.similarity_split_settings.update_reference_settings.use_logarithmic_scaling = True
    sim_split.to_yaml_file(output)

    print(f"PLIF recall evaluator settings written to {output} (suffix: {suffix})")


if __name__ == "__main__":
    main()
