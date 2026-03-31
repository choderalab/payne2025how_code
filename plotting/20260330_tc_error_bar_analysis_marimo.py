import marimo

__generated_with = "0.21.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    from harbor.analysis.cross_docking import EvaluatorFactory, DockingDataModel
    from pathlib import Path

    return DockingDataModel, EvaluatorFactory, Path, mo


@app.cell
def _(DockingDataModel, Path):
    # Load data
    DATA_PATH = Path("/Users/apayne/Downloads/combined_docking_results")
    data = DockingDataModel.deserialize(DATA_PATH / "ALL_1_poses.json")
    return (data,)


@app.cell
def _():
    # make sure I remember how to create the correct evaluator
    return


@app.cell
def _(EvaluatorFactory, Path):
    rmsd_cutoff = 2
    output = Path("test_evaluator_factory")
    output.mkdir(exist_ok=True, parents=True)

    ef = EvaluatorFactory(name="tc_similarity")
    ef.success_rate_evaluator_settings.use = True
    ef.success_rate_evaluator_settings.success_rate_column = "PoseData_RMSD"
    ef.success_rate_evaluator_settings.rmsd_cutoff = rmsd_cutoff

    ef.scorer_settings.rmsd_scorer_settings.use = True
    ef.scorer_settings.rmsd_scorer_settings.rmsd_column_name = "PoseData_RMSD"

    ef.scorer_settings.posit_scorer_settings.use = True
    ef.scorer_settings.posit_scorer_settings.posit_score_column_name = (
        "PoseData_docking-confidence-POSIT"
    )
    sim_split = ef.__deepcopy__()
    sim_split.name = "increasing_similarity_tanimoto_combo_aligned"
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
    sim_split.pairwise_split_settings.similarity_split_settings.update_reference_settings.use = (
        True
    )
    sim_split.pairwise_split_settings.similarity_split_settings.update_reference_settings.use_logarithmic_scaling = (
        True
    )
    sim_split.to_yaml_file(output)
    return (sim_split,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Create Evaluators
    """)
    return


@app.cell
def _(data, sim_split):
    evs = sim_split.create_evaluators(data)
    return (evs,)


@app.cell
def _(data, evs):
    ev = evs[50]
    ev.n_bootstraps = 50
    results = ev.run(data)
    return ev, results


@app.cell
def _(results):
    results.ci_lower
    return


@app.cell
def _(results):
    results.ci_upper
    return


@app.cell
def _(data, ev):
    dataset_split = ev.run_dataset_split(data)
    return (dataset_split,)


@app.cell
def _(dataset_split):
    dataset_split.dataframe
    return


@app.cell
def _(dataset_split, ev):
    pose_selected = ev.run_pose_selector([dataset_split])
    return (pose_selected,)


@app.cell
def _(pose_selected):
    pose_selected[0].dataframe
    return


@app.cell
def _(ev, pose_selected):
    similarity_split = ev.run_similarity_split([pose_selected[0].copy() for n in range(10)])
    return (similarity_split,)


@app.cell
def _(similarity_split):
    similarity_split
    return


@app.cell
def _(similarity_split):
    ref_sets = [list(split.dataframe["Reference_Structure"].unique()) for split in similarity_split]
    return (ref_sets,)


@app.cell
def _(ref_sets):
    ref_sets
    return


@app.cell
def _(ev, similarity_split):
    scored = ev.run_scorer(similarity_split)
    return (scored,)


@app.cell
def _(scored):
    scored[0].dataframe
    return


@app.cell
def _(ev, scored):
    ev.evaluator.run(scored[0])
    return


@app.cell
def _(scored):
    sum(scored[0].dataframe["PoseData_RMSD"] <= 2) / len(scored[0].dataframe)
    return


@app.cell
def _(ev):
    ev.similarity_split
    return


@app.cell
def _():
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Run similarity split analysis from scratch
    """)
    return


@app.cell
def _():
    import numpy as np
    import pandas as pd

    return np, pd


@app.cell
def _(data):
    df = data.dataframe
    return (df,)


@app.cell
def _(df):
    df
    return


@app.cell
def _(df, np):
    n_bootstraps = 10
    manual_results = [] 
    for n_ref in [1, 5,10,100,403]:
        for tc_threshold in np.linspace(0,1,10):
            for i in range(n_bootstraps):
                _df = df[df["TanimotoComboData_Tanimoto"] <= tc_threshold]

                def _sample_references(group):
                    if len(group) <= n_ref:
                        return group
                    return group.sample(n=n_ref)

                def _apply_sample(df_, groupby_col):
                    grouped = df_.groupby(groupby_col).apply(
                        _sample_references, include_groups=False
                    )
                    # include_groups=False drops the groupby key; restore it from the index
                    grouped[groupby_col] = grouped.index.get_level_values(
                        0
                    )
                    return grouped.reset_index(drop=True)

                _df = _apply_sample(_df, "Query_Ligand")

                def _apply_scorer(_df):
                    return _df.sort_values("PoseData_docking-confidence-POSIT", ascending=False).groupby("Query_Ligand").head(1)

                _df = _apply_scorer(_df)

                def _calculate_success_rate(_df):
                    return sum(_df["PoseData_RMSD"] <= 2) / 403

                manual_results.append({"N_Ref": n_ref,
                                      "TC_Threshold":tc_threshold,
                                      "Success_Rate": _calculate_success_rate(_df)})
    return (manual_results,)


@app.cell
def _(df, np):
    from joblib import Parallel, delayed

    _n_bootstraps = 10
    tc_thresholds = np.linspace(0, 1, 10)
    n_refs = [1, 5, 10, 100, 403]

    # Cache filtered DFs — TC filter is identical across all bootstraps
    tc_filtered = {
      tc: df[df["TanimotoComboData_Tanimoto"] <= tc]
      for tc in tc_thresholds
    }

    def run_bootstrap(n_ref, tc_threshold, df_base):
        # Vectorized sampling: shuffle rows, then keep first n_ref per group
        shuffled = df_base.sample(frac=1)
        shuffled = shuffled.copy()
        shuffled["_rank"] = shuffled.groupby("Query_Ligand").cumcount()
        sampled = shuffled[shuffled["_rank"] < n_ref].drop(columns="_rank")

        # # Top-1 by POSIT confidence per ligand (faster than sort_values + head)
        # idx = sampled.groupby("Query_Ligand")["PoseData_docking-confidence-POSIT"].idxmax()
        # top1 = sampled.loc[idx]

        def _apply_scorer(_df):
            return _df.sort_values("PoseData_docking-confidence-POSIT", ascending=False).groupby("Query_Ligand").head(1)

        top1 = _apply_scorer(sampled)

        success_rate = (top1["PoseData_RMSD"] <= 2).sum() / 403
        return {"N_Ref": n_ref, "TC_Threshold": tc_threshold, "Success_Rate": success_rate}

    manual_results_v2 = Parallel(n_jobs=-1)(
        delayed(run_bootstrap)(n_ref, tc, tc_filtered[tc])
        for n_ref in n_refs
        for tc in tc_thresholds
        for _ in range(_n_bootstraps)
    )
    return (manual_results_v2,)


@app.cell
def _():
    return


@app.cell
def _():
    # collect results
    return


@app.cell
def _(manual_results, pd):
    manual_results_df = pd.DataFrame.from_records(manual_results)
    return


@app.cell
def _(manual_results_v2):
    manual_results_v2
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Average references available per query vs TC threshold
    """)
    return


@app.cell
def _(df, np, pd):
    _tc_thresholds = np.linspace(0, 1, 50)

    _rows = []
    for _tc in _tc_thresholds:
        _sub = df[df["TanimotoComboData_Tanimoto"] <= _tc]
        _sub = _sub.groupby(["Query_Ligand", "Reference_Ligand"]).head(1)
        _total = _sub.groupby("Query_Ligand").size()
        _good = (
            _sub[_sub["PoseData_RMSD"] <= 2]
            .groupby("Query_Ligand")
            .size()
            .reindex(_total.index, fill_value=0)
        )
        _rows.append(
            {
                "TC_Threshold": _tc * 2,  # scale to TanimotoCombo [0,2]
                "Total_Pairs": len(_sub),
                "Avg_Total_Refs": _total.mean(),
                "Std_Total_Refs": _total.std(),
                "Avg_Good_Refs": _good.mean(),
                "Std_Good_Refs": _good.std(),
            }
        )

    ref_counts_df = pd.DataFrame(_rows)
    ref_counts_df
    return (ref_counts_df,)


@app.cell
def _(ref_counts_df):
    import matplotlib.pyplot as plt
    import seaborn as sns

    fig_refs, axes = plt.subplots(1, 2, figsize=(13, 4))

    # --- Left: avg refs per query with std bands ---
    ax = axes[0]
    x = ref_counts_df["TC_Threshold"]

    for col_mean, col_std, label, color in [
        ("Avg_Total_Refs", "Std_Total_Refs", "All references", "steelblue"),
        ("Avg_Good_Refs", "Std_Good_Refs", "References with RMSD ≤ 2 Å", "coral"),
    ]:
        ax.plot(x, ref_counts_df[col_mean], label=label, color=color)
        ax.fill_between(
            x,
            ref_counts_df[col_mean] - ref_counts_df[col_std],
            ref_counts_df[col_mean] + ref_counts_df[col_std],
            color=color,
            alpha=0.2,
        )

    ax.set_xlabel("TanimotoCombo Threshold (Aligned)")
    ax.set_ylabel("Avg. references per query ligand")
    ax.set_title("References available per query vs TC threshold")
    ax.legend()
    sns.despine(ax=ax)

    # --- Right: total pairs as fraction of max possible (403 * 402) ---
    ax2 = axes[1]
    total_possible = 403 * 402
    ax2.plot(
        x,
        ref_counts_df["Total_Pairs"] / total_possible,
        color="mediumpurple",
    )
    ax2.set_xlabel("TanimotoCombo Threshold (Aligned)")
    ax2.set_ylabel("Fraction of total possible pairs")
    ax2.set_title(f"Query–reference pairs vs TC threshold\n(max = 403 × 402 = {total_possible:,})")
    ax2.set_ylim(0, 1)
    sns.despine(ax=ax2)

    plt.tight_layout()
    plt.show()
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
