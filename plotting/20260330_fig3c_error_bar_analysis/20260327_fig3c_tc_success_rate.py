import marimo

__generated_with = "0.21.1"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Fig 3C — TanimotoCombo Success Rate (minimal repro)

    Reproduces the TC similarity-split success rate plot from Figure 3C of the
    SARS-CoV-2 docking retrospective paper.

    Set `DATA_DIR` to wherever the `*_tc_combined_results.csv` files live.
    """)
    return


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    from matplotlib.colors import LogNorm
    from pathlib import Path

    return LogNorm, Path, pd, plt, sns


@app.cell
def _(Path):
    DATA_DIR = Path(
        "/Users/apayne/Downloads/analyzed_results/"
    )
    return (DATA_DIR,)


@app.cell
def _(DATA_DIR, pd):
    posit_tc = pd.read_csv(DATA_DIR / "ALL_1_poses_tc_combined_results.csv")
    posit_tc["Method"] = "POSIT"
    fred_tc = pd.read_csv(DATA_DIR / "FRED_1_poses_tc_combined_results.csv")
    fred_tc["Method"] = "FRED"
    raw = pd.concat([posit_tc, fred_tc], ignore_index=True)
    raw
    return (raw,)


@app.cell
def _(raw):
    df = raw[
        (raw["TanimotoComboData_Type"] == "TanimotoCombo")
        & (raw["TanimotoComboData_Aligned"] == True)
        & (raw["Total"] == 403)
    ].copy()
    # Scale Similarity_Threshold from [0,1] to TanimotoCombo range [0,2]
    df["TC_Threshold"] = df["Similarity_Threshold"].astype(float) * 2
    # CI error bars (clipped at 0)
    df["Error_Lower"] = (df["Fraction"] - df["CI_Lower"]).clip(lower=0)
    df["Error_Upper"] = (df["CI_Upper"] - df["Fraction"]).clip(lower=0)
    df["CI_Width"] = df["CI_Upper"] - df["CI_Lower"]
    print(f"Rows: {len(df)}")
    print(f"N_Reference_Structures values: {sorted(df['N_Reference_Structures'].unique())}")
    print(f"Score values: {df['Scoring_Method'].unique() if 'Scoring_Method' in df.columns else df['Score'].unique()}")
    df
    return (df,)


@app.cell
def _(df):
    # Quick sanity check on columns present
    df.dtypes
    return


@app.cell
def _(df):
    # Detect the score column name (may differ between runs)
    score_col = "Score" if "Score" in df.columns else "Scoring_Method"
    label_map = {
        "tc_threshold": "TanimotoCombo (Aligned)",
        "fraction": "Fraction of Ligands Posed <2Å RMSD from Crystal Pose",
        "ci_lower": "CI Lower",
        "ci_upper": "CI Upper",
        "n_ref": "Number of Reference Structures",
        "score": score_col,
        "RMSD": "RMSD (Positive Control)",
        "POSIT_Probability": "POSIT Probability",
    }
    return label_map, score_col


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Fig 3C — faceted TC success rate (no CI bands)
    """)
    return


@app.cell
def _(LogNorm, df, label_map, plt, score_col, sns):
    sns.set_style("white")
    fig_3c = sns.relplot(
        data=df,
        x="TC_Threshold",
        y="Fraction",
        hue="N_Reference_Structures",
        col="Method",
        row=score_col,
        kind="line",
        palette="viridis",
        hue_norm=LogNorm(),
        height=3,
        aspect=1.25,
        legend="full",
    )
    fig_3c.set_axis_labels(label_map["tc_threshold"], label_map["fraction"])
    fig_3c.set_titles("{col_name} | {row_name}")
    fig_3c.legend.set_title("N Refs")
    fig_3c._legend.set_bbox_to_anchor((0.775, 0.7))
    plt.suptitle("Fig 3C — TanimotoCombo success rate", y=1.02)
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## With CI bands — one N_refs at a time
    """)
    return


@app.cell
def _(df, mo):
    n_ref_options = sorted(df["N_Reference_Structures"].unique())
    n_ref_picker = mo.ui.dropdown(
        options={str(n): n for n in n_ref_options},
        value=str(n_ref_options[len(n_ref_options) // 2]),
        label="N Reference Structures",
    )
    n_ref_picker
    return (n_ref_picker,)


@app.cell
def _(df, label_map, n_ref_picker, plt, score_col, sns):
    n_sel = n_ref_picker.value
    sub = df[df["N_Reference_Structures"] == n_sel].sort_values("TC_Threshold")

    score_vals = sorted(sub[score_col].unique())
    method_vals = sorted(sub["Method"].unique())
    palette = sns.color_palette(n_colors=len(method_vals))
    method_color = dict(zip(method_vals, palette))
    style_map = dict(zip(score_vals, ["-", "--"]))

    fig_bands, ax = plt.subplots(figsize=(7, 4))
    for (method, score), grp in sub.groupby(["Method", score_col]):
        grp = grp.sort_values("TC_Threshold")
        color = method_color[method]
        ls = style_map.get(score, "-")
        ax.plot(grp["TC_Threshold"], grp["Fraction"], color=color, ls=ls, label=f"{method} / {score}")
        ax.fill_between(grp["TC_Threshold"], grp["CI_Lower"], grp["CI_Upper"], color=color, alpha=0.15)

    ax.set_xlabel(label_map["tc_threshold"])
    ax.set_ylabel(label_map["fraction"])
    ax.set_title(f"TC success rate — N refs = {n_sel}")
    ax.legend(fontsize=8, loc="lower right")
    sns.despine()
    plt.tight_layout()
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## CI width diagnostic (R2.7 / R3.9)

    If CI widths are roughly constant across N_refs, something is off.
    Expected: widths should narrow as N_refs increases (more references → less variance).
    """)
    return


@app.cell
def _(df, plt, score_col, sns):
    def plot_ci_summary():
        ci_summary = (
            df.groupby(["N_Reference_Structures", "Method", score_col])["CI_Width"]
            .mean()
            .reset_index()
            .rename(columns={"CI_Width": "Mean_CI_Width"})
        )
    
        fig_ci, ax_ci = plt.subplots(figsize=(7, 4))
        for (method, score), grp in ci_summary.groupby(["Method", score_col]):
            grp = grp.sort_values("N_Reference_Structures")
            ax_ci.plot(grp["N_Reference_Structures"], grp["Mean_CI_Width"], marker="o", label=f"{method} / {score}")
    
        ax_ci.set_xscale("log")
        ax_ci.set_xlabel("N Reference Structures (log scale)")
        ax_ci.set_ylabel("Mean CI Width (CI_Upper - CI_Lower)")
        ax_ci.set_title("Does CI width shrink with more references?")
        ax_ci.legend(fontsize=8)
        sns.despine()
        plt.tight_layout()
        plt.show()
        return ci_summary
    ci_summary = plot_ci_summary()
    return (ci_summary,)


@app.cell
def _(ci_summary):
    ci_summary.sort_values("N_Reference_Structures")
    return


if __name__ == "__main__":
    app.run()
