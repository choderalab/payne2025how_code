import marimo

__generated_with = "0.21.1"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # RMSD Cutoff Sensitivity — 1.5 / 2.0 / 2.5 Å (R3.7)

    Compares the main paper analyses at three RMSD success thresholds to assess
    sensitivity of conclusions to the 2.0 Å cutoff choice.

    **Data location:** `DATA_DIR` below — expects combined results CSVs with this naming:
    - `*_combined_results.csv` (2.0 Å, already collected)
    - `*_rmsd1p5_combined_results.csv` (1.5 Å, collected when pipeline finishes)
    - `*_rmsd2p5_combined_results.csv` (2.5 Å, collected when pipeline finishes)

    The `EvaluationMetric_Cutoff` column is used to tag each row's cutoff.
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
    from pathlib import Path

    return Path, pd, plt, sns


@app.cell
def _(Path):
    DATA_DIR = Path("/Users/apayne/Downloads/analyzed_results")
    return (DATA_DIR,)


@app.cell
def _(DATA_DIR):
    # Each entry: (label, cutoff_float, csv_path)
    # Add 2.5 Å paths when those files are downloaded.
    ANALYSES = {
        "datesplit": [
            (2.0, DATA_DIR / "ALL_1_poses_datesplit_combined_results.csv"),
            (1.5, DATA_DIR / "ALL_1_poses_datesplit_rmsd1p5_combined_results.csv"),
            (2.5, DATA_DIR / "ALL_1_poses_datesplit_rmsd2p5_combined_results.csv"),
        ],
        "x_to_y": [
            (2.0, DATA_DIR / "ALL_1_poses_x_to_y_combined_results.csv"),
            (1.5, DATA_DIR / "ALL_1_poses_x_to_y_rmsd1p5_combined_results.csv"),
            (2.5, DATA_DIR / "ALL_1_poses_x_to_y_rmsd2p5_combined_results.csv"),
        ],
        "x_to_x": [
            (2.0, DATA_DIR / "ALL_1_poses_x_to_x_combined_results.csv"),
            (1.5, DATA_DIR / "ALL_1_poses_x_to_x_rmsd1p5_combined_results.csv"),
            (2.5, DATA_DIR / "ALL_1_poses_x_to_x_rmsd2p5_combined_results.csv"),
        ],
    }
    return (ANALYSES,)


@app.cell
def _(ANALYSES, pd):
    def _load(analysis_name, entries):
        frames = []
        for cutoff, path in entries:
            if path.exists():
                df = pd.read_csv(path)
                df["_cutoff"] = cutoff
                df["_analysis"] = analysis_name
                frames.append(df)
                print(f"  loaded {path.name}: {len(df)} rows")
            else:
                print(f"  MISSING {path.name}")
        return pd.concat(frames, ignore_index=True) if frames else None

    all_frames = {}
    for _name, _entries in ANALYSES.items():
        print(f"=== {_name} ===")
        all_frames[_name] = _load(_name, _entries)

    CUTOFF_PALETTE = {1.5: "#e377c2", 2.0: "#1f77b4", 2.5: "#2ca02c"}
    CUTOFF_LABELS  = {1.5: "1.5 Å", 2.0: "2.0 Å (paper)", 2.5: "2.5 Å"}
    return CUTOFF_LABELS, CUTOFF_PALETTE, all_frames


@app.cell
def _(CUTOFF_LABELS, CUTOFF_PALETTE, plt, sns):
    def plot_cutoff_comparison(df, x_col, y_col="Fraction", score_filter=None,
                               ci_lower="CI_Lower", ci_upper="CI_Upper",
                               xlabel=None, ylabel="Fraction posed <RMSD threshold",
                               title="", ax=None):
        """Line plot with CI bands, one line per RMSD cutoff."""
        if score_filter is not None:
            df = df[df["Score"] == score_filter]
        if ax is None:
            fig, ax = plt.subplots(figsize=(7, 4))
        else:
            fig = ax.figure

        for cutoff, grp in df.groupby("_cutoff"):
            grp = grp.sort_values(x_col)
            color = CUTOFF_PALETTE.get(cutoff, "grey")
            label = CUTOFF_LABELS.get(cutoff, f"{cutoff} Å")
            ax.plot(grp[x_col], grp[y_col], color=color, marker="o", ms=4, label=label)
            ax.fill_between(grp[x_col], grp[ci_lower], grp[ci_upper],
                            color=color, alpha=0.15)

        ax.set_xlabel(xlabel or x_col)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(title="RMSD cutoff", fontsize=8)
        sns.despine(ax=ax)
        return fig, ax

    return (plot_cutoff_comparison,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Date split — N_Reference_Structures vs success rate
    Direct analogue of Figs 4-5. Shows how success rate grows as more
    (date-ordered) reference structures are available, across all three cutoffs.
    """)
    return


@app.cell
def _(all_frames, mo):
    ds = all_frames.get("datesplit")
    if ds is None:
        mo.stop(True, "No datesplit data loaded.")
    ds
    return (ds,)


@app.cell
def _(ds, plot_cutoff_comparison, plt, sns):
    sns.set_style("white")
    _scores = sorted(ds["Score"].unique())
    _fig_ds, _axes_ds = plt.subplots(1, len(_scores), figsize=(5 * len(_scores), 4),
                                     sharey=False, squeeze=False)
    for _ax, _score in zip(_axes_ds[0], _scores):
        plot_cutoff_comparison(
            ds, x_col="N_Reference_Structures", score_filter=_score,
            xlabel="N Reference Structures", title=f"Date split — {_score}", ax=_ax
        )
    plt.suptitle("Date-split sensitivity to RMSD cutoff", y=1.02)
    plt.tight_layout()
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Cross-scaffold (x_to_y) — aggregate success rate per scaffold-pair size

    Each point is a unique (query scaffold, reference scaffold) pair.
    Grouped and averaged by Total (ligands per pair) to get a comparable summary.
    Lines show the mean fraction per Total bin, across RMSD cutoffs.
    """)
    return


@app.cell
def _(all_frames, mo):
    xy = all_frames.get("x_to_y")
    if xy is None:
        mo.stop(True, "No x_to_y data loaded.")
    # Only keep reasonably sized pairs (Total > 1) for summary stats
    xy_agg = xy[xy["Total"] > 1].copy()
    xy_agg
    return (xy_agg,)


@app.cell
def _(CUTOFF_LABELS, CUTOFF_PALETTE, plt, sns, xy_agg):
    _scores_xy = sorted(xy_agg["Score"].unique())
    _fig_xy, _axes_xy = plt.subplots(1, len(_scores_xy), figsize=(5 * len(_scores_xy), 4),
                                     sharey=True, squeeze=False)
    sns.set_style("white")
    for _ax, _score in zip(_axes_xy[0], _scores_xy):
        for _cutoff, _grp in xy_agg[xy_agg["Score"] == _score].groupby("_cutoff"):
            _summary = _grp.groupby("Total")["Fraction"].mean().reset_index()
            _ax.plot(_summary["Total"], _summary["Fraction"],
                     color=CUTOFF_PALETTE.get(_cutoff, "grey"), marker="o", ms=4,
                     label=CUTOFF_LABELS.get(_cutoff, f"{_cutoff} Å"))
        _ax.set_xscale("log")
        _ax.set_xlabel("N query ligands in scaffold pair (log)")
        _ax.set_ylabel("Mean success rate")
        _ax.set_title(f"x_to_y — {_score}")
        _ax.legend(title="RMSD cutoff", fontsize=8)
        sns.despine(ax=_ax)
    plt.suptitle("Cross-scaffold sensitivity to RMSD cutoff", y=1.02)
    plt.tight_layout()
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Self-scaffold (x_to_x) — aggregate success rate per Total

    Same as above but for within-scaffold docking (x_to_x).
    Expected to show higher success rate than x_to_y.
    """)
    return


@app.cell
def _(all_frames, mo):
    xx = all_frames.get("x_to_x")
    if xx is None:
        mo.stop(True, "No x_to_x data loaded.")
    xx_agg = xx[xx["Total"] > 1].copy()
    xx_agg
    return (xx_agg,)


@app.cell
def _(CUTOFF_LABELS, CUTOFF_PALETTE, plt, sns, xx_agg):
    _scores_xx = sorted(xx_agg["Score"].unique())
    _fig_xx, _axes_xx = plt.subplots(1, len(_scores_xx), figsize=(5 * len(_scores_xx), 4),
                                     sharey=True, squeeze=False)
    sns.set_style("white")
    for _ax, _score in zip(_axes_xx[0], _scores_xx):
        for _cutoff, _grp in xx_agg[xx_agg["Score"] == _score].groupby("_cutoff"):
            _summary = _grp.groupby("Total")["Fraction"].mean().reset_index()
            _ax.plot(_summary["Total"], _summary["Fraction"],
                     color=CUTOFF_PALETTE.get(_cutoff, "grey"), marker="o", ms=4,
                     label=CUTOFF_LABELS.get(_cutoff, f"{_cutoff} Å"))
        _ax.set_xscale("log")
        _ax.set_xlabel("N query ligands in scaffold (log)")
        _ax.set_ylabel("Mean success rate")
        _ax.set_title(f"x_to_x — {_score}")
        _ax.legend(title="RMSD cutoff", fontsize=8)
        sns.despine(ax=_ax)
    plt.suptitle("Self-scaffold sensitivity to RMSD cutoff", y=1.02)
    plt.tight_layout()
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Numeric summary — mean success rate by cutoff and analysis
    """)
    return


@app.cell
def _(all_frames, pd):
    _rows = []
    for _analysis, _df in all_frames.items():
        if _df is None:
            continue
        _sub = _df[_df["Total"] > 10] if "Total" in _df.columns else _df
        for (_cutoff, _score), _grp in _sub.groupby(["_cutoff", "Score"]):
            _rows.append({
                "Analysis": _analysis,
                "RMSD Cutoff (Å)": _cutoff,
                "Score": _score,
                "N pairs": len(_grp),
                "Mean Fraction": _grp["Fraction"].mean().round(3),
                "Median Fraction": _grp["Fraction"].median().round(3),
            })
    summary_table = pd.DataFrame(_rows).sort_values(["Analysis", "Score", "RMSD Cutoff (Å)"])
    summary_table
    return


if __name__ == "__main__":
    app.run()
