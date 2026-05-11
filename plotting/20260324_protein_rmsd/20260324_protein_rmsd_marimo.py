import marimo

__generated_with = "0.21.1"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import pandas as pd
    from pathlib import Path

    return Path, pd


@app.cell
def _(Path):
    fig_path = Path("/Users/apayne/science/sars-cov-2-retro-paper-asap/plotting/20260324_protein_rmsd")
    fig_path.mkdir(parents=True, exist_ok=True)

    def save_fig(fig, filename, dpi=200, suffix=".pdf"):
        figpath = Path(fig_path / f"{filename}")
        fig.savefig(figpath.with_suffix(suffix), bbox_inches="tight", dpi=dpi)

    label_map = {
        "Scope": "Protein Region",
        "Full chain": "Full Chain A",
        "Binding site only": "Binding Site Only",
        "ProteinRMSDData_RMSD": "Protein RMSD (Å)",
        "Query_Date": "Query Structure Collection Date",
    }
    return (save_fig,)


@app.cell
def _(pd):
    data = pd.read_parquet(
        "/Users/apayne/Downloads/combined_docking_results/ALL_1_poses.parquet"
    )
    return (data,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Top 4 scaffolds only
    """)
    return


@app.cell
def _(data):
    data_1 = data[
        data.QueryData_Scaffold_ID.isin([1, 2, 3, 4])
        & data.RefData_Scaffold_ID.isin([1, 2, 3, 4])
    ]
    return (data_1,)


@app.cell
def _():
    # Protein RMSD histograms — one per scaffold
    # For each scaffold, use only the earliest reference structure (by RefData_Date).
    # Deduplicate to one row per (Reference_Structure, Query_Structure, Binding_Site_Only, Atom_Selection)
    # to avoid inflating counts from the chemical-similarity parameter fan-out.
    return


@app.cell
def _(data_1, pd):
    import matplotlib.pyplot as plt
    import seaborn as sns

    data_1["RefData_Date"] = pd.to_datetime(data_1["RefData_Date"])
    # Parse dates so we can sort them
    rmsd_cols = [
        "Reference_Structure",
        "Query_Structure",
        "ProteinRMSDData_RMSD",
        "ProteinRMSDData_Binding_Site_Only",
        "ProteinRMSDData_Atom_Selection",
        "RefData_Scaffold_ID",
        "QueryData_Scaffold_ID",
        "RefData_Date",
    ]
    rmsd_df = (
        data_1[rmsd_cols]
        .drop_duplicates()
        .dropna(subset=["ProteinRMSDData_RMSD"])
        .query("ProteinRMSDData_Atom_Selection == 'heavy_atom'")
    )
    rmsd_df = rmsd_df[rmsd_df["RefData_Scaffold_ID"] == rmsd_df["QueryData_Scaffold_ID"]]
    rmsd_df["Query_Structure_Series"] = rmsd_df["Query_Structure"].apply(lambda x: x[5])
    # Deduplicate: one RMSD value per unique structure pair + RMSD type
    earliest_refs = (
        rmsd_df.sort_values("RefData_Date")
        .groupby("RefData_Scaffold_ID")["Reference_Structure"]
        .first()
    )
    print("Earliest reference structure per scaffold:")
    print(earliest_refs)
    # For each scaffold, find the single earliest reference structure by date
    # Map bool to readable label for the legend
    rmsd_df["Scope"] = rmsd_df["ProteinRMSDData_Binding_Site_Only"].map(
        {False: "Full chain", True: "Binding site only"}
    )
    return earliest_refs, plt, rmsd_df, sns


@app.cell
def _():
    max_rmsd_val = 2.1
    return (max_rmsd_val,)


@app.cell
def _(rmsd_df):
    rmsd_df
    return


@app.cell
def _(sns):
    sns.set_style("white")
    palette = {"Full chain": "#4878d0", "Binding site only": "#ee854a"}
    return (palette,)


@app.cell
def _(earliest_refs, max_rmsd_val, palette, plt, rmsd_df, save_fig, sns):
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    fig.suptitle(
        "Protein RMSD to earliest reference structure (heavy atom)", fontsize=13
    )
    for _ax, _scaffold_id in zip(axes.flat, [1, 2, 3, 4]):
        _ref_structure = earliest_refs[_scaffold_id]
        _subset = rmsd_df[
            (rmsd_df["QueryData_Scaffold_ID"] == _scaffold_id)
            & (rmsd_df["RefData_Scaffold_ID"] == _scaffold_id)
            & (rmsd_df["Reference_Structure"] == _ref_structure)
        ]
        sns.histplot(
            data=_subset,
            x="ProteinRMSDData_RMSD",
            hue="Scope",
            palette=palette,
            bins=30,
            alpha=0.6,
            element="step",
            ax=_ax,
        )
        for _scope, _color in palette.items():
            med = _subset.loc[
                _subset["Scope"] == _scope, "ProteinRMSDData_RMSD"
            ].median()
            _ax.axvline(
                med,
                color=_color,
                linestyle="--",
                linewidth=1.2,
                label=f"{_scope} median: {med:.2f} Å",
            )
        _ax.set_title(f"Scaffold {_scaffold_id}\n(ref: {_ref_structure})", fontsize=10)
        _ax.set_xlabel("Protein RMSD (Å)")
        _ax.set_ylabel("Count")
        _ax.set_xlim(0, max_rmsd_val)
        if _ax.get_legend():
            _ax.get_legend().remove()
    handles = [
        plt.Line2D(
            [0],
            [0],
            color=palette["Full chain"],
            linewidth=6,
            alpha=0.6,
            label="Full chain",
        ),
        plt.Line2D(
            [0],
            [0],
            color=palette["Binding site only"],
            linewidth=6,
            alpha=0.6,
            label="Binding site only",
        ),
        plt.Line2D(
            [0],
            [0],
            color=palette["Full chain"],
            linestyle="--",
            linewidth=1.2,
            label="Full chain median",
        ),
        plt.Line2D(
            [0],
            [0],
            color=palette["Binding site only"],
            linestyle="--",
            linewidth=1.2,
            label="Binding site median",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=4,
        bbox_to_anchor=(0.5, -0.04),
        fontsize=9,
    )
    save_fig(fig, "protein_rmsd_per_scaffold")
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # RMSD vs date of collection
    """)
    return


@app.cell
def _(rmsd_df):
    # Build a structure -> date lookup from the deduplicated rmsd_df
    # (RefData_Date is the date for whichever structure appears as Reference_Structure)
    structure_date = (
        rmsd_df[["Reference_Structure", "RefData_Date"]]
        .drop_duplicates()
        .set_index("Reference_Structure")["RefData_Date"]
    )
    # Add the query structure's collection date by mapping Query_Structure through the same lookup
    # Return a new dataframe rather than mutating rmsd_df (marimo cells must be pure)
    rmsd_df_dated = rmsd_df.copy()
    rmsd_df_dated["Query_Date"] = rmsd_df_dated["Query_Structure"].map(structure_date)
    return (rmsd_df_dated,)


@app.cell
def _(earliest_refs, max_rmsd_val, palette, plt, rmsd_df_dated, save_fig, sns):
    marker_map = {"x": "o", "P": "^"}

    fig2, axes2 = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    fig2.suptitle("Protein RMSD vs query structure date (heavy atom)", fontsize=13)
    for _ax, _scaffold_id in zip(axes2.flat, [1, 2, 3, 4]):
        _ref_structure = earliest_refs[_scaffold_id]
        _subset = rmsd_df_dated[
            (rmsd_df_dated["QueryData_Scaffold_ID"] == _scaffold_id)
            & (rmsd_df_dated["RefData_Scaffold_ID"] == _scaffold_id)
            & (rmsd_df_dated["Reference_Structure"] == _ref_structure)
        ].dropna(subset=["Query_Date"])
        sns.scatterplot(
            data=_subset,
            x="Query_Date",
            y="ProteinRMSDData_RMSD",
            hue="Scope",
            style="Query_Structure_Series",
            palette=palette,
            markers=marker_map,
            alpha=0.6,
            s=20,
            ax=_ax,
            legend=False,
        )
        _ax.set_title(f"Scaffold {_scaffold_id}\n(ref: {_ref_structure})", fontsize=10)
        _ax.set_xlabel("Query structure date")
        _ax.set_ylabel("Protein RMSD (Å)")
        _ax.set_ylim(0, max_rmsd_val)
        _ax.tick_params(axis="x", rotation=30)

    # Shared legend: colour = Scope, marker = Query_Structure_Series
    handles2 = [
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=palette["Full chain"], markersize=8, alpha=0.7,
                   label="Full chain"),
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=palette["Binding site only"], markersize=8, alpha=0.7,
                   label="Binding site only"),
        plt.Line2D([0], [0], marker="o", color="grey", markersize=8, alpha=0.7,
                   linestyle="None", label="Monoclinic (x-series)"),
        plt.Line2D([0], [0], marker="^", color="grey", markersize=8, alpha=0.7,
                   linestyle="None", label="Orthorhombic (p-series)"),
    ]
    fig2.legend(handles=handles2, loc="lower center", ncol=4,
                bbox_to_anchor=(0.5, -0.04), fontsize=9)

    save_fig(fig2, "protein_rmsd_vs_date_per_scaffold")
    plt.show()
    return


@app.cell
def _(earliest_refs, max_rmsd_val, palette, plt, rmsd_df_dated, save_fig, sns):
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

    fig3 = plt.figure(figsize=(14, 8))
    fig3.suptitle(
        "Protein RMSD vs query structure date with RMSD marginal distributions (heavy atom)",
        fontsize=13,
    )

    _outer = GridSpec(2, 2, figure=fig3, hspace=0.45, wspace=0.4)
    _marker_map = {"x": "o", "P": "^"}
    _hist_axes = []

    _all_dates = rmsd_df_dated["Query_Date"].dropna()
    _date_pad = (_all_dates.max() - _all_dates.min()) * 0.02
    _date_min = _all_dates.min() - _date_pad
    _date_max = _all_dates.max() + _date_pad

    for _i, _scaffold_id in enumerate([1, 2, 3, 4]):
        _row, _col = divmod(_i, 2)
        _inner = GridSpecFromSubplotSpec(
            1, 2, subplot_spec=_outer[_row, _col], width_ratios=[4, 1], wspace=0.05
        )
        _ax_main = fig3.add_subplot(_inner[0, 0])
        _ax_hist = fig3.add_subplot(_inner[0, 1], sharey=_ax_main)

        _ref = earliest_refs[_scaffold_id]
        _subset = rmsd_df_dated[
            (rmsd_df_dated["QueryData_Scaffold_ID"] == _scaffold_id)
            & (rmsd_df_dated["RefData_Scaffold_ID"] == _scaffold_id)
            & (rmsd_df_dated["Reference_Structure"] == _ref)
        ].dropna(subset=["Query_Date"])

        sns.scatterplot(
            data=_subset,
            x="Query_Date",
            y="ProteinRMSDData_RMSD",
            hue="Scope",
            style="Query_Structure_Series",
            palette=palette,
            markers=_marker_map,
            alpha=0.6,
            s=20,
            ax=_ax_main,
            legend=False,
        )
        _ax_main.set_title(f"Scaffold {_scaffold_id}\n(ref: {_ref})", fontsize=10)
        _ax_main.set_xlabel("Query structure date")
        _ax_main.set_ylabel("Protein RMSD (Å)")
        _ax_main.set_ylim(0, max_rmsd_val)
        _ax_main.set_xlim(_date_min, _date_max)
        _ax_main.tick_params(axis="x", rotation=30)

        _bins = [max_rmsd_val * i / 20 for i in range(21)]
        for _scope, _color in palette.items():
            _vals = _subset.loc[
                _subset["Scope"] == _scope, "ProteinRMSDData_RMSD"
            ].dropna()
            _ax_hist.hist(
                _vals,
                bins=_bins,
                orientation="horizontal",
                color=_color,
                alpha=0.7,
                density=False,
                histtype="stepfilled",
            )

        _ax_hist.set_xlabel("Count", fontsize=8)
        plt.setp(_ax_hist.get_yticklabels(), visible=False)
        _ax_hist.tick_params(axis="y", length=0)
        sns.despine(ax=_ax_hist, left=True)
        _hist_axes.append(_ax_hist)

    _max_hist_x = max(ax.get_xlim()[1] for ax in _hist_axes)
    for _ax_h in _hist_axes:
        _ax_h.set_xlim(0, _max_hist_x)

    handles3 = [
        plt.Line2D(
            [0], [0], marker="o", color="w",
            markerfacecolor=palette["Full chain"], markersize=8, alpha=0.7,
            label="Full chain",
        ),
        plt.Line2D(
            [0], [0], marker="o", color="w",
            markerfacecolor=palette["Binding site only"], markersize=8, alpha=0.7,
            label="Binding site only",
        ),
        plt.Line2D(
            [0], [0], marker="o", color="grey", markersize=8, alpha=0.7,
            linestyle="None", label="Monoclinic (x-series)",
        ),
        plt.Line2D(
            [0], [0], marker="^", color="grey", markersize=8, alpha=0.7,
            linestyle="None", label="Orthorhombic (p-series)",
        ),
    ]
    fig3.legend(
        handles=handles3, loc="lower center", ncol=4,
        bbox_to_anchor=(0.5, -0.04), fontsize=9,
    )

    save_fig(fig3, "protein_rmsd_vs_date_marginal")
    save_fig(fig3, "protein_rmsd_vs_date_marginal", suffix=".png")
    plt.show()
    return


@app.cell
def _(rmsd_df):
    rmsd_df
    return


@app.cell
def _():
    # Combine both plots
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


@app.cell
def _(rmsd_df_dated):
    high_rmsd = rmsd_df_dated[rmsd_df_dated["ProteinRMSDData_RMSD"] > 1.0].sort_values(
        "ProteinRMSDData_RMSD", ascending=False
    )
    high_rmsd
    return (high_rmsd,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Structures with Protein RMSD > 1 Å

    - **Total entries:** {len(high_rmsd)}
    - **Unique query structures:** {high_rmsd["Query_Structure"].nunique()}
    - **Unique reference structures:** {high_rmsd["Reference_Structure"].nunique()}
    - **Scaffolds represented:** {sorted(high_rmsd["QueryData_Scaffold_ID"].unique())}
    - **RMSD range:** {high_rmsd["ProteinRMSDData_RMSD"].min():.2f} – {high_rmsd["ProteinRMSDData_RMSD"].max():.2f} Å
    """)
    return


@app.cell
def _(high_rmsd):
    # Summary counts per scaffold and scope
    high_rmsd_summary = (
        high_rmsd.groupby(["QueryData_Scaffold_ID", "Scope"])
        .agg(
            count=("ProteinRMSDData_RMSD", "size"),
            mean_rmsd=("ProteinRMSDData_RMSD", "mean"),
            median_rmsd=("ProteinRMSDData_RMSD", "median"),
            max_rmsd=("ProteinRMSDData_RMSD", "max"),
        )
        .reset_index()
        .rename(columns={"QueryData_Scaffold_ID": "Scaffold"})
        .sort_values(["Scaffold", "Scope"])
    )
    high_rmsd_summary
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
