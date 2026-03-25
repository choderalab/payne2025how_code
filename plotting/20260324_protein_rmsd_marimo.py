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
    fig_path = Path("./20260324_protein_rmsd")
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
def _(sns):
    sns.set_style("white")
    palette = {"Full chain": "#4878d0", "Binding site only": "#ee854a"}
    return (palette,)


@app.cell
def _(earliest_refs, palette, plt, rmsd_df, save_fig, sns):
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
        _ax.set_xlim(0, 5)
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
def _(earliest_refs, palette, plt, rmsd_df_dated, save_fig):
    fig2, axes2 = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    fig2.suptitle("Protein RMSD vs query structure date (heavy atom)", fontsize=13)
    for _ax, _scaffold_id in zip(axes2.flat, [1, 2, 3, 4]):
        _ref_structure = earliest_refs[_scaffold_id]
        _subset = rmsd_df_dated[
            (rmsd_df_dated["QueryData_Scaffold_ID"] == _scaffold_id)
            & (rmsd_df_dated["RefData_Scaffold_ID"] == _scaffold_id)
            & (rmsd_df_dated["Reference_Structure"] == _ref_structure)
        ].dropna(subset=["Query_Date"])
        for _scope, _color in palette.items():
            s = _subset[_subset["Scope"] == _scope]
            _ax.scatter(
                s["Query_Date"],
                s["ProteinRMSDData_RMSD"],
                color=_color,
                alpha=0.5,
                s=12,
                label=_scope,
            )
        _ax.set_title(f"Scaffold {_scaffold_id}\n(ref: {_ref_structure})", fontsize=10)
        _ax.set_xlabel("Query structure date")
        _ax.set_ylabel("Protein RMSD (Å)")
        _ax.set_ylim(0, 5)
        _ax.tick_params(axis="x", rotation=30)
        if _ax.get_legend():
            _ax.get_legend().remove()
    handles2 = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=palette["Full chain"],
            markersize=8,
            alpha=0.7,
            label="Full chain",
        ),
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=palette["Binding site only"],
            markersize=8,
            alpha=0.7,
            label="Binding site only",
        ),
    ]
    fig2.legend(
        handles=handles2,
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.5, -0.04),
        fontsize=9,
    )
    save_fig(fig2, "protein_rmsd_vs_date_per_scaffold")
    plt.show()
    return


if __name__ == "__main__":
    app.run()
