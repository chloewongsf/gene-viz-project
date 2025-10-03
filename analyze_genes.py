import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def add_significance_score(df, how="neglog10"):
    """
    how:
      - 'neglog10' -> score = -log10(padj)
      - 'combo'    -> score = -log10(padj) * abs(log2FoldChange)
    """
    neglog10 = -np.log10(df["padj"].replace(0, np.nextafter(0,1)))
    if how == "neglog10":
        df = df.assign(score=neglog10)
    elif how == "combo":
        df = df.assign(score=neglog10 * df["log2FoldChange"].abs())
    else:
        raise ValueError("how must be 'neglog10' or 'combo'")
    return df

def pick_topN(long_df, N=30, per_cell_line=False, how="neglog10"):
    """
    Returns a dataframe of the top-N rows to label.
    - If per_cell_line=True, picks top N within each cell line.
    - Else picks top N overall.
    """
    scored = add_significance_score(long_df.copy(), how=how)
    if per_cell_line:
        return (scored
                .sort_values(["cell_line","score"], ascending=[True, False])
                .groupby("cell_line", group_keys=False)
                .head(N))
    else:
        return scored.sort_values("score", ascending=False).head(N)

file_path = "data.xlsx"
sheets = ["GSK3 C2BBe1", "GSK3 HCT116", "GSK3 HT55"]

# Load all sheets into dataframes
dfs = {sheet: pd.read_excel(file_path, sheet_name=sheet) for sheet in sheets}

# Filter, sort, and clean each sheet
for sheet, df in dfs.items():
    df = df[df["regulation"].isin(["up_regulated", "down_regulated"])]
    df = df.sort_values(by=["padj", "log2FoldChange"])
    dfs[sheet] = df

# Find common genes across all 3 sheets with same regulation direction
common = dfs[sheets[0]][["gene_id", "regulation"]]
for sheet in sheets[1:]:
    common = common.merge(
        dfs[sheet][["gene_id", "regulation"]],
        on="gene_id",
        suffixes=("", f"_{sheet}")
    )

# Keep only those with consistent regulation across all 3
def consistent(row):
    return len(set(row)) == 1  # all the same direction

common = common[common.apply(lambda row: consistent(row[1:]), axis=1)]

print(f"Found {len(common)} common consistently regulated genes")

# Collect p-adjusted values for heatmap
heatmap_data = pd.DataFrame({
    sheet: dfs[sheet].set_index("gene_id").loc[common["gene_id"], "padj"]
    for sheet in sheets
}, index=common["gene_id"])

# Convert padj → -log10(padj), handling NaNs
heatmap_data_log = heatmap_data.applymap(lambda x: np.nan if pd.isna(x) else -np.log10(x))

# 1) Sort rows by average significance
row_avg = heatmap_data_log.mean(axis=1)
heatmap_data_log_sorted = heatmap_data_log.loc[row_avg.sort_values(ascending=False).index]
plt.figure(figsize=(10, 6))
sns.heatmap(heatmap_data_log_sorted, cmap="viridis", cbar_kws={"label": "-log10(padj)"})
plt.title("Common genes — sorted by average significance")
plt.xlabel("Cell Line"); plt.ylabel("Gene ID")
plt.tight_layout(); plt.savefig("heatmap_padj_sorted.png", dpi=300); plt.show()

# 2) Quick table export (for scanning)
summary = heatmap_data.copy()
summary["avg_-log10padj"] = heatmap_data.apply(lambda r: (-np.log10(r)).mean(), axis=1)
summary = summary.loc[row_avg.sort_values(ascending=False).index]
summary.to_csv("common_genes_summary.csv")

# 3) Bar chart of top N genes by average significance
TOP_N = 20
top_genes = row_avg.sort_values(ascending=False).head(TOP_N).index
stack_df = (-np.log10(heatmap_data.loc[top_genes])).reset_index().melt(id_vars="gene_id",
                                                                      var_name="cell_line",
                                                                      value_name="neglog10_padj")
plt.figure(figsize=(12, 6))
sns.barplot(data=stack_df, x="gene_id", y="neglog10_padj", hue="cell_line")
plt.xticks(rotation=75, ha="right"); plt.ylabel("-log10(padj)")
plt.title(f"Top {TOP_N} overlapping genes by average significance")
plt.tight_layout(); plt.savefig("bar_top_genes.png", dpi=300); plt.show()

# 4) Log2FC heatmap (effect sizes), same genes & order
# build a matrix of log2FC for the same overlapping genes
fc_mat = pd.DataFrame({
    s: dfs[s].set_index("gene_id").loc[heatmap_data.index, "log2FoldChange"]
    for s in sheets
})
fc_mat_sorted = fc_mat.loc[heatmap_data_log_sorted.index]  # keep same row order

plt.figure(figsize=(10, 6))
sns.heatmap(fc_mat_sorted, cmap="RdBu_r", center=0, cbar_kws={"label":"log2FC"})
plt.title("Common genes — log2 fold change (same order)")
plt.xlabel("Cell Line"); plt.ylabel("Gene ID")
plt.tight_layout(); plt.savefig("heatmap_log2fc_sorted.png", dpi=300); plt.show()

# 5) Volcano plots (one per sheet) for QA
for s in sheets:
    df = dfs[s].copy()
    df = df[df["regulation"].isin(["up_regulated","down_regulated"])]
    df["neglog10_padj"] = -np.log10(df["padj"].replace(0, np.nextafter(0,1)))
    plt.figure(figsize=(6,5))
    sns.scatterplot(data=df, x="log2FoldChange", y="neglog10_padj",
                    hue="regulation", palette={"up_regulated":"tab:orange","down_regulated":"tab:blue"},
                    alpha=0.7, s=12)
    plt.axvline(0, ls="--", c="gray", lw=1); plt.ylabel("-log10(padj)")
    plt.title(f"Volcano — {s}")
    plt.tight_layout(); plt.savefig(f"volcano_{s}.png", dpi=300); plt.show()

# 6) Combined volcano plot with all three cell lines
def prep_sheet(file_path, sheet):
    df = pd.read_excel(file_path, sheet_name=sheet)
    df = df[df["regulation"].isin(["up_regulated","down_regulated"])]
    df = df.copy()
    df["cell_line"] = sheet.replace("GSK3 ","")
    df["neglog10_padj"] = -np.log10(df["padj"].replace(0, np.nextafter(0,1)))
    return df[["gene_id","log2FoldChange","padj","neglog10_padj","regulation","cell_line"]]

long = pd.concat([prep_sheet(file_path, s) for s in sheets], ignore_index=True)

def volcano_all_matplotlib(long_df, N_labels=30, label_per_cell=False, score_how="neglog10",
                           padj_thresh=0.05, lfc_thresh=1.0, outfile_base="volcano_all"):
    # base scatter
    xmax = np.nanmax(np.abs(long_df["log2FoldChange"]))
    xlim = (-max(2, xmax), max(2, xmax))
    plt.figure(figsize=(8,6))
    ax = sns.scatterplot(
        data=long_df,
        x="log2FoldChange", y="neglog10_padj",
        hue="cell_line", style="regulation",
        alpha=0.55, s=18, edgecolor="none"
    )
    ax.axhline(-np.log10(padj_thresh), ls="--", lw=1, c="gray")
    ax.axvline( lfc_thresh,  ls="--", lw=1, c="gray")
    ax.axvline(-lfc_thresh,  ls="--", lw=1, c="gray")
    ax.set_xlim(xlim)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10(padj)")
    ax.set_title("Volcano — all cell lines")
    ax.legend(frameon=False, ncols=2)
    plt.tight_layout()
    plt.savefig(f"{outfile_base}.png", dpi=300)

    # add labels (top N)
    if N_labels and N_labels > 0:
        top = pick_topN(long_df, N=N_labels, per_cell_line=label_per_cell, how=score_how)
        for _, r in top.iterrows():
            ax.text(r["log2FoldChange"], r["neglog10_padj"], str(r["gene_id"]),
                    fontsize=7, ha="left", va="bottom", alpha=0.9)
        plt.tight_layout()
        plt.savefig(f"{outfile_base}_labeled.png", dpi=300)
    plt.show()

# Create volcano plots with different labeling strategies
# Overall top 30 by -log10(padj)
volcano_all_matplotlib(long, N_labels=30, label_per_cell=False, score_how="neglog10",
                       outfile_base="volcano_all_overall")

# Top 10 per cell line using the combo score
volcano_all_matplotlib(long, N_labels=10, label_per_cell=True, score_how="combo",
                       outfile_base="volcano_all_perline")

# 7) Interactive volcano plot (plotly)
try:
    import plotly.express as px
    
    fig = px.scatter(
        long, x="log2FoldChange", y="neglog10_padj",
        color="cell_line", symbol="regulation",
        hover_data=["gene_id","padj"],
        opacity=0.6
    )
    fig.add_vline(x=1, line_dash="dash", opacity=0.4)
    fig.add_vline(x=-1, line_dash="dash", opacity=0.4)
    fig.add_hline(y=-np.log10(0.05), line_dash="dash", opacity=0.4)
    fig.update_layout(
        title="Volcano — all cell lines",
        xaxis_title="log2 fold change",
        yaxis_title="-log10(padj)"
    )
    fig.write_image("volcano_all_cell_lines_interactive.png")  # optional static export
    fig.show()
    print("Interactive volcano plot created successfully!")
except ImportError:
    print("Plotly not available - skipping interactive plot")

# 8) Faceted volcano plots (if overplotting is an issue)
g = sns.FacetGrid(long, col="cell_line", sharex=True, sharey=True, height=4)
g.map_dataframe(sns.scatterplot, x="log2FoldChange", y="neglog10_padj",
                hue="regulation", alpha=0.6, s=15, edgecolor="none")
g.add_legend(); g.set_xlabels("log2FC"); g.set_ylabels("-log10(padj)")
plt.tight_layout(); plt.savefig("volcano_facets.png", dpi=300)
plt.show()
