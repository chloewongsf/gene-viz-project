# streamlit app: gene-expression comparator (3 cell lines)
# -------------------------------------------------------
# features:
# - upload Excel with 3 sheets (one per cell line)
# - filter by padj and |log2FC|
# - find overlap genes regulated in same direction across all 3 (or any 2 of 3)
# - heatmap of -log10(padj) (sorted by average significance)
# - heatmap of log2FC (same order)
# - volcano: all three cell lines on one plot
#   * with toggle to label top-N genes (overall or per cell line)
# - per-sheet volcanoes
# - download overlap table as CSV

import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import warnings
import io, zipfile
from datetime import datetime
from pandas import ExcelWriter

# Suppress deprecation warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# -------------------- settings --------------------
BUNDLED = Path("data/data.xlsx")  # bundled sample
SHEET_NAMES = ["GSK3 C2BBe1", "GSK3 HCT116", "GSK3 HT55"]  # reorder if needed

# -------------------- helpers --------------------
def _neglog10(x):
    x = np.asarray(x, dtype=float)
    x = np.where(x <= 0, np.nextafter(0, 1), x)
    return -np.log10(x)

def read_sheet(excel_source, sheet):
    df = pd.read_excel(excel_source, sheet_name=sheet)
    df.columns = [c.strip() for c in df.columns]
    df = df[df["regulation"].isin(["up_regulated","down_regulated"])].copy()
    df["neglog10_padj"] = _neglog10(df["padj"])
    df["direction"] = np.where(df["log2FoldChange"]>0, 1,
                        np.where(df["log2FoldChange"]<0, -1, 0))
    df["cell_line"] = sheet.replace("GSK3 ","")
    return df

def filter_table(df, padj_thr, lfc_thr):
    out = df[df["padj"] <= padj_thr].copy()
    if lfc_thr > 0:
        out = out[out["log2FoldChange"].abs() >= lfc_thr]
    return out.sort_values(["padj","log2FoldChange"], ascending=[True, False])

def overlap_same_direction(tables, require_all=True):
    names = list(tables.keys())
    slim = {
        n: tables[n][["gene_id","gene_name","log2FoldChange","padj","direction"]].rename(
            columns={"log2FoldChange": f"log2fc_{n}",
                     "padj": f"padj_{n}",
                     "direction": f"dir_{n}"}
        )
        for n in names
    }
    merged_all = slim[names[0]].merge(slim[names[1]], on="gene_id", how="inner").merge(
                 slim[names[2]], on="gene_id", how="inner")
    if require_all:
        dir_cols = [f"dir_{n}" for n in names]
        keep = (merged_all[dir_cols].sum(axis=1).abs() == len(dir_cols))
        return merged_all[keep], names
    else:
        # any 2-of-3
        outs = []
        for i in range(3):
            for j in range(i+1,3):
                a,b = names[i], names[j]
                m = slim[a].merge(slim[b], on="gene_id", how="inner")
                same = (np.sign(m[f"log2fc_{a}"]) == np.sign(m[f"log2fc_{b}"])) & (m[f"log2fc_{a}"]!=0) & (m[f"log2fc_{b}"]!=0)
                outs.append(m[same])
        out = pd.concat(outs, ignore_index=True).drop_duplicates("gene_id")
        return out, names

def build_long(df_dict):
    parts = []
    for name, df in df_dict.items():
        tmp = df[["gene_id","gene_name","log2FoldChange","padj","neglog10_padj","regulation"]].copy()
        tmp["cell_line"] = name
        parts.append(tmp)
    return pd.concat(parts, ignore_index=True)

def add_score(df, how="combo"):
    neglog = _neglog10(df["padj"])
    if how == "neglog10":
        return df.assign(score=neglog)
    return df.assign(score=neglog * df["log2FoldChange"].abs())

def volcano_all_plot(long_df, N_labels=0, per_cell=False, score_how="combo",
                     padj_thr=0.05, lfc_thr=1.0):
    fig = px.scatter(long_df, x="log2FoldChange", y="neglog10_padj",
                     color="cell_line", symbol="regulation", opacity=0.55,
                     hover_data=["gene_name","padj"])
    fig.add_vline(x= lfc_thr,  line_dash="dash", opacity=0.4)
    fig.add_vline(x=-lfc_thr, line_dash="dash", opacity=0.4)
    fig.add_hline(y=-np.log10(padj_thr), line_dash="dash", opacity=0.4)

    if N_labels and N_labels > 0:
        if per_cell:
            top = (long_df.groupby("cell_line", group_keys=False)
                         .apply(lambda g: add_score(g, score_how).nlargest(N_labels, "score")))
        else:
            top = add_score(long_df, score_how).nlargest(N_labels, "score")
        fig.add_scatter(x=top["log2FoldChange"], y=top["neglog10_padj"],
                        mode="markers+text", text=top["gene_name"],
                        textposition="top center",
                        marker=dict(size=10, line=dict(width=0.5, color="black")),
                        showlegend=False)
    fig.update_layout(title="Volcano — All Cell Lines",
                      xaxis_title="log2 fold change",
                      yaxis_title="-log10(padj)")
    return fig



# -------------------- UI --------------------
st.set_page_config(layout="wide")
st.title("Gene Expression Dashboard (3 Cell Lines)")

xl = pd.ExcelFile("data/data.xlsx")

# sidebar controls
st.sidebar.header("Controls")
padj_thr = st.sidebar.number_input("padj ≤", 0.0, 1.0, 0.05, 0.005)
lfc_thr  = st.sidebar.number_input("|log2FC| ≥ (optional)", 0.0, 10.0, 0.0, 0.1)
overlap_mode = st.sidebar.selectbox("Overlap", ["All 3 (same direction)","Any 2 of 3"])
label_N = st.sidebar.number_input("Volcano: label top-N genes", 0, 200, 30, 5)
label_per_line = st.sidebar.checkbox("Label per cell line", value=False)
label_score = st.sidebar.selectbox("Label score", ["combo (|log2FC| × −log10padj)","neglog10(padj)"])
score_key = "combo" if label_score.startswith("combo") else "neglog10"

# load & filter
raw = {s.replace("GSK3 ",""): read_sheet("data/data.xlsx", s) for s in SHEET_NAMES}
filtered = {name: filter_table(df, padj_thr, lfc_thr) for name, df in raw.items()}

st.subheader("Filtered Tables (top rows)")
for name, df in filtered.items():
    st.write(f"**{name}** — rows: {len(df)}")
    st.dataframe(df[["gene_name","log2FoldChange","padj","regulation"]].head(25), use_container_width=True)

# overlap views
require_all = overlap_mode.startswith("All")
overlap_df, order_names = overlap_same_direction(filtered, require_all=require_all)
st.markdown("---")
st.subheader(f"Overlap Genes — {overlap_mode}")
st.write(f"**Count:** {len(overlap_df)}")

# Create a more readable display with gene names
if len(overlap_df) > 0:
    # Reorder columns to put gene_name first, then gene_id, then the rest
    display_cols = ["gene_name", "gene_id"]
    # Add all other columns except gene_name and gene_id
    other_cols = [col for col in overlap_df.columns if col not in ["gene_name", "gene_id"]]
    display_cols.extend(other_cols)
    
    # Create display dataframe with reordered columns
    display_df = overlap_df[display_cols].head(100)
    st.dataframe(display_df, use_container_width=True)
else:
    st.dataframe(overlap_df.head(100), use_container_width=True)

if len(overlap_df) > 0:
    padj_cols = [f"padj_{n}" for n in order_names if f"padj_{n}" in overlap_df.columns]
    z = overlap_df.set_index("gene_id")[padj_cols].rename(columns=lambda c: c.replace("padj_",""))
    z_log = pd.DataFrame(_neglog10(z), index=z.index, columns=z.columns)
    z_log["avg"] = z_log.mean(axis=1)
    z_log = z_log.sort_values("avg", ascending=False).drop(columns="avg")

    # Add top N genes slider
    topN = st.sidebar.slider("Show top N genes in heatmaps", 5, len(z_log), min(30, len(z_log)))
    z_log_subset = z_log.head(topN)
    
    # log2FC heatmap in the same order
    fc_cols = [f"log2fc_{n}" for n in order_names if f"log2fc_{n}" in overlap_df.columns]
    fc = overlap_df.set_index("gene_id")[fc_cols].reindex(z_log.index)
    fc = fc.rename(columns=lambda c: c.replace("log2fc_",""))
    fc_subset = fc.reindex(z_log_subset.index)

    st.markdown("### Heatmap — −log10(padj) (sorted by average)")
    # Get gene names for y-axis labels
    gene_names = overlap_df.set_index("gene_id").loc[z_log_subset.index, "gene_name"]
    fig_hm = go.Figure(data=go.Heatmap(
        z=z_log_subset.values,
        x=z_log_subset.columns.tolist(),
        y=gene_names.tolist(),
        colorscale="Viridis",
        colorbar=dict(title="-log10(padj)")
    ))
    fig_hm.update_layout(height=900, title="Heatmap — −log10(padj)")
    st.plotly_chart(fig_hm, use_container_width=True, config={"displayModeBar": True})

    st.markdown("### Heatmap — log2 fold change (same order)")
    fig_fc = go.Figure(data=go.Heatmap(
        z=fc_subset.values,
        x=fc_subset.columns.tolist(),
        y=gene_names.tolist(),
        colorscale="RdBu",
        colorbar=dict(title="log2FC")
    ))
    fig_fc.update_layout(height=900, title="Heatmap — log2 fold change (same order)")
    st.plotly_chart(fig_fc, use_container_width=True, config={"displayModeBar": True})

    st.download_button("Download overlap table (CSV)",
        overlap_df.to_csv(index=False).encode(), "overlap_genes.csv", "text/csv")
else:
    st.info("No overlapping genes with current thresholds. Try relaxing filters or switch to 'Any 2 of 3'.")

# volcano (all lines)
st.markdown("---")
st.subheader("Volcano — All Cell Lines (with label toggle)")
long_df = build_long(filtered)
fig_volc = volcano_all_plot(long_df,
                            N_labels=int(label_N),
                            per_cell=label_per_line,
                            score_how=score_key,
                            padj_thr=padj_thr,
                            lfc_thr=lfc_thr if lfc_thr>0 else 1.0)
st.plotly_chart(fig_volc, use_container_width=True)

# per-cell-line volcanoes
with st.expander("Per-cell-line Volcanoes"):
    for name, df in filtered.items():
        fig = px.scatter(df, x="log2FoldChange", y="neglog10_padj",
                         color="regulation", opacity=0.65,
                         hover_data=["gene_name","padj"])
        thr = lfc_thr if lfc_thr>0 else 1.0
        fig.add_vline(x= thr, line_dash="dash", opacity=0.4)
        fig.add_vline(x=-thr, line_dash="dash", opacity=0.4)
        fig.add_hline(y=-np.log10(padj_thr), line_dash="dash", opacity=0.4)
        fig.update_layout(title=f"{name}")
        st.plotly_chart(fig, use_container_width=True)

# Downloads section
st.markdown("---")
st.header("Downloads")

# ---------- Excel Export ----------
if st.button("📊 Export Excel workbook (.xlsx)"):
    out_xlsx = io.BytesIO()
    with ExcelWriter(out_xlsx, engine="openpyxl") as xw:
        # parameters
        params = pd.DataFrame({
            "param": ["padj_threshold", "log2FC_threshold", "overlap_mode", "generated_at"],
            "value": [padj_thr, lfc_thr, "all_3_same_direction" if require_all else "any_2_of_3", f"{datetime.now():%Y-%m-%d %H:%M:%S}"]
        })
        params.to_excel(xw, index=False, sheet_name="params")

        # overlap table (if available)
        try:
            overlap_df.to_excel(xw, index=False, sheet_name="overlap_genes")
        except Exception:
            pass

        # filtered tables per cell line
        for name, df in filtered.items():
            df_out = df[["gene_id","gene_name","log2FoldChange","padj","regulation"]].copy()
            df_out.to_excel(xw, index=False, sheet_name=f"filtered_{name}")

    out_xlsx.seek(0)
    st.download_button(
        "⬇️ Download Excel",
        out_xlsx,
        file_name="genentech_results.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        use_container_width=True
    )
