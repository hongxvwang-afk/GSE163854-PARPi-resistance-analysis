#!/usr/bin/env python3
"""DESeq2-style core-HR expression-change plot for GSE163854.

Requires pydeseq2 0.5.2. The model is fitted to all non-zero genes in the raw
count matrix; BH-adjusted P values therefore reflect genome-wide testing.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


ROOT = Path(__file__).resolve().parents[1]
COUNTS = ROOT / "data" / "raw" / "GSE163854_counts.txt.gz"
OUT_STEM = ROOT / "figures" / "GSE163854_core_HR_DESeq2_barplot_final"

SENSITIVE = ["P2A1", "P2A2", "P2A3", "F1A3", "F1A4"]
RESISTANT = ["F351", "F352", "F451", "F452", "F453"]
SAMPLES = SENSITIVE + RESISTANT
CORE_HR = [
    "BLM", "TOP3A", "RMI1", "RMI2", "BRCA1", "BRCA2", "PALB2",
    "BARD1", "BRIP1", "RAD51", "RAD51B", "RAD51C", "RAD51D",
    "XRCC2", "XRCC3",
]
TOP_N = 8

UP = "#C65D32"
DOWN = "#2E6F93"
INK = "#1F272B"
MUTED = "#5B666B"
GRID = "#D9E0E3"


def significance_label(fdr: float) -> str:
    if pd.isna(fdr) or fdr >= 0.05:
        return "ns"
    if fdr < 0.001:
        return "***"
    if fdr < 0.01:
        return "**"
    return "*"


def run_deseq2() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(COUNTS, sep="\t", index_col=0)
    missing_samples = [sample for sample in SAMPLES if sample not in raw.columns]
    if missing_samples:
        raise ValueError(f"Missing samples: {missing_samples}")

    counts = raw.loc[:, SAMPLES]
    if counts.index.duplicated().any():
        counts = counts.groupby(level=0).sum()
    counts = counts.loc[counts.sum(axis=1).gt(0)].T.astype(int)

    metadata = pd.DataFrame(
        {"condition": ["sensitive"] * len(SENSITIVE) + ["resistant"] * len(RESISTANT)},
        index=SAMPLES,
    )
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design="~condition",
        refit_cooks=True,
        n_cpus=4,
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(
        dds,
        contrast=["condition", "resistant", "sensitive"],
        alpha=0.05,
        cooks_filter=True,
        independent_filter=True,
        n_cpus=4,
        quiet=True,
    )
    stats.summary()

    results = stats.results_df.rename(
        columns={
            "log2FoldChange": "resistant_vs_sensitive_log2FC",
            "lfcSE": "log2FC_SE",
            "stat": "wald_statistic",
            "pvalue": "wald_p_value",
            "padj": "BH_FDR",
        }
    )
    results.index.name = "gene"
    missing_genes = [gene for gene in CORE_HR if gene not in results.index]
    if missing_genes:
        raise ValueError(f"Missing core-HR genes: {missing_genes}")

    core = results.loc[CORE_HR, [
        "baseMean", "resistant_vs_sensitive_log2FC", "log2FC_SE",
        "wald_statistic", "wald_p_value", "BH_FDR",
    ]].reset_index()
    core["absolute_log2FC"] = core["resistant_vs_sensitive_log2FC"].abs()
    core["direction"] = np.where(
        core["resistant_vs_sensitive_log2FC"].ge(0),
        "Higher in resistant",
        "Lower in resistant",
    )
    core["significance"] = core["BH_FDR"].map(significance_label)
    core["nominal_p_lt_0.05"] = core["wald_p_value"].lt(0.05)
    core["BH_FDR_lt_0.05"] = core["BH_FDR"].lt(0.05)
    core = core.sort_values("absolute_log2FC", ascending=False).reset_index(drop=True)
    core["absolute_change_rank"] = np.arange(1, len(core) + 1)
    plotted = core.head(TOP_N).copy()
    return core, plotted


def draw_chart(plotted: pd.DataFrame) -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 12,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

    fig, (ax, ax_far) = plt.subplots(
        1,
        2,
        figsize=(8.6, 6.3),
        sharey=True,
        gridspec_kw={"width_ratios": [6.2, 1.55], "wspace": 0.035},
    )
    fig.patch.set_facecolor("white")

    y = np.arange(len(plotted))
    values = plotted["resistant_vs_sensitive_log2FC"].to_numpy()
    colors = np.where(values >= 0, UP, DOWN)
    labels = [
        f"{gene}  {significance}"
        for gene, significance in zip(plotted["gene"], plotted["significance"])
    ]

    for axis in (ax, ax_far):
        axis.barh(y, values, height=0.58, color=colors, edgecolor="white", linewidth=1.0)
        axis.grid(axis="x", color=GRID, linewidth=0.8)
        axis.set_axisbelow(True)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.spines["left"].set_visible(False)
        axis.tick_params(axis="x", labelsize=15, colors=MUTED)
        axis.tick_params(axis="y", length=0, pad=10, colors=INK, labelsize=12)

    ax.axvline(0, color=INK, linewidth=1.0, zorder=0)
    ax.set_xlim(-0.60, 0.68)
    ax.set_xticks([-0.5, 0.0, 0.5])
    ax.set_yticks(y, labels=labels, fontsize=12, fontstyle="italic")
    ax.invert_yaxis()

    ax_far.set_xlim(8.80, 9.55)
    ax_far.set_xticks([9.0, 9.5])
    ax_far.tick_params(axis="y", labelleft=False)

    d = 0.012
    left_kwargs = dict(transform=ax.transAxes, color=INK, clip_on=False, linewidth=1.15)
    ax.plot((1 - d, 1 + d), (-d, +d), **left_kwargs)
    ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **left_kwargs)
    right_kwargs = dict(transform=ax_far.transAxes, color=INK, clip_on=False, linewidth=1.15)
    ax_far.plot((-d, +d), (-d, +d), **right_kwargs)
    ax_far.plot((-d, +d), (1 - d, 1 + d), **right_kwargs)

    for i, (gene, value) in enumerate(zip(plotted["gene"], values)):
        if gene == "RAD51C":
            ax_far.text(value + 0.035, i, f"{value:+.2f}", ha="left", va="center", fontsize=12, fontweight="bold", color=INK)
        elif value >= 0:
            ax.text(value + 0.022, i, f"{value:+.2f}", ha="left", va="center", fontsize=12, color=INK)
        else:
            ax.text(value - 0.022, i, f"{value:+.2f}", ha="right", va="center", fontsize=12, color=INK)

    fig.suptitle(
        "Core HR-gene expression changes during acquired PARPi resistance",
        x=0.12,
        y=0.975,
        ha="left",
        va="top",
        fontsize=16,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.12,
        0.925,
        "GSE163854 PH039 HGSOC PDX · DESeq2 · 5 sensitive vs 5 niraparib-resistant tumors",
        ha="left",
        va="top",
        fontsize=12,
        color=MUTED,
    )

    legend = [
        Patch(facecolor=UP, edgecolor="none", label="Higher in resistant"),
        Patch(facecolor=DOWN, edgecolor="none", label="Lower in resistant"),
    ]
    fig.legend(
        handles=legend,
        loc="upper left",
        bbox_to_anchor=(0.12, 0.865),
        frameon=False,
        ncol=2,
        fontsize=12,
        handlelength=1.3,
        columnspacing=1.7,
    )

    fig.text(
        0.12,
        0.050,
        "ns, BH FDR ≥ 0.05; * < 0.05; ** < 0.01; *** < 0.001. DESeq2 Wald test; genome-wide BH correction.",
        ha="left",
        va="bottom",
        fontsize=10.5,
        color=MUTED,
    )

    fig.subplots_adjust(left=0.22, right=0.92, top=0.79, bottom=0.15)
    fig.savefig(OUT_STEM.with_suffix(".png"), dpi=600, facecolor="white", bbox_inches="tight")
    fig.savefig(OUT_STEM.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(OUT_STEM.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT_STEM.parent.mkdir(parents=True, exist_ok=True)
    core, plotted = run_deseq2()
    processed = ROOT / "data" / "processed"
    processed.mkdir(parents=True, exist_ok=True)
    core.to_csv(processed / "core_HR_DESeq2_all15_statistics.csv", index=False)
    plotted.to_csv(processed / "core_HR_DESeq2_top8_plotted.csv", index=False)
    draw_chart(plotted)


if __name__ == "__main__":
    main()
