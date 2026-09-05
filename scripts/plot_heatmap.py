#!/usr/bin/env python3
"""Plot the final pathway-annotated GSE163854 heatmap.

The script consumes the processed row-Z-score matrix committed to the
repository. The companion preprocessing description in docs/METHODS.md
explains how the matrix was generated from the GEO raw count matrix.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
Z_MATRIX = ROOT / "data" / "processed" / "heatmap_row_Z_scores.csv"
OUT_DIR = ROOT / "figures"
OUT_STEM = OUT_DIR / "GSE163854_pathway_heatmap_final"

SENSITIVE = ["P2A1", "P2A2", "P2A3", "F1A3", "F1A4"]
RESISTANT = ["F351", "F352", "F451", "F452", "F453"]
SAMPLES = SENSITIVE + RESISTANT

MODULES = {
    "BLM–BTRR / HR": [
        "BLM", "TOP3A", "RMI1", "RMI2", "BRCA1", "BRCA2", "PALB2",
        "BARD1", "BRIP1", "RAD51", "RAD51B", "RAD51C", "RAD51D",
        "XRCC2", "XRCC3",
    ],
    "ATR checkpoint": [
        "ATR", "ATRIP", "CHEK1", "CLSPN", "TOPBP1", "TIMELESS", "TIPIN",
        "WEE1", "CDC25A", "RAD17",
    ],
    "Fork protection / FA": [
        "FANCD2", "FANCI", "FANCM", "SMARCAL1", "ZRANB3", "HLTF", "WRN",
        "RECQL", "PARP1",
    ],
    "End resection": ["MRE11A", "RAD50", "NBN", "RBBP8", "EXO1", "DNA2"],
    "R-loop processing": ["DHX9", "RNASEH1", "RNASEH2A", "SETX", "DDX5", "DDX21", "AQR"],
    "Cell cycle / RS": [
        "CCNE1", "CDK2", "CDC6", "MCM2", "MCM4", "MCM6", "CDC45", "PCNA",
        "RRM1", "RRM2",
    ],
    "NHEJ / alt-EJ": ["TP53BP1", "RIF1", "MAD2L2", "XRCC5", "XRCC6", "LIG4", "POLQ", "PRKDC"],
    "PI3K–AKT": ["AKT1", "PIK3CA", "PIK3R1", "MTOR", "PTEN"],
}

GRID_COLOR = "#DCE4E7"
INK = "#20282C"
MUTED = "#4E595E"
CMAP = LinearSegmentedColormap.from_list(
    "blue_white_orange", ["#2E6F93", "#F7F7F4", "#C65D32"], N=256
)
NORM = Normalize(vmin=-2.4, vmax=2.4, clip=True)


def main() -> None:
    genes = [gene for module_genes in MODULES.values() for gene in module_genes]
    annotated_z = pd.read_csv(Z_MATRIX, index_col=0)
    missing = [gene for gene in genes if gene not in annotated_z.index]
    if missing:
        raise ValueError(f"Missing genes from Z-score matrix: {missing}")
    z = annotated_z.loc[genes, SAMPLES].astype(float)

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 12,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

    n_rows = len(genes)
    fig = plt.figure(figsize=(17.2, 19.6), facecolor="white")
    ax = fig.add_axes([0.075, 0.075, 0.73, 0.845])
    ax.set_xlim(-4.25, len(SAMPLES) + 0.15)
    ax.set_ylim(n_rows + 2.35, -2.65)
    ax.axis("off")

    fig.text(
        0.055, 0.979,
        "GSE163854 transcriptional landscape during acquired PARPi resistance in HGSOC",
        fontsize=18, fontweight="bold", color="#1F272B", ha="left", va="top",
    )
    fig.text(
        0.055, 0.952,
        "PH039 HGSOC PDX; 5 sensitive and 5 niraparib-resistant tumors; mixed-response tumors excluded; row-standardized RNA-seq expression",
        fontsize=12, color=MUTED, ha="left", va="top",
    )

    ax.text(2.5, -1.55, "Sensitive tumors", ha="center", va="center", fontsize=15, color="#263034")
    ax.text(7.5, -1.55, "Niraparib-resistant tumors", ha="center", va="center", fontsize=15, color="#263034")

    for y, gene in enumerate(genes):
        for x, sample in enumerate(SAMPLES):
            ax.add_patch(
                Rectangle(
                    (x + 0.035, y + 0.07), 0.93, 0.86,
                    facecolor=CMAP(NORM(z.loc[gene, sample])),
                    edgecolor=GRID_COLOR, linewidth=0.55,
                )
            )

    ax.plot([5, 5], [-2.05, n_rows + 0.95], color="#31383C", lw=1.5, solid_capstyle="butt")

    emphasized = {"BLM", "RAD51C", "ATR", "CHEK1", "FANCD2", "CCNE1", "POLQ"}
    row_start = 0
    for module, module_genes in MODULES.items():
        row_end = row_start + len(module_genes)
        center = (row_start + row_end) / 2
        bracket_x = -1.62
        ax.plot([bracket_x, bracket_x], [row_start + 0.08, row_end - 0.08], color="#263034", lw=1.45)
        ax.plot([bracket_x, bracket_x + 0.20], [row_start + 0.08, row_start + 0.08], color="#263034", lw=1.45)
        ax.plot([bracket_x, bracket_x + 0.20], [row_end - 0.08, row_end - 0.08], color="#263034", lw=1.45)
        ax.text(-2.02, center, module, ha="right", va="center", fontsize=12,
                fontweight="bold", color="#30383C")
        if row_end < n_rows:
            ax.plot([0, len(SAMPLES)], [row_end, row_end], color="white", lw=1.4)
        row_start = row_end

    for y, gene in enumerate(genes):
        ax.text(
            -0.28, y + 0.51, gene, ha="right", va="center", fontsize=12,
            fontstyle="italic", fontweight="bold" if gene in emphasized else "normal",
            color=INK,
        )

    for x, sample in enumerate(SAMPLES):
        ax.text(x + 0.50, n_rows + 0.35, sample, ha="right", va="top",
                rotation=45, rotation_mode="anchor", fontsize=11.5, color=INK)

    cbar_ax = fig.add_axes([0.845, 0.44, 0.025, 0.34])
    scalar_map = mpl.cm.ScalarMappable(norm=NORM, cmap=CMAP)
    scalar_map.set_array([])
    cbar = fig.colorbar(scalar_map, cax=cbar_ax, orientation="vertical")
    cbar.set_ticks([-2, 0, 2])
    cbar.ax.tick_params(length=4, width=0.8, labelsize=15)
    cbar.outline.set_visible(False)
    cbar.set_label("Expression (row Z-score)", fontsize=12.5, labelpad=9)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_STEM.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(OUT_STEM.with_suffix(".pdf"), facecolor="white")
    fig.savefig(OUT_STEM.with_suffix(".svg"), facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
