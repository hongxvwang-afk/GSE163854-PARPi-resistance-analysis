#!/usr/bin/env python3
"""Normalize GSE163854 counts and generate the row-Z-score heatmap matrix."""

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
COUNTS = ROOT / "data" / "raw" / "GSE163854_counts.txt.gz"
OUT_DIR = ROOT / "data" / "processed"

SAMPLES = ["P2A1", "P2A2", "P2A3", "F1A3", "F1A4", "F351", "F352", "F451", "F452", "F453"]

MODULES = {
    "BLM–BTRR / HR": ["BLM", "TOP3A", "RMI1", "RMI2", "BRCA1", "BRCA2", "PALB2", "BARD1", "BRIP1", "RAD51", "RAD51B", "RAD51C", "RAD51D", "XRCC2", "XRCC3"],
    "ATR checkpoint": ["ATR", "ATRIP", "CHEK1", "CLSPN", "TOPBP1", "TIMELESS", "TIPIN", "WEE1", "CDC25A", "RAD17"],
    "Fork protection / FA": ["FANCD2", "FANCI", "FANCM", "SMARCAL1", "ZRANB3", "HLTF", "WRN", "RECQL", "PARP1"],
    "End resection": ["MRE11A", "RAD50", "NBN", "RBBP8", "EXO1", "DNA2"],
    "R-loop processing": ["DHX9", "RNASEH1", "RNASEH2A", "SETX", "DDX5", "DDX21", "AQR"],
    "Cell cycle / RS": ["CCNE1", "CDK2", "CDC6", "MCM2", "MCM4", "MCM6", "CDC45", "PCNA", "RRM1", "RRM2"],
    "NHEJ / alt-EJ": ["TP53BP1", "RIF1", "MAD2L2", "XRCC5", "XRCC6", "LIG4", "POLQ", "PRKDC"],
    "PI3K–AKT": ["AKT1", "PIK3CA", "PIK3R1", "MTOR", "PTEN"],
}


def median_ratio_normalize(counts: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    positive_logs = np.log(counts.where(counts > 0))
    geometric_means = np.exp(positive_logs.mean(axis=1, skipna=True))
    valid = geometric_means.gt(0) & np.isfinite(geometric_means)
    size_factors = {}
    for sample in counts.columns:
        ratios = counts.loc[valid, sample] / geometric_means.loc[valid]
        ratios = ratios.replace([0, np.inf, -np.inf], np.nan).dropna()
        size_factors[sample] = ratios.median()
    size_factors = pd.Series(size_factors, dtype=float)
    size_factors /= np.exp(np.log(size_factors).mean())
    return counts.div(size_factors, axis=1), size_factors


def main() -> None:
    raw = pd.read_csv(COUNTS, sep="\t", index_col=0)
    if raw.index.duplicated().any():
        raw = raw.groupby(level=0).sum()
    missing_samples = [sample for sample in SAMPLES if sample not in raw.columns]
    if missing_samples:
        raise ValueError(f"Missing samples: {missing_samples}")

    genes = [gene for module_genes in MODULES.values() for gene in module_genes]
    missing_genes = [gene for gene in genes if gene not in raw.index]
    if missing_genes:
        raise ValueError(f"Missing genes: {missing_genes}")

    normalized, size_factors = median_ratio_normalize(raw.astype(float))
    log_expression = np.log2(normalized.loc[genes, SAMPLES] + 1)
    row_mean = log_expression.mean(axis=1)
    row_sd = log_expression.std(axis=1, ddof=1).replace(0, np.nan)
    z = log_expression.sub(row_mean, axis=0).div(row_sd, axis=0).fillna(0)
    z.insert(0, "module", [module for module, module_genes in MODULES.items() for _ in module_genes])

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    z.to_csv(OUT_DIR / "heatmap_row_Z_scores.csv")
    log_expression.to_csv(OUT_DIR / "heatmap_log2_normalized_expression.csv")
    size_factors.rename("size_factor").to_csv(OUT_DIR / "heatmap_size_factors.csv")


if __name__ == "__main__":
    main()
