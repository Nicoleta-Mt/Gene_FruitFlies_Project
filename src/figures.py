"""Create the figures for the first co-expression network."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import (
    FIGURES_DIR,
    LOG_EXPRESSION_PATH,
    MIN_CORRELATION,
    MIN_FPKM,
    REPLICATE_MATRIX_PATH,
)

BLUE = "#2a78d6"
GRAY = "#52514e"
LIGHT_GRAY = "#c9c8c2"
TISSUE_PATTERN = r"^(?:adult|larval)_(?:male|female|both)_(.+)_r\d$"
STYLE = {
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.bbox": "tight",
    "savefig.dpi": 200,
}


def save(fig, name):
    fig.savefig(FIGURES_DIR / f"{name}.png")
    plt.close(fig)


def correlate(expression):
    return np.corrcoef(expression, dtype=np.float32)


def plot_expression(fpkm):
    values = fpkm.to_numpy()
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.hist(np.log2(values[values > 0] + 1), bins=60, color=BLUE)
    ax.axvline(np.log2(MIN_FPKM + 1), color=GRAY, linestyle="--")
    ax.set(
        title=f"Detected expression, with the {MIN_FPKM} FPKM threshold",
        xlabel="log2(FPKM + 1)",
        ylabel="Measurements",
    )
    save(fig, "expression-distribution")


def plot_sample_similarity(log_expression):
    conditions = log_expression.columns.str.replace(r"_r\d$", "", regex=True)
    replicates = conditions.value_counts(sort=False)
    centers = replicates.cumsum() - replicates / 2 - 0.5

    fig, ax = plt.subplots(figsize=(11, 10))
    image = ax.imshow(correlate(log_expression.T), cmap="Blues", vmin=0.4, vmax=1)
    ax.set_xticks(centers, replicates.index, rotation=90, fontsize=7)
    ax.set_yticks(centers, replicates.index, fontsize=7)
    ax.set_title("Replicates of the same condition are most similar")
    fig.colorbar(image, label="Pearson correlation", shrink=0.8, extend="min")
    save(fig, "sample-correlation")


def plot_gene_correlations(correlation):
    pairs = correlation[np.triu_indices_from(correlation, k=1)]
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.hist(pairs, bins=100, range=(-1, 1), color=BLUE)
    for cutoff in (-MIN_CORRELATION, MIN_CORRELATION):
        ax.axvline(cutoff, color=GRAY, linestyle="--")
    ax.set(
        title=f"Only pairs beyond |r| = {MIN_CORRELATION} become connections",
        xlabel="Pearson correlation (r)",
        ylabel="Gene pairs",
    )
    save(fig, "correlation-distribution")


def plot_tissue_influence(log_expression, correlation):
    tissues = log_expression.columns.str.extract(TISSUE_PATTERN, expand=False)
    connected = np.abs(correlation) >= MIN_CORRELATION
    n_genes = len(connected)
    n_connections = np.count_nonzero(connected) - n_genes

    share_lost = {}
    for tissue in tissues.unique():
        without = correlate(log_expression.loc[:, tissues != tissue])
        n_kept = np.count_nonzero(connected & (np.abs(without) >= MIN_CORRELATION))
        share_lost[tissue] = 100 * (1 - (n_kept - n_genes) / n_connections)
    share_lost = pd.Series(share_lost).sort_values()

    top_tissue = share_lost.idxmax()
    colors = [BLUE if name == top_tissue else LIGHT_GRAY for name in share_lost.index]
    fig, ax = plt.subplots(figsize=(8, 7))
    bars = ax.barh(share_lost.index.str.replace("_", " "), share_lost, color=colors)
    ax.bar_label(bars, fmt="%.1f%%", padding=3, fontsize=8)
    ax.set(
        title=f"Removing {top_tissue} breaks the most connections",
        xlabel="Connections lost when the tissue is removed (%)",
    )
    save(fig, "tissue-influence")


def run():
    fpkm = pd.read_parquet(REPLICATE_MATRIX_PATH)
    log_expression = pd.read_parquet(LOG_EXPRESSION_PATH)
    correlation = correlate(log_expression)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    with plt.rc_context(STYLE):
        plot_expression(fpkm)
        plot_sample_similarity(log_expression)
        plot_gene_correlations(correlation)
        plot_tissue_influence(log_expression, correlation)


if __name__ == "__main__":
    run()