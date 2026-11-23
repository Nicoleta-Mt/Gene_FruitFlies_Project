"""Build the co-expression network from log-transformed expression."""

import numpy as np
import pandas as pd

from config import LOG_EXPRESSION_PATH, MIN_CORRELATION, NETWORK_EDGES_PATH


def run():
    log_expression = pd.read_parquet(LOG_EXPRESSION_PATH)
    genes = log_expression.index.to_numpy()

    correlation = np.corrcoef(log_expression, dtype=np.float32)
    is_edge = np.triu(np.abs(correlation) >= MIN_CORRELATION, k=1)
    gene_a, gene_b = np.nonzero(is_edge)

    edges = pd.DataFrame(
        {
            "gene_a": genes[gene_a],
            "gene_b": genes[gene_b],
            "correlation": correlation[gene_a, gene_b],
        }
    )
    edges.to_parquet(NETWORK_EDGES_PATH, index=False)

    connected_genes = pd.concat([edges["gene_a"], edges["gene_b"]]).nunique()

    print(f"Connections: {len(edges):,}")
    print(f"Positive: {(edges['correlation'] > 0).sum():,}")
    print(f"Negative: {(edges['correlation'] < 0).sum():,}")
    print(f"Connected genes: {connected_genes:,}")
    print(f"Isolated genes: {len(genes) - connected_genes:,}")


if __name__ == "__main__":
    run()