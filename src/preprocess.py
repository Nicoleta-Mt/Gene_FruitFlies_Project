"""Filter low-expression genes and log2-transform the replicate matrix."""

import numpy as np
import pandas as pd

from config import LOG_EXPRESSION_PATH, MIN_FPKM, MIN_SAMPLES, REPLICATE_MATRIX_PATH


def run():
    fpkm = pd.read_parquet(REPLICATE_MATRIX_PATH)

    is_kept = (fpkm >= MIN_FPKM).sum(axis="columns") >= MIN_SAMPLES
    log_expression = np.log2(fpkm[is_kept] + 1)

    log_expression.to_parquet(LOG_EXPRESSION_PATH)

    print(f"Genes kept: {is_kept.sum():,}")
    print(f"Genes removed: {(~is_kept).sum():,}")
    print(f"Samples: {log_expression.shape[1]}")


if __name__ == "__main__":
    run()