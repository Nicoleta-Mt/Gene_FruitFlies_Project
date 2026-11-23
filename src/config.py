"""Shared paths, database settings, and analysis parameters."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

PROCESSED_DIR = ROOT / "data" / "processed"
FIGURES_DIR = ROOT / "results" / "figures"

REPLICATE_MATRIX_PATH = PROCESSED_DIR / "replicate-matrix.parquet"
LOG_EXPRESSION_PATH = PROCESSED_DIR / "log-expression.parquet"
NETWORK_EDGES_PATH = PROCESSED_DIR / "network-edges.parquet"

DB_URL = "mysql+pymysql://fly:atlas@127.0.0.1:3307/FlyAtlas2"

# A gene counts as expressed in a sample at this FPKM or higher
MIN_FPKM = 1
# A gene is kept if it is expressed in at least this many samples
MIN_SAMPLES = 12
# Two genes are connected if their absolute correlation is at least this
MIN_CORRELATION = 0.8