"""Build the replicate matrix (genes x samples) from the FlyAtlas2 database."""

import pandas as pd
from sqlalchemy import create_engine

from config import DB_URL, PROCESSED_DIR, REPLICATE_MATRIX_PATH

TISSUE_COLUMNS = {
    "TissueID": "tissue_id",
    "Stage": "stage",
    "Sex": "sex",
    "TissueName": "tissue_name",
}
EXPRESSION_COLUMNS = {
    "FBgn": "gene_id",
    "TissueID": "tissue_id",
    "Replicate1": "r1",
    "Replicate2": "r2",
    "Replicate3": "r3",
}


def read_table(engine, name, columns):
    return pd.read_sql_table(name, engine, columns=list(columns)).rename(
        columns=columns
    )


def run():
    engine = create_engine(DB_URL)
    tissues = read_table(engine, "Tissue", TISSUE_COLUMNS).set_index("tissue_id")
    expression = read_table(engine, "GeneFPKM", EXPRESSION_COLUMNS)

    conditions = (
        tissues["stage"]
        .str.cat(tissues[["sex", "tissue_name"]], sep="_")
        .str.lower()
        .str.replace(" ", "_")
    )

    replicate_matrix = (
        expression.assign(condition=expression["tissue_id"].map(conditions))
        .drop(columns="tissue_id")
        .melt(
            id_vars=["gene_id", "condition"],
            var_name="replicate",
            value_name="fpkm",
        )
        .dropna(subset="fpkm")
        .assign(sample=lambda df: df["condition"] + "_" + df["replicate"])
        .pivot(index="gene_id", columns="sample", values="fpkm")
    )

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    replicate_matrix.to_parquet(REPLICATE_MATRIX_PATH)

    print(f"Genes: {replicate_matrix.shape[0]:,}")
    print(f"Samples: {replicate_matrix.shape[1]}")
    print(f"Missing values: {replicate_matrix.isna().sum().sum()}")


if __name__ == "__main__":
    run()