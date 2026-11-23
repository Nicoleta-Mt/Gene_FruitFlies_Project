# Fruit Fly Gene Co-expression Networks

Compares gene co-expression networks across life stages, tissues, and sexes in
*Drosophila melanogaster*, using [FlyAtlas2](https://motif.mvls.gla.ac.uk/FlyAtlas2/) data.

## Setup

Requires [uv](https://docs.astral.sh/uv/) and Docker.

```bash
uv sync
mkdir -p data/raw
curl -L -o data/raw/flyatlas2-2025-03-03.sql \
  https://motif.mvls.gla.ac.uk/downloads/FlyAtlas2_2025.03.03.sql
docker compose up -d --wait
```

## Run

```bash
uv run src/main.py
```

Settings are in `src/config.py`. Output goes to `data/processed` and `results/figures`.