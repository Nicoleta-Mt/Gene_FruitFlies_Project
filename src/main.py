"""Run the full pipeline from the database to the figures."""

import extract
import figures
import network
import preprocess

STEPS = [extract, preprocess, network, figures]

if __name__ == "__main__":
    for step in STEPS:
        print(f"\n--- {step.__name__} ---")
        step.run()