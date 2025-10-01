#!/usr/bin/env bash
set -euo pipefail
conda activate rnaseq-pipeline
python -m src.main all --config-path config/config.yaml --samples-csv samples/samples.csv
