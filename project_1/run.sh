#!/usr/bin/env bash
# Project‐local run script: reproducible, no global tweaks needed.

# Clean old Snakemake state
rm -rf .snakemake

# Invoke Snakemake with:
#  - use conda
#  - conda frontend (plain conda)
#  - a project‐local prefix for envs
#  - a fixed core count
snakemake \
  --use-conda \
  --conda-frontend conda \
  --conda-prefix $(pwd)/.conda_envs \
  --cores 4 \
  "$@"
