#!/usr/bin/env bash
# Re-download the AutoDock Vina Windows binaries (gitignored — not committed).
# Run from the repo root:  bash vendor/vina_bin/setup.sh
set -e
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VER="1.2.7"
BASE="https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v${VER}"
echo "Downloading AutoDock Vina ${VER} Windows binaries -> $DIR"
curl -sL -o "$DIR/vina.exe"       "$BASE/vina_${VER}_win.exe"
curl -sL -o "$DIR/vina_split.exe" "$BASE/vina_split_${VER}_win.exe"
"$DIR/vina.exe" --version && echo "OK"

# The docking chem stack lives in a separate micromamba env; create it once with:
#   MAMBA_ROOT_PREFIX=.qc .qc/Library/bin/micromamba.exe create -y -n drugdisc-agent \
#     -f vendor/AtomisticSkills/conda-envs/drugdisc-agent/core_env.yaml
# On Windows, drop the `vina` and `fpocket` conda deps (no win-64 build) — this
# binary replaces the Vina Python API. See services/vina_engine.py.
