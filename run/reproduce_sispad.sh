#!/usr/bin/env bash
# Reproduce the SISPAD 2026 abstract figures (SISPAD_2026_Abhineet-9.pdf) from
# scratch: rank-1 projected NEGF-SCBA sweep + figure regeneration.
#
# Authoritative device/method (see docs/REPRODUCE_SISPAD.md):
#   ZnO/Mg0.3Zn0.7O symmetric RTD (2/3/2 nm), 10 um x 10 um pixel,
#   0-800 mV, Anderson mixing (DIIS, M=8) beta=0.3, max_iter=100, tol=1e-4,
#   Ef=20 meV, a=0.2 nm (Np=135), dE=2 meV, D2_bulk=0.001 eV2, mol D0^2=0.1 eV2.
#   Expected: Baseline peak |I_R| ~ 69 nA at 560 mV (300 K); Mol_A 10 K ~ 110 nA.
#
# Grid: the published Fig 3 figures used the 51-bias-point (16 mV) sweep — this
# is the clean grid whose I-V and raw-second-difference d2I/dV2 match the paper
# (the 201-pt grid is denser but its raw d2I is jitter-dominated). Override the
# point count with the first arg, e.g. `bash run/reproduce_sispad.sh 201`.
#
# Usage:
#   bash run/reproduce_sispad.sh [V_POINTS] [OUT_DIR]
# Defaults: V_POINTS=51, OUT_DIR=results/reproduce_<today>
#
# A 51-pt run is ~30 min/molecule (4 species ~= 2 h). For a fast confidence
# check that the code still reproduces the numbers without any sweep, run
# instead:  python run/verify_reproduction.py
set -euo pipefail
cd "$(dirname "$0")/.."
# shellcheck disable=SC1091
source venv/bin/activate 2>/dev/null || true

VPTS="${1:-51}"
OUTDIR="${2:-results/reproduce_$(date +%F)}"
DEVICE="ZnO_MgZnO_symmetric"
COMMON=(--device "$DEVICE" --V-min 0 --V-max 0.8 --V-points "$VPTS"
        --scba-mix 0.3 --scba-max-iter 100 --scba-tol 1e-4 --out-dir "$OUTDIR")

mkdir -p "$OUTDIR"
echo "==> Reproducing SISPAD -9 into $OUTDIR  (V_points=$VPTS)"

# Fig 3 (300 K): Baseline + 3 analytes.
python run/run_rank1_sweep.py "${COMMON[@]}" --T 300 \
    --molecules Baseline Mol_A Mol_B Mol_AB

# Fig 4 (temperature sweep): Mol_A at 10/77/150 K (300 K already done above).
for T in 10 77 150; do
    python run/run_rank1_sweep.py "${COMMON[@]}" --T "$T" --molecules Mol_A
done

# Paper data figures (Fig 3a/b/c + ΔI): fig1_IV, fig2_d2IdV2, fig3_deltaI, fig4_deltaD.
# Recovered generator; d2I/dV2 = raw np.diff(I,2)/dV2, no smoothing.
python run/make_paper_figs.py --data-dir "$OUTDIR" --out-dir "$OUTDIR"

echo "==> Done. Figures + .npz in $OUTDIR"
echo "    Compare against golden: results/sispad_scba_2026-04-14/"
echo "    NOTE: paper Fig 4 (temperature, temp.png) is not regenerated here — pending."
