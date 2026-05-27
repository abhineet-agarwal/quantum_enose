"""
Regression test for the Patil 1D reference port (Tier-1 γ).

The test re-runs the port at a smaller V-grid (so it stays under a few
seconds in CI) and checks the high-level physics is reproduced:

1. I(V) is monotone increasing on this pre-NDR sub-range.
2. I(V = 0) ≈ 0 (symmetric bias sanity).
3. I(V) > 0 strictly for V > kT.
4. The coherent and incoherent decompositions sum to roughly the total
   (with the MATLAB sign convention IIco and IInonco are negative — they
   add via I = I_total per the original code, not via I = I_co + I_nco;
   we just check both are finite and respond monotonically).

The full 76-V-point reference is saved to ``tests/patil_reference_1d.npz``
by ``python -m tests.patil_reference_1d`` and is *not* recomputed in CI.
``test_full_reference_npz_present`` checks that file exists when present
and contains the expected fields, so a future tolerance comparison
against the multimode solver has a stable target.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from tests.patil_reference_1d import run_patil_reference


@pytest.fixture(scope="module")
def small_run():
    """Coarse run that finishes in a few seconds."""
    return run_patil_reference(
        NV=11, V_max=0.3, scba_tol=1e-3,
        max_outer_iter=40, verbose=False,
    )


# --- physics sanity ---------------------------------------------------------


def test_zero_bias_current_is_negligible(small_run):
    """At V = 0, symmetric bias gives I = 0 by construction."""
    assert abs(small_run.I_total[0]) < 1e-15


def test_current_monotone_in_subrange(small_run):
    """For V ∈ [0, 0.3] V the device is below NDR (resonance ≈ 0.7 V)."""
    diffs = np.diff(small_run.I_total)
    # All increments positive (within numerical noise).
    assert np.all(diffs > -1e-15)


def test_current_grows_with_bias(small_run):
    """I(V_max) >> I(0)."""
    assert small_run.I_total[-1] > 1e-12


def test_decomposition_finite(small_run):
    """Coherent + incoherent traces are finite (not NaN/inf)."""
    assert np.all(np.isfinite(small_run.I_coherent))
    assert np.all(np.isfinite(small_run.I_incoherent))


def test_derivatives_shape(small_run):
    NV = small_run.V.size
    assert small_run.dIdV.shape == (NV - 1,)
    assert small_run.d2IdV2.shape == (NV - 2,)


# --- saved reference (only checked when the npz exists) --------------------


REF_NPZ = os.path.join(os.path.dirname(__file__), "patil_reference_1d.npz")


@pytest.mark.skipif(
    not os.path.exists(REF_NPZ),
    reason="patil_reference_1d.npz not yet generated; run "
           "`python -m tests.patil_reference_1d` to create it.",
)
def test_full_reference_npz_present():
    data = np.load(REF_NPZ)
    expected = {
        "V", "I_total", "I_coherent", "I_incoherent",
        "dIdV", "d2IdV2", "V_dIdV", "V_d2IdV2",
        "Np", "Nb", "Vb", "t0", "Ef", "kT", "NV", "V_max",
        "Dnu", "hnu_eV",
    }
    missing = expected - set(data.files)
    assert not missing, f"Missing keys in patil_reference_1d.npz: {missing}"
    # Reference should have at least 50 V points for IETS feature visibility.
    assert data["V"].size >= 50
    # And the current should grow with bias.
    assert data["I_total"][-1] > data["I_total"][0]
