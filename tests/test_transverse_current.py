"""Tests for the transverse-mode (Tsu-Esaki) device current.

``run_rank1_keldysh_single_bias`` returns the current of a single spin in a
single transverse mode (prefactor q/h). ``tsu_esaki_current`` performs the
transverse-mode sum analytically for a sensor pixel, which is the treatment
documented in ``docs/STACK_DECISION.md`` Sec. 3 and
``docs/METHOD_DERIVATION.md`` Sec. 9.4.
"""
from __future__ import annotations

import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from core.scba_rank1_keldysh import tsu_esaki_current

_QE = 1.602176634e-19
_HBAR = 1.054571817e-34
_M_EFF = 0.28 * 9.10938356e-31
_KT = 0.02585
_AREA = 1e-10  # (10 um)^2


def _transmission(E):
    """A resonance on a small background."""
    return 0.9 * np.exp(-0.5 * ((E - 0.296) / 0.004) ** 2) + 0.02


class TestTsuEsakiCurrent(unittest.TestCase):
    def test_matches_explicit_transverse_integration(self) -> None:
        """The analytic supply function must equal the explicit E_t integral.

        That reduction is the entire content of the formula, so it is the one
        thing worth checking numerically rather than by inspection.
        """
        E = np.linspace(-0.05, 0.60, 4001)
        T = _transmission(E)
        mu_L, mu_R = 0.30, -0.26
        analytic = tsu_esaki_current(E, T, mu_L, mu_R, _KT, _M_EFF, _AREA)

        f = lambda x: 1.0 / (1.0 + np.exp(np.clip(x / _KT, -600, 600)))
        E_t = np.linspace(0.0, 2.0, 20001)
        trapz = getattr(np, "trapezoid", None) or np.trapz
        inner = np.array([trapz(f(e + E_t - mu_L) - f(e + E_t - mu_R), E_t) for e in E])
        # A * (2q/h) * (m/(2 pi hbar^2)) * int dE int dE_t ... ; q^2 converts both eV->J
        pre = (_AREA * (2 * _QE / (2 * np.pi * _HBAR))
               * (_M_EFF / (2 * np.pi * _HBAR ** 2)) * _QE * _QE)
        numeric = pre * trapz(T * inner, E)
        self.assertLess(abs(analytic - numeric) / abs(numeric), 1e-5,
                        f"analytic={analytic:.6e} numeric={numeric:.6e}")

    def test_zero_at_zero_bias(self) -> None:
        E = np.linspace(-0.05, 0.60, 2001)
        I = tsu_esaki_current(E, _transmission(E), 0.02, 0.02, _KT, _M_EFF, _AREA)
        self.assertLess(abs(I), 1e-30, f"I={I}")

    def test_linear_in_area(self) -> None:
        E = np.linspace(-0.05, 0.60, 2001)
        T = _transmission(E)
        a = tsu_esaki_current(E, T, 0.30, -0.26, _KT, _M_EFF, _AREA)
        b = tsu_esaki_current(E, T, 0.30, -0.26, _KT, _M_EFF, 4.0 * _AREA)
        self.assertAlmostEqual(b / a, 4.0, places=10)

    def test_antisymmetric_in_bias(self) -> None:
        E = np.linspace(-0.05, 0.60, 2001)
        T = _transmission(E)
        fwd = tsu_esaki_current(E, T, 0.30, -0.26, _KT, _M_EFF, _AREA)
        rev = tsu_esaki_current(E, T, -0.26, 0.30, _KT, _M_EFF, _AREA)
        self.assertAlmostEqual(fwd / rev, -1.0, places=10)


if __name__ == "__main__":
    unittest.main()
