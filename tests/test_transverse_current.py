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

from core.scba_rank1_keldysh import (landauer_current_1mode, transverse_integrated_current,
                                     tsu_esaki_current)

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


class TestControlVariate(unittest.TestCase):
    """The E_t quadrature with the elastic control subtracted.

    Geometry mimics the device at the NDR peak: V = 576 mV, emitter band edge
    at +V/2, Ef = 20 meV, a narrow resonance just above the edge. The E_t
    window is mu_L + 12 kT, as in the sweep.
    """
    V, EF = 0.576, 0.020

    def _setup(self, kT):
        E = np.linspace(-0.25, 0.50, 3001)
        edge = self.V / 2
        T_ctrl = np.where(E > edge, 0.9 * np.exp(-0.5 * ((E - 0.296) / 0.004) ** 2) + 0.02, 0.0)
        # stand-in for the SCBA result: the control, rescaled, plus a sideband
        T_true = 1.04 * T_ctrl + np.where(
            E > edge, 0.05 * np.exp(-0.5 * ((E - 0.330) / 0.010) ** 2), 0.0)
        mk = lambda T: (lambda ef: landauer_current_1mode(E, T, ef + self.V / 2, ef - self.V / 2, kT))
        exact = lambda T: tsu_esaki_current(E, T, self.EF + self.V / 2, self.EF - self.V / 2,
                                            kT, _M_EFF, _AREA)
        e_max = self.EF + self.V / 2 + 12 * kT
        return mk, exact, T_ctrl, T_true, e_max

    def test_landauer_integrates_to_tsu_esaki(self) -> None:
        """Convention check: prefactors, spin and E-rule must all line up."""
        kT = 0.02585
        mk, exact, T_ctrl, _, e_max = self._setup(kT)
        I, _, _ = transverse_integrated_current(mk(T_ctrl), self.EF, kT, _M_EFF, _AREA,
                                                e_max, n_nodes=400)
        ref = exact(T_ctrl)
        self.assertLess(abs(I - ref) / ref, 1e-4, f"quad={I:.6e} tsu_esaki={ref:.6e}")

    def test_control_fixes_low_temperature(self) -> None:
        """At 10 K the plain 8-node rule is badly off; with the control it is not."""
        kT = 0.02585 * 10 / 300
        mk, exact, T_ctrl, T_true, e_max = self._setup(kT)
        ref = exact(T_true)
        plain, _, _ = transverse_integrated_current(mk(T_true), self.EF, kT, _M_EFF, _AREA,
                                                    e_max, n_nodes=8)
        cv, _, _ = transverse_integrated_current(mk(T_true), self.EF, kT, _M_EFF, _AREA,
                                                 e_max, n_nodes=8, control=mk(T_ctrl),
                                                 control_exact=exact(T_ctrl))
        err_plain, err_cv = abs(plain - ref) / ref, abs(cv - ref) / ref
        self.assertGreater(err_plain, 0.05, f"setup no longer reproduces the problem: {err_plain:.3%}")
        self.assertLess(err_cv, 0.01, f"control variate error {err_cv:.3%} (plain {err_plain:.3%})")

    def test_control_requires_exact(self) -> None:
        with self.assertRaises(ValueError):
            transverse_integrated_current(lambda ef: 0.0, 0.02, 0.025, _M_EFF, _AREA, 0.3,
                                          n_nodes=4, control=lambda ef: 0.0)


if __name__ == "__main__":
    unittest.main()
