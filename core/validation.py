"""Numerical diagnostics. These are not statistical or experimental errors."""
from __future__ import annotations

import numpy as np


def current_diagnostics(I_left, I_right, converged, rtol=2e-3, atol=1e-14):
    """Check both solver status and Kirchhoff balance, including at zero bias.

    Both contact currents use the solver's outward sign convention, so their
    sum must vanish. Absolute tolerance handles equilibrium/vanishing current.
    Passing this check does not establish energy/bias-grid convergence.
    """
    left, right = np.broadcast_arrays(np.asarray(I_left), np.asarray(I_right))
    defect = np.abs(left + right)
    scale = np.maximum(np.abs(left), np.abs(right))
    relative = np.divide(defect, scale, out=np.zeros_like(defect, dtype=float), where=scale > 0)
    valid = (np.asarray(converged, dtype=bool) & np.isfinite(left) & np.isfinite(right)
             & (defect <= atol + rtol * scale))
    return dict(conservation_abs_A=defect, conservation_rel=relative, current_valid=valid)


def central_second_difference(values, step):
    """Second derivative at the centre of an equally spaced three-point stencil."""
    values = np.asarray(values)
    if values.shape[-1] != 3 or step <= 0:
        raise ValueError("Need three samples and a positive voltage step")
    return (values[..., 0] - 2 * values[..., 1] + values[..., 2]) / step**2


def refinement_agrees(coarse, fine, rtol, atol):
    """Elementwise refinement test with an explicit floor in observable units."""
    coarse, fine = np.asarray(coarse), np.asarray(fine)
    return (np.isfinite(coarse) & np.isfinite(fine)
            & (np.abs(coarse - fine) <= atol + rtol * np.abs(fine)))


def require_plot_data(data, field, allow_unvalidated=False):
    """Reject missing currents and uncertified spectra before making figures.

    Diagnostic/historical plots may explicitly bypass numerical certification,
    but cannot turn an uncomputed device current (old zero placeholder) into data.
    """
    if field not in data:
        raise ValueError(f"Missing current field {field}")
    if field == "I_device" and not bool(data.get("transverse", False)):
        raise ValueError("I_device was not calculated; rerun with --transverse")
    values = np.asarray(data[field])
    if values.shape != np.asarray(data["V"]).shape or not np.all(np.isfinite(values)):
        raise ValueError(f"{field} has invalid or uncomputed values")
    V = np.asarray(data["V"])
    if V.size < 3 or not np.all(np.diff(V) > 0) or not np.allclose(np.diff(V), V[1]-V[0]):
        raise ValueError("Second-difference plots require an increasing uniform bias grid")
    if allow_unvalidated:
        return
    diag = current_diagnostics(data["I_L"], data["I_R"], data["converged"])
    valid = diag["current_valid"]
    for key in ("numerical_valid", "poisson_converged", "physical_valid"):
        if key in data:
            valid = valid & np.asarray(data[key])
    if field == "I_device":
        valid = valid & np.asarray(data.get("transverse_valid", False))
    if not np.all(valid):
        raise ValueError("Solver, conservation, or Poisson checks failed")
    if not bool(data.get("spectrum_grid_validated", False)):
        raise ValueError("Spectrum lacks energy/bias-grid validation; use an explicit diagnostic plot option")
    if not np.all(data.get("physical_valid", False)):
        raise ValueError("SCBA rates/populations lack physical validation; rerun with the current solver")
