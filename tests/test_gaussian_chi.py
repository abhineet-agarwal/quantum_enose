"""
Unit tests for the Gaussian longitudinal form factor χ(z).

This replaces the Fix 6 `mol_radius = barrier_half_width` top-hat heuristic
with the physical σ_mol = 3 Å Gaussian derived from METHOD_DERIVATION.md §2
Eq. (1). The rank-1 projected SCBA uses ``gaussian_chi_projector`` upstream
of ``projected_dyson_step`` to build the spatial envelope of σ̃^R; these
tests cover the sampling, the three normalization modes, and the behaviour
in two physically meaningful limits.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from core.self_energy import gaussian_chi_projector


# --- fixtures ----------------------------------------------------------------


@pytest.fixture
def sispad_grid():
    """Typical SISPAD 2026 grid: 1.5 Å spacing, 80 sites, molecule centred."""
    return dict(Np=80, mol_site=40, sigma_mol_m=3e-10, dz_m=1.5e-10)


# --- shape / value checks ----------------------------------------------------


def test_chi_shape_and_peak(sispad_grid):
    chi = gaussian_chi_projector(**sispad_grid)
    assert chi.shape == (sispad_grid["Np"],)
    # Peak at the molecule site.
    assert np.argmax(chi) == sispad_grid["mol_site"]
    # Default normalization is peak → 1.
    assert np.isclose(chi.max(), 1.0)


def test_chi_gaussian_profile(sispad_grid):
    """χ(z) = exp(-(z-z_mol)^2 / 2σ²) up to normalization."""
    chi = gaussian_chi_projector(**sispad_grid)
    sigma = sispad_grid["sigma_mol_m"]
    dz = sispad_grid["dz_m"]
    z0 = sispad_grid["mol_site"]
    i = np.arange(sispad_grid["Np"])
    z = (i - z0) * dz
    expected = np.exp(-0.5 * (z / sigma) ** 2)
    expected /= expected.max()  # match peak normalization
    np.testing.assert_allclose(chi, expected, atol=1e-14, rtol=0)


def test_chi_symmetry(sispad_grid):
    chi = gaussian_chi_projector(**sispad_grid)
    z0 = sispad_grid["mol_site"]
    Np = sispad_grid["Np"]
    # Symmetric about the centre over the common half-width.
    half = min(z0, Np - 1 - z0)
    left = chi[z0 - half:z0][::-1]
    right = chi[z0 + 1:z0 + 1 + half]
    np.testing.assert_allclose(left, right, atol=1e-14, rtol=0)


# --- normalization modes -----------------------------------------------------


def test_peak_normalization(sispad_grid):
    chi = gaussian_chi_projector(**sispad_grid, normalization="peak")
    assert np.isclose(chi.max(), 1.0, atol=1e-14)


def test_l1_normalization(sispad_grid):
    chi = gaussian_chi_projector(**sispad_grid, normalization="l1")
    # Σ χ_i · dz ≈ 1 (Riemann sum of the continuum integral).
    riemann = chi.sum() * sispad_grid["dz_m"]
    assert np.isclose(riemann, 1.0, atol=1e-12)


def test_l2_normalization(sispad_grid):
    chi = gaussian_chi_projector(**sispad_grid, normalization="l2")
    riemann = (chi ** 2).sum() * sispad_grid["dz_m"]
    assert np.isclose(riemann, 1.0, atol=1e-12)


def test_l1_continuum_limit():
    """A well-resolved Gaussian should match the analytic ∫ χ dz = 1 target."""
    sigma = 3e-10
    dz = 0.1e-10                  # 20× oversampled
    Np = 601
    chi = gaussian_chi_projector(Np=Np, mol_site=Np // 2,
                                 sigma_mol_m=sigma, dz_m=dz,
                                 normalization="l1")
    # The analytic continuum peak value for a normalized Gaussian is
    #   χ_peak = 1 / (sqrt(2π) · σ)
    expected_peak = 1.0 / (np.sqrt(2 * np.pi) * sigma)
    assert np.isclose(chi.max(), expected_peak, rtol=1e-4)


# --- physical-limit sanity checks -------------------------------------------


def test_narrow_gaussian_approaches_delta(sispad_grid):
    """σ_mol → 0 should localize χ to a single site at the molecule centre."""
    narrow = {**sispad_grid, "sigma_mol_m": 1e-12}  # ≪ dz
    chi = gaussian_chi_projector(**narrow)
    # Peak site is ~1, neighbours are astronomically small.
    assert chi[sispad_grid["mol_site"]] == pytest.approx(1.0)
    neighbours = np.concatenate([chi[:sispad_grid["mol_site"]],
                                 chi[sispad_grid["mol_site"] + 1:]])
    assert np.max(neighbours) < 1e-30


def test_wide_gaussian_spans_grid(sispad_grid):
    """σ_mol ≫ grid extent → χ nearly flat (continuum limit of top-hat)."""
    wide = {**sispad_grid, "sigma_mol_m": 1e-7}  # 1 µm, grid is 12 nm
    chi = gaussian_chi_projector(**wide)
    # Every site is within << σ of the centre, so χ ≈ 1 everywhere.
    assert np.all(chi > 0.997)


# --- sandwich application ----------------------------------------------------


def test_sandwich_against_diag_matrix(sispad_grid):
    """chi[:, None] * Σ * chi[None, :] == diag(chi) @ Σ @ diag(chi)."""
    rng = np.random.default_rng(20260408)
    Np = sispad_grid["Np"]
    Sigma = (rng.standard_normal((Np, Np))
             + 1j * rng.standard_normal((Np, Np)))
    chi = gaussian_chi_projector(**sispad_grid)

    sandwich_outer = chi[:, None] * Sigma * chi[None, :]
    sandwich_diag = np.diag(chi) @ Sigma @ np.diag(chi)
    np.testing.assert_allclose(sandwich_outer, sandwich_diag,
                               atol=1e-14, rtol=0)


# --- input validation --------------------------------------------------------


@pytest.mark.parametrize("bad_sigma", [0, -1e-10])
def test_rejects_nonpositive_sigma(bad_sigma):
    with pytest.raises(ValueError):
        gaussian_chi_projector(Np=10, mol_site=5,
                               sigma_mol_m=bad_sigma, dz_m=1e-10)


@pytest.mark.parametrize("bad_dz", [0, -1e-10])
def test_rejects_nonpositive_dz(bad_dz):
    with pytest.raises(ValueError):
        gaussian_chi_projector(Np=10, mol_site=5,
                               sigma_mol_m=3e-10, dz_m=bad_dz)


@pytest.mark.parametrize("bad_site", [-1, 10, 11])
def test_rejects_out_of_range_mol_site(bad_site):
    with pytest.raises(ValueError):
        gaussian_chi_projector(Np=10, mol_site=bad_site,
                               sigma_mol_m=3e-10, dz_m=1e-10)


def test_rejects_unknown_normalization():
    with pytest.raises(ValueError):
        gaussian_chi_projector(Np=10, mol_site=5, sigma_mol_m=3e-10,
                               dz_m=1e-10, normalization="l42")
