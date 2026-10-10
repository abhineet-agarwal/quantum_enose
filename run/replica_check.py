"""Are the extra NDR features LO-phonon replicas? Energy-resolved check.

Elastic resonant tunneling: electrons enter from the emitter and leave at the
collector at the same energy (the resonance). A phonon-emission replica:
they enter ~hw_LO above the resonance and leave at the resonance energy.
At each bias (converged coupled up-branch potential) this compares the
emitter injection spectrum and the collector spectrum with the resonance
energies from the ballistic transmission.
"""
from __future__ import annotations
import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np
import run.sc_scatter_sweep as sw
from core.poisson import ballistic_transmission

bc = sw.bc
HW = 0.072
d = np.load(os.path.join(os.path.dirname(bc.OUT), "..", "2026-10-07", "sc_coupled_D20.015_up.npz"), allow_pickle=True)
ph = sw.phonons(0.015)


def spectra(V, U):
    E = sw.window(V, 0.25)
    r = sw.scba(V, U, E, ph, None, tol=1e-5, max_iter=400)
    Ec = E + 1e-12j
    sig = lambda u, ub: -bc.t0 * np.exp(1j * np.arccos(1 - (Ec - u - ub) / (2 * bc.t0)))
    gL = np.real(1j * (sig(U[0], bc.UB[0]) - np.conj(sig(U[0], bc.UB[0]))))
    gR = np.real(1j * (sig(U[-1], bc.UB[-1]) - np.conj(sig(U[-1], bc.UB[-1]))))
    fL = 1 / (1 + np.exp((E - bc.EF - V / 2) / bc.kT)); fR = 1 / (1 + np.exp((E - bc.EF + V / 2) / bc.kT))
    inj = -(gL * ((1 - fL) * r.n_diag[:, 0] - fL * r.p_diag[:, 0]))      # entering from emitter
    out = gR * ((1 - fR) * r.n_diag[:, -1] - fR * r.p_diag[:, -1])        # leaving into collector
    return E, inj, out, r


def resonances(V, U):
    E = np.arange(V / 2 - 0.25, bc.EF + V / 2 + 12 * bc.kT, 0.00005)
    T = ballistic_transmission(E, bc.H_z, bc.UB, U, bc.t0, method="rgf")
    pk = [k for k in range(1, len(T) - 1) if T[k] > T[k - 1] and T[k] >= T[k + 1] and T[k] > 1e-3]
    return [(E[k], T[k]) for k in pk]


def weighted(E, w, lo, hi):
    m = (E >= lo) & (E < hi)
    return float(w[m].sum() / max(w.sum(), 1e-300))


if __name__ == "__main__":
    for Vm in (224, 480, 512, 1216, 1408, 1440):
        V = Vm / 1e3; U = d[str(Vm)].item()["U"]
        E, inj, out, r = spectra(V, U)
        res = resonances(V, U)
        Ein = float((E * inj).sum() / inj.sum()); Eout = float((E * out).sum() / out.sum())
        print(f"V = {Vm} mV  (I {d[str(Vm)].item()['I_scba_mode']*1e9:.3f} nA, emitter band edge {V/2*1e3:.0f} meV, mu_L {(bc.EF+V/2)*1e3:.0f} meV)")
        print(f"   resonances (T > 1e-3): " + ", ".join(f"{e*1e3:.1f} meV (T {t:.2g})" for e, t in res))
        print(f"   current-weighted energy: enters {Ein*1e3:.1f} meV, leaves {Eout*1e3:.1f} meV, "
              f"drop {(Ein-Eout)*1e3:.1f} meV ({(Ein-Eout)/HW:.2f} hw_LO)")
        for e0, _ in res:
            print(f"   share of collected current within 4 meV of {e0*1e3:.1f} meV: {weighted(E, out, e0-0.004, e0+0.004):.0%}; "
                  f"of injected current within 4 meV of it: {weighted(E, inj, e0-0.004, e0+0.004):.0%}, "
                  f"of injected within 4 meV of it + hw: {weighted(E, inj, e0+HW-0.004, e0+HW+0.004):.0%}", flush=True)
