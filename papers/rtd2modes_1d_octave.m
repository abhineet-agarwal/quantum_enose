% rtd2modes_1d_octave.m
%
% Headless Octave port of papers/rtd2modes_1d.m for the SISPAD Tier-1 γ
% numerical reference. The physics is byte-identical to Patil's original
% (lines 22-26 + 44-48 of rtd2modes_1d.m). The differences from the
% original are exactly:
%
%   1. VV = linspace(0, 0.8, 51) instead of linspace(0, 0.3, 76).
%      Patil's published Fig. 2 plots out to V = 0.8 V, but the committed
%      .m file only sweeps to 0.3 V — the published figure must have used
%      a longer sweep. We use 51 points so the run is tractable in Octave
%      (~1-2 hours for the full SCBA loop).
%
%   2. The per-iV save() call inside the loop is removed (slow on Octave
%      and we only need the final arrays).
%
%   3. All figure() / plot() / saveas() calls are removed — this script
%      is run headless and the comparison plot is generated in Python by
%      run/plot_patil_overlay.py after loading the saved .mat.
%
%   4. The signal package is loaded explicitly with `pkg load signal` so
%      `butter`/`filtfilt` are available under Octave (under MATLAB they
%      live in the Signal Processing Toolbox).
%
%   5. The final `save` writes to a fixed path: tests/patil_octave_reference.mat
%      The variables saved are exactly those needed for the Python overlay:
%      VV, II, II3, II4, IIco, IInonco, G1, G2, g2, IETS, IETS2, VG, VETS,
%      plus the parameter block (Ef, kT, t0, Np, Vb, Nb, Dnu, hnu, dE, NV).
%
% Run with:
%
%     octave --no-gui --eval "run('papers/rtd2modes_1d_octave.m')"
%
% from the project root (so the relative save path resolves correctly).

clear all
close all
pkg load signal

%Constants (all MKS, except energy which is in eV)
hbar = 1.06e-34; q = 1.6e-19; m = 0.067 * 9.1e-31;
IE = (q * q) / (2 * pi * hbar);
Ef = 0.02; kT = 0.025;
t0 = 5.2;

a = sqrt((hbar^2) / (2 * m * (t0) * q));

%%%%%%%%%%%%%%%%%%%%%% one peak in TM vs E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NS = 1; NC = 40; ND = 1; Np = NS + NC + ND;
Nb = 15; Vb = 0.6;
UB = [zeros(NS, 1); Vb * ones(Nb, 1); zeros(NC - 2 * Nb, 1); ...
      Vb * ones(Nb, 1); zeros(ND, 1)];

T = (2 * t0 * diag(ones(1, Np))) ...
    - (t0 * diag(ones(1, Np - 1), 1)) ...
    - (t0 * diag(ones(1, Np - 1), -1));
T = T + diag(UB);

%Bias  -- PATCH 1: extend sweep to 0.8 V to match Patil Fig. 2
NV = 51;
VV = linspace(0, 0.8, NV);
dV = VV(2) - VV(1);
dE = 0.005;
D0 = 0;
Dnu = [0.1, 0.1];
Nph = size(Dnu, 2);
hnu = [18, 35];
Nhnu = 1 ./ ((exp(dE * hnu ./ kT)) - 1);

II = zeros(1, NV);
II3 = zeros(1, NV);
II4 = zeros(1, NV);
IIco = zeros(1, NV);
IInonco = zeros(1, NV);

% PATCH 6: write to an ABSOLUTE path so save() never depends on cwd, and
% checkpoint after every bias point so a transient FS hiccup at the end
% does not lose an hour of compute. The earlier 2026-04-08 run reached
% iV=51/51 then crashed at the final save.
out_path = '/Users/abhineet/Downloads/quantum_enose/tests/patil_octave_reference.mat';

t_start = tic;
for iV = 1:NV
    V = VV(iV); mu1 = Ef + (V / 2); mu2 = Ef - (V / 2);
    U1 = V * [.5 * ones(1, NS) linspace(0.5, -0.5, NC) -.5 * ones(1, ND)];
    U1 = U1';
    E = [-0.2:dE:0.8]; NE = size(E, 2);
    zplus = i * 1e-12;
    fprintf('iV=%d/%d  V=%.4f  elapsed=%.0fs\n', iV, NV, V, toc(t_start));
    fflush(stdout);

    f1 = 1 ./ (1 + exp((E - mu1) ./ kT));
    f2 = 1 ./ (1 + exp((E - mu2) ./ kT));

    sigin1 = zeros(Np, Np, NE); sigout1 = zeros(Np, Np, NE);
    sigin2 = zeros(Np, Np, NE); sigout2 = zeros(Np, Np, NE);
    sigin = 0 * ones(Np, Np, NE); sigout = 0 * ones(Np, Np, NE);
    n = zeros(Np, Np, NE);
    p = zeros(Np, Np, NE);
    gamp = zeros(Np, Np, NE);
    gam1 = zeros(Np, Np, NE);
    gam2 = zeros(Np, Np, NE);
    G = zeros(Np, Np, NE);

    change = 1;
    it = 1 / 2;
    iter_count = 0;
    while change > 1e-5 && iter_count < 200
        iter_count = iter_count + 1;
        for k = 1:NE
            sig1 = zeros(Np); sig2 = zeros(Np);
            ck = 1 - ((E(k) + zplus - U1(1) - UB(1)) / (2 * t0));
            ka = acos(ck);
            sig1(1, 1) = -t0 * exp(i * ka);
            gam1(:, :, k) = i * (sig1 - sig1');
            ck = 1 - ((E(k) + zplus - U1(Np) - UB(Np)) / (2 * t0));
            ka = acos(ck);
            sig2(Np, Np) = -t0 * exp(i * ka);
            gam2(:, :, k) = i * (sig2 - sig2');
            sigin1(:, :, k) = f1(k) * gam1(:, :, k);
            sigin2(:, :, k) = f2(k) * gam2(:, :, k);
            sigout1(:, :, k) = (1 - f1(k)) * gam1(:, :, k);
            sigout2(:, :, k) = (1 - f2(k)) * gam2(:, :, k);
            gamp(:, :, k) = sigin(:, :, k) + sigout(:, :, k);
            G(:, :, k) = inv(sparse(((E(k) + zplus) * eye(Np)) - T ...
                                    - diag(U1) - sig1 - sig2 ...
                                    + (i * 0.5 * gamp(:, :, k))));
            A = i * (G(:, :, k) - G(:, :, k)');
            n(:, :, k) = real(G(:, :, k) * ((f1(k) * gam1(:, :, k)) ...
                              + (f2(k) * gam2(:, :, k)) ...
                              + sigin(:, :, k)) * G(:, :, k)');
            p(:, :, k) = A - n(:, :, k);
        end
        siginnew = D0 * n; sigoutnew = D0 * p;
        for iph = 1:Nph
            ne = zeros(Np, Np, NE);
            na = zeros(Np, Np, NE);
            pe = zeros(Np, Np, NE);
            pa = zeros(Np, Np, NE);
            inu = hnu(iph);
            if inu < NE
                for nn = 1:Np
                    ne(nn, nn, [1:NE - inu]) = n(nn, nn, [inu + 1:NE]);
                    na(nn, nn, [inu + 1:NE]) = n(nn, nn, [1:NE - inu]);
                    pe(nn, nn, [1:NE - inu]) = p(nn, nn, [inu + 1:NE]);
                    pa(nn, nn, [inu + 1:NE]) = p(nn, nn, [1:NE - inu]);
                end
                siginnew = siginnew + ((Nhnu(iph) + 1) * Dnu(iph) * ne) ...
                                    + (Nhnu(iph) * Dnu(iph) * na);
                sigoutnew = sigoutnew + (Nhnu(iph) * Dnu(iph) * pe) ...
                                      + ((Nhnu(iph) + 1) * Dnu(iph) * pa);
            end
        end
        change = sum(sum(sum(abs(siginnew - sigin)))) ...
               + sum(sum(sum(abs(sigoutnew - sigout))));
        sigin = ((1 - it) * sigin) + (it * siginnew);
        sigout = ((1 - it) * sigout) + (it * sigoutnew);
    end

    I1 = 0; I3 = 0; Ico = 0; Inco = 0; I2 = 0;
    for k = 1:NE
        I1 = I1 + real(trace((sigout2(:, :, k) * n(:, :, k)) ...
                             - (sigin2(:, :, k) * p(:, :, k))));
        I3 = I3 + real(trace((sigout(:, :, k) * n(:, :, k)) ...
                             - (sigin(:, :, k) * p(:, :, k))));
        I2 = I2 + real(trace((sigout1(:, :, k) * n(:, :, k)) ...
                             - (sigin1(:, :, k) * p(:, :, k))));
        Ico = Ico + real(trace((sigin2(:, :, k) * G(:, :, k) * gam1(:, :, k) * G(:, :, k)') ...
                               - (gam2(:, :, k) * G(:, :, k) * sigin1(:, :, k) * G(:, :, k)')));
        Inco = Inco + real(trace((sigin2(:, :, k) * G(:, :, k) * gamp(:, :, k) * G(:, :, k)') ...
                                 - (gam2(:, :, k) * G(:, :, k) * sigin(:, :, k) * G(:, :, k)')));
    end
    II(iV) = sum(I1) * dE * IE;
    II3(iV) = sum(I3) * dE * IE;
    II4(iV) = sum(I2) * dE * IE;
    IIco(iV) = sum(Ico) * dE * IE;
    IInonco(iV) = sum(Inco) * dE * IE;

    % PATCH 6 cont'd: incremental checkpoint. Saves the current arrays
    % every bias point so a crash leaves a usable partial reference.
    try
        save('-v7', out_path, 'VV', 'II', 'II3', 'II4', 'IIco', 'IInonco', ...
             'iV', 'Ef', 'kT', 't0', 'Np', 'Vb', 'Nb', 'Dnu', 'hnu', 'dE', 'NV');
    catch err
        fprintf('checkpoint save failed at iV=%d: %s\n', iV, err.message);
    end
end

G1 = diff(II) ./ dV; VG = VV([2:NV]);
IETS = diff(G1) ./ dV; VETS = VV([3:NV]);

[bf, af] = butter(2, 0.1);
II2 = filtfilt(bf, af, II);
G2 = diff(II2) ./ dV;
g2 = filtfilt(bf, af, G2);
IETS2 = diff(g2) ./ dV;

% PATCH 5: final save with derived quantities (G1, G2, g2, IETS, IETS2)
% added on top of the per-iV checkpoints. Uses the same absolute path.
save('-v7', out_path, ...
     'VV', 'II', 'II3', 'II4', 'IIco', 'IInonco', ...
     'G1', 'G2', 'g2', 'IETS', 'IETS2', 'VG', 'VETS', ...
     'Ef', 'kT', 't0', 'Np', 'Vb', 'Nb', 'Dnu', 'hnu', 'dE', 'NV');
fprintf('saved %s\n', out_path);
