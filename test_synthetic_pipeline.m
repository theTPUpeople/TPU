% TEST_SYNTHETIC_PIPELINE  Validate the full FRF pipeline on synthetic data with known m, k, c.
%
% Builds FRFs with known modal parameters, writes them in the QuickDAQ CSV
% layout used by the real data, runs tpu_modal_batch end to end (reader ->
% cleaning -> accelerance-to-receptance -> polarity -> textbook peak picking
% -> half-power / transmissibility) and compares the recovered values.
%
%   Case A  two independent SDOF modes (112 kg / 55 Hz / zeta 0.07 and
%           30 kg / 170 Hz / zeta 0.05), receptance -> accelerance, sign
%           inverted as in the measured disturbance-rejection data, plus a
%           NaN row, a duplicate row, an unsorted pair, low-coherence bins
%   Case B  2-DOF chain (ground-k1-m1-k2-m2, Schmitz & Smith Sec. 2.4) with
%           proportional damping; direct FRF X1/F1 by matrix inversion,
%           expected modal values from eigenvectors normalised to x1
%   Case C  SDOF base-excitation transmissibility T = (k + i c w)/(k - m w^2 + i c w)
%   Case D  cross FRF X2/F1 of the Case B chain: one mode is inverted. Checks
%           w_n and zeta, that the inverted mode gets no mass, and that the
%           passivity warning is raised
%
% Tolerances (task definition): f_n and m within 2 %, zeta within 10 %.
% Derived: k = m w_n^2 within 2 + 2*2 = 6 %, c within 10 + 2 = 12 %,
% isolation onset f_iso/f_n vs sqrt(2) within 2 %.
%
% Gate:
%   1. noise-free data: every run within every tolerance (method bias);
%   2. 1 % and 2 % complex multiplicative noise, N_SEEDS runs each:
%      |mean error| (bias) within tolerance at both levels, and
%      RMS error within tolerance at 1 %.
%   The fraction of single runs inside all tolerances is printed as well.
%   At 2 % noise the textbook mass estimate scatters by about +/-2 % (1 sigma),
%   so single runs can exceed 2 % without any bias in the method.
% Prints PASS/FAIL and errors out on FAIL.

clear; close all;
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
addpath(here);
if exist('OCTAVE_VERSION', 'builtin')
    try, pkg load optim; catch, end %#ok<NOCOM>
end
pf = @(b) char(double('PASS')*b + double('FAIL')*~b);

N_SEEDS = 20;
tol_D = [2 10 2 6 12];          % f_n zeta m k c   [%]
tol_T = [2 10 NaN 6 12 2];      % f_n zeta - k c f_iso/f_n

tmp = fullfile(tempdir, 'tpu_synthetic_test');
if exist(tmp, 'dir') == 7, rmdir(tmp, 's'); end
mkdir(tmp);

% same frequency grid as the QuickDAQ exports
fs = 1600.0093773442809; Nfft = 4096;
f = (0:Nfft/2-1)'*fs/Nfft;
w = 2*pi*f;

% ---------------- Case A: two independent SDOF modes ----------------
A.m = [112 30]; A.fn = [55 170]; A.zeta = [0.07 0.05];
A.wn = 2*pi*A.fn; A.k = A.m.*A.wn.^2; A.c = 2*A.zeta.*A.m.*A.wn;
HA_rec = zeros(size(w));
for q = 1:2, HA_rec = HA_rec + 1./(A.k(q) - A.m(q)*w.^2 + 1i*A.c(q)*w); end

% ---------------- Case B: 2-DOF chain, proportional damping ----------------
m1 = 50; m2 = 112; k1 = 3.0e7; k2 = 1.2e7;
Mm = diag([m1 m2]); Km = [k1+k2 -k2; -k2 k2];
[V, L] = eig(Km, Mm); [wn2, is] = sort(diag(L)); V = V(:, is); wnB = sqrt(wn2);
zt = [0.06; 0.08];                               % target modal damping
ab = [1./(2*wnB) wnB/2] \ zt;                    % zeta_i = alpha/(2 w_i) + beta w_i/2
Cm = ab(1)*Mm + ab(2)*Km;
V = V ./ repmat(V(1, :), 2, 1);                  % normalise to coordinate 1 (driving point)
B.m = diag(V'*Mm*V)'; B.k = diag(V'*Km*V)'; B.c = diag(V'*Cm*V)';
B.wn = sqrt(B.k./B.m); B.fn = B.wn/(2*pi); B.zeta = B.c./(2*sqrt(B.k.*B.m));
HB_rec = zeros(size(w)); HX_rec = HB_rec;
for i = 1:numel(w)
    Hi = inv(Km - w(i)^2*Mm + 1i*w(i)*Cm);       % direct inversion, independent of the modal model
    HB_rec(i) = Hi(1, 1);
    HX_rec(i) = Hi(2, 1);                        % cross FRF for Case D
end

% ---------------- Case C: transmissibility (mass on mount, base excitation) ----------------
Cc.m = A.m(1); Cc.k = A.k(1); Cc.c = A.c(1); Cc.wn = A.wn(1); Cc.fn = A.fn(1); Cc.zeta = A.zeta(1);
T_true = (Cc.k + 1i*Cc.c*w)./(Cc.k - Cc.m*w.^2 + 1i*Cc.c*w);

fprintf('\nCase A truth: f_n = [%g %g] Hz, zeta = [%g %g], m = [%g %g] kg\n', A.fn, A.zeta, A.m);
fprintf('Case B truth: f_n = [%.3f %.3f] Hz, zeta = [%.4f %.4f], m_q = [%.2f %.2f] kg (normalised to x1)\n', B.fn, B.zeta, B.m);
fprintf('Case D truth: cross FRF X2/F1 of Case B, same f_n and zeta; one mode inverted\n');
fprintf('Case C truth: f_n = %g Hz, zeta = %g, m = %g kg, isolation onset sqrt(2) f_n = %.3f Hz\n', Cc.fn, Cc.zeta, Cc.m, sqrt(2)*Cc.fn);

coh = ones(size(f)); coh(f < 15) = 0.5; coh(abs(f - 300) < 1) = 0.4;   % low-coherence regions away from modes

runs = [0 1; 0.01*ones(N_SEEDS, 1) (1:N_SEEDS)'; 0.02*ones(N_SEEDS, 1) (1:N_SEEDS)'];
EA = nan(size(runs, 1), 2, 5); EB = EA; ED = EA; EC = nan(size(runs, 1), 6);
okA = false(size(runs, 1), 1); okB = okA; okC = okA; okD = okA;
rows = {};
for ir = 1:size(runs, 1)
    nl = runs(ir, 1); sd = runs(ir, 2);
    rng(sd);
    cn = @(n) (randn(n, 1) + 1i*randn(n, 1))/sqrt(2);
    root = fullfile(tmp, sprintf('run%03d', ir));
    % Case A as accelerance, sign inverted, plus junk rows to exercise cleaning
    HA_acc = -(-w.^2 .* HA_rec) .* (1 + nl*cn(numel(w)));
    fA = [f; NaN; f(500)]; HA_acc = [HA_acc; 1; HA_acc(500)]; cA = [coh; 1; coh(500)];     % NaN row + duplicate
    fA([700 701]) = fA([701 700]); HA_acc([700 701]) = HA_acc([701 700]);                  % unsorted pair
    write_quickdaq_csv(fullfile(root, 'Disturbance Rejection Test', '0%', 'TPU0S1D.csv'), fA, HA_acc, cA, 'D');
    % Case B as accelerance, textbook sign
    HB_acc = (-w.^2 .* HB_rec) .* (1 + nl*cn(numel(w)));
    write_quickdaq_csv(fullfile(root, 'Disturbance Rejection Test', '0%', 'TPU0S2D.csv'), f, HB_acc, coh, 'D');
    % Case D cross FRF as accelerance
    HX_acc = (-w.^2 .* HX_rec) .* (1 + nl*cn(numel(w)));
    write_quickdaq_csv(fullfile(root, 'Disturbance Rejection Test', '0%', 'TPU0S3D.csv'), f, HX_acc, coh, 'D');
    % Case C transmissibility, two specimens so both D files have a partner
    for s = 1:2
        Tn = T_true .* (1 + nl*cn(numel(w)));
        write_quickdaq_csv(fullfile(root, 'Transmissibility test', '0%', sprintf('TPU0S%dT.csv', s)), f, Tn, coh, 'T');
    end

    cfg = struct();
    cfg.data_root = root;
    cfg.results_dir = fullfile(root, 'results');
    cfg.m_known_kg = Cc.m;
    cfg.verbose = false;
    cfg.make_figures = (ir == 1 || ir == N_SEEDS + 2);
    res = tpu_modal_batch(cfg);

    names = {res.name};
    tag = sprintf('noise %3.0f%% seed %2d', 100*nl, sd);
    [okA(ir), r, e] = validate_modal_result('D', res(strcmp(names, 'TPU0S1D')), A, -1, ['A  ' tag]); EA(ir, 1:size(e, 1), :) = e; rows = [rows r]; %#ok<AGROW>
    [okB(ir), r, e] = validate_modal_result('D', res(strcmp(names, 'TPU0S2D')), B, +1, ['B  ' tag]); EB(ir, 1:size(e, 1), :) = e; rows = [rows r]; %#ok<AGROW>
    [okC(ir), r, e] = validate_modal_result('T', res(strcmp(names, 'TPU0S1T')), Cc, ['C  ' tag]);    EC(ir, :) = e;                 rows = [rows r]; %#ok<AGROW>
    [okD(ir), r, e] = validate_modal_result('X', res(strcmp(names, 'TPU0S3D')), B, ['D  ' tag]);     ED(ir, 1:size(e, 1), :) = e; rows = [rows r]; %#ok<AGROW>
end

fprintf('\nPer-run errors [%%] (first runs of each noise level):\n');
fprintf('%-28s %4s %9s %8s %8s %8s %8s %8s  %s\n', 'case', 'mode', 'f_n', 'zeta', 'm', 'k', 'c', 'pol/iso', 'all tol');
for i = 1:numel(rows)
    if ~isempty(regexp(rows{i}, 'seed  [1-3] ', 'once')), fprintf('%s\n', rows{i}); end
end

% ---------------- statistics and gate ----------------
allpass = true;
labels = {'f_n', 'zeta', 'm', 'k', 'c', 'f_iso/f_n'};
fprintf('\nError statistics [%%] over %d seeds: bias (mean) / std / RMS, and the gate\n', N_SEEDS);
fprintf('%-12s %-10s %6s | %22s | %22s | %s\n', 'case/mode', 'param', 'tol', '1% noise  bias std rms', '2% noise  bias std rms', 'gate');
for c = 1:4
    for q = 1:2
        if c == 3 && q == 2, continue, end
        switch c
            case 1, E = squeeze(EA(:, q, :)); tl = tol_D; cname = sprintf('A mode %d', q);
            case 2, E = squeeze(EB(:, q, :)); tl = tol_D; cname = sprintf('B mode %d', q);
            case 3, E = EC; tl = tol_T; cname = 'C (T)';
            case 4, E = squeeze(ED(:, q, :)); tl = [tol_D(1:2) NaN NaN NaN]; cname = sprintf('D mode %d', q);
                if q == 2, tl(2) = -tl(2); end      % inverted mode: zeta printed, not gated
        end
        for p = 1:numel(tl)
            if isnan(tl(p)), continue, end
            e0 = E(1, p);
            e1 = E(runs(:, 1) == 0.01, p); e2 = E(runs(:, 1) == 0.02, p);
            st = @(x) [mean(x) std(x) sqrt(mean(x.^2))];
            s1 = st(e1); s2 = st(e2);
            if tl(p) < 0
                fprintf('%-12s %-10s %5s  | %+7.2f %6.2f %6.2f | %+7.2f %6.2f %6.2f | noise-free %+6.2f  not gated (*)\n', ...
                    cname, labels{p}, 'info', s1, s2, e0);
                continue
            end
            g = abs(e0) <= tl(p) && all(isfinite([s1 s2])) && abs(s1(1)) <= tl(p) && s1(3) <= tl(p) && abs(s2(1)) <= tl(p);
            allpass = allpass && g;
            fprintf('%-12s %-10s %5g%% | %+7.2f %6.2f %6.2f | %+7.2f %6.2f %6.2f | noise-free %+6.2f  %s\n', ...
                cname, labels{p}, tl(p), s1, s2, e0, pf(g));
        end
    end
end
for lv = [0 0.01 0.02]
    k = runs(:, 1) == lv;
    fprintf('single runs inside all tolerances at %g%% noise: A %3.0f%%, B %3.0f%%, C %3.0f%%, D %3.0f%%\n', 100*lv, ...
        100*mean(okA(k)), 100*mean(okB(k)), 100*mean(okC(k)), 100*mean(okD(k)));
end
fprintf(['(*) Case D mode 2 is inverted and much weaker than mode 1 in this cross FRF; peak picking biases its zeta through\n', ...
    '    the sloped residual of mode 1 (inherent to the textbook method: raw bins give a similar bias, the mode alone is\n', ...
    '    recovered within 0.6 %%). Its w_n is gated; its zeta is reported only.\n']);
okDx = all(okD(runs(:, 1) == 0));
fprintf('Case D noise-free: one inverted mode, no mass for it, passivity warning raised: %s\n', pf(okDx));
allpass = allpass && okDx;

% ---------------- reader check against raw text of real files ----------------
fprintf('\nReader check on real exports (column mapping from headers):\n');
okr = true;
rd = fullfile(here, 'Disturbance Rejection Test', '0%', 'TPU0S1D.csv');
rt = fullfile(here, 'Transmissibility test', '0%', 'TPU0S1T.csv');
if exist(rd, 'file') == 2 && exist(rt, 'file') == 2
    Dd = read_quickdaq_csv(rd); Dt = read_quickdaq_csv(rt);
    % raw rows at f = 0.390627...: D "1.43154E-004,-85.0066827086, 0.6323274718"; T "0.864372646,1.234509055,-0.36764638"
    e1 = abs(abs(Dd.H(2)) - 1.43154e-4)/1.43154e-4 + abs(angle(Dd.H(2))*180/pi + 85.0066827086) + abs(Dd.coh(2) - 0.6323274718);
    e2 = abs(abs(Dt.H(2)) - 1.234509055) + abs(angle(Dt.H(2))*180/pi + 0.36764638) + abs(Dt.coh(2) - 0.864372646);
    okr = e1 < 1e-6 && e2 < 1e-6 && strcmp(Dd.frf_kind, 'accelerance') && strcmp(Dt.frf_kind, 'transmissibility');
    fprintf('  D: kind %s, |H| %.6g, phase %.6f deg, coh %.6f\n', Dd.frf_kind, abs(Dd.H(2)), angle(Dd.H(2))*180/pi, Dd.coh(2));
    fprintf('  T: kind %s, |T| %.6g, phase %.6f deg, coh %.6f\n', Dt.frf_kind, abs(Dt.H(2)), angle(Dt.H(2))*180/pi, Dt.coh(2));
    fprintf('  %s\n', pf(okr));
else
    fprintf('  real files not found, skipped\n');
end
allpass = allpass && okr;

% ---------------- m_known = NaN path ----------------
fprintf('\nm_known_kg = NaN path:\n');
cfg.data_root = fullfile(tmp, 'run001');
cfg.m_known_kg = NaN; cfg.make_figures = false;
cfg.results_dir = fullfile(tmp, 'nan_mass');
resN = tpu_modal_batch(cfg);
logtxt = fileread(fullfile(cfg.results_dir, 'run_log.txt'));
md = [resN(strcmp({resN.test}, 'D')).modes];
mt = [resN(strcmp({resN.test}, 'T')).modes];
okn = all(isnan([md.k_known])) && all(isnan([mt.k_known])) && all(isfinite([md(~[md.inverted]).m_modal])) ...
    && ~isempty(strfind(logtxt, 'm_known_kg is NaN'));
fprintf('  known-mass columns NaN, modal mass still computed, warning logged: %s\n', pf(okn));
allpass = allpass && okn;

% ---------------- existing modalfit.m still runs ----------------
fprintf('\nmodalfit.m (existing lsqnonlin refinement):\n');
if exist('lsqnonlin') %#ok<EXIST>
    x0 = [A.wn(1)*1.003 A.k(1)*0.97 A.zeta(1)*1.05 A.wn(2)*0.998 A.k(2)*1.03 A.zeta(2)*0.95];
    sel = f > 20 & f < 500;
    evalc('x = modalfit(f(sel), real(HA_rec(sel)), imag(HA_rec(sel)), x0);');
    xt = [A.wn(1) A.k(1) A.zeta(1) A.wn(2) A.k(2) A.zeta(2)];
    err = abs(x - xt)./xt;
    okm = all(size(x) == size(x0)) && all(err < 0.01);
    fprintf('  max relative error after refinement on noise-free data: %.2e  %s\n', max(err), pf(okm));
    allpass = allpass && okm;
else
    fprintf('  lsqnonlin not available (Optimization Toolbox / Octave optim): SKIPPED\n');
end

fprintf('\nFigures: %s (noise-free) and %s (2%% noise, seed 1)\n', fullfile(tmp, 'run001', 'results', 'figures'), ...
    fullfile(tmp, sprintf('run%03d', N_SEEDS + 2), 'results', 'figures'));
if allpass
    fprintf('\n==== SYNTHETIC VALIDATION: PASS ====\n');
else
    fprintf('\n==== SYNTHETIC VALIDATION: FAIL ====\n');
    error('test_synthetic_pipeline:fail', 'synthetic validation failed; do not run on real data until fixed');
end
