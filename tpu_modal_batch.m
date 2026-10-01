function [results, xcheck] = tpu_modal_batch(cfg)
%TPU_MODAL_BATCH Textbook peak-picking modal identification for all TPU specimens.
%   [RESULTS, XCHECK] = TPU_MODAL_BATCH(CFG) processes every QuickDAQ CSV
%   export under CFG.data_root/<test folder>/<group>/ for the two tests
%
%     D  Disturbance rejection: impact hammer on the mass, force input,
%        acceleration output (accelerance). Converted to receptance
%        H_rec = H_acc/(-w^2), then peak picking per Schmitz & Smith,
%        Machining Dynamics 2nd ed., Sec. 2.5.1, Eqs. 2.58-2.63 -> w_n, zeta,
%        A/B, k_q, m_q (modal mass), c_q.
%     T  Transmissibility: base acceleration input, mass acceleration output
%        (dimensionless). w_n from peak |T|, zeta from the half-power
%        bandwidth of |T|, isolation onset where |T| < 1. No force -> no mass;
%        k and c use CFG.m_known_kg only ("mass-assumed").
%
%   Writes to CFG.results_dir:
%     modal_results.csv      one row per specimen x test x mode
%     crosscheck_D_vs_T.csv  w_n and zeta compared between the two tests
%     modal_results.mat      RESULTS, XCHECK, CFG (incl. [K_q], [M_q], [C_q])
%     run_log.txt            parameters, conversions, assumptions, warnings
%     figures/<file>.png     Re/Im stacked plots after textbook Fig. 2.21
%
%   See FRF_Modal_Parameter_Fitting_Optimization.m for the CFG fields.

cfg = apply_defaults(cfg);
if exist(cfg.results_dir, 'dir') ~= 7, mkdir(cfg.results_dir); end
figdir = fullfile(cfg.results_dir, 'figures');
if cfg.make_figures && exist(figdir, 'dir') ~= 7, mkdir(figdir); end
fid = fopen(fullfile(cfg.results_dir, 'run_log.txt'), 'w');
if fid < 0, error('tpu_modal_batch:log', 'cannot open run_log.txt in %s', cfg.results_dir); end
lastwarn('');

try
    [results, xcheck] = run_all(cfg, fid, figdir);
catch err
    lg(fid, cfg, 'ERROR', '%s', err.message);
    fclose(fid);
    rethrow(err);
end
fclose(fid);
end

% =========================================================================
function [results, xcheck] = run_all(cfg, fid, figdir)
if exist('OCTAVE_VERSION', 'builtin'), env = ['GNU Octave ' OCTAVE_VERSION]; else, env = ['MATLAB ' version]; end
lg(fid, cfg, '', '=== TPU modal identification run: %s (%s) ===', datestr(now, 'yyyy-mm-dd HH:MM:SS'), env);

% ---------------- parameters ----------------
lg(fid, cfg, '', '');
lg(fid, cfg, '', '--- PARAMETERS ---');
lg(fid, cfg, '', 'data_root          = %s', cfg.data_root);
lg(fid, cfg, '', 'm_known_kg         = %g', cfg.m_known_kg);
lg(fid, cfg, '', 'frf_type           = %s   (force-input dataset)', cfg.frf_type);
lg(fid, cfg, '', 'n_modes            = %d', cfg.n_modes);
if isempty(cfg.min_prominence)
    lg(fid, cfg, '', 'min_prominence     = []  -> %.3g x max|Im| per FRF (min_prominence_rel)', cfg.min_prominence_rel);
else
    lg(fid, cfg, '', 'min_prominence     = %g (absolute, FRF units)', cfg.min_prominence);
end
lg(fid, cfg, '', 'polarity           = %s', num2str(cfg.polarity));
lg(fid, cfg, '', 'f_band_Hz          = [%g %g]', cfg.f_band_Hz);
lg(fid, cfg, '', 'min_coherence      = %g', cfg.min_coherence);
lg(fid, cfg, '', 'smooth_points      = %d (D, box), smooth_points_T = %d (T, %s)', cfg.smooth_points, cfg.smooth_points_T, cfg.smooth_method_T);
lg(fid, cfg, '', 'peak_refine        = %s', num2str(cfg.peak_refine));
lg(fid, cfg, '', 'real_search_frac   = %g', cfg.real_search_frac);
lg(fid, cfg, '', 'disagree_pct       = %g %%', cfg.disagree_pct);
lg(fid, cfg, '', 'phase_units        = %s', cfg.phase_units);

% ---------------- method ----------------
lg(fid, cfg, '', '');
lg(fid, cfg, '', '--- METHOD (Schmitz & Smith, Machining Dynamics 2nd ed., Sec. 2.5.1, pp. 42-43) ---');
lg(fid, cfg, '', 'Fig. 2.21 points: 1,2 = minima of Im(H_rec) -> w_n1, w_n2; 3/4 = Re max/min around mode 1 (w3, w4);');
lg(fid, cfg, '', '  5/6 = Re max/min around mode 2 (w5, w6); A, B = Im(H_rec) at points 1, 2. Columns f_a/f_b = w3,w5 / w4,w6.');
lg(fid, cfg, '', 'Eq. 2.58/2.59: zeta_q = (w_b - w_a)/(2 w_n);  Eq. 2.60/2.61: k_q = -1/(2 zeta_q A);');
lg(fid, cfg, '', 'Eq. 2.62: m_q = k_q/w_n^2  (= 1/(2 zeta w_n^2 |A|) = m_modal);  Eq. 2.63: c_q = 2 zeta_q sqrt(k_q m_q).');
lg(fid, cfg, '', 'Known-mass values: k = m_known w_n^2, c = 2 zeta m_known w_n;  %%diff = 100 (m_modal - m_known)/m_known.');
lg(fid, cfg, '', 'Cross-check (D): half-power bandwidth of |H_rec|.  T: real/imag peak picking on G = (T-1)/w^2,');
lg(fid, cfg, '', '  which for a base-excited viscous SDOF equals m*H_rec (T = 1 + m w^2 H_rec), so w_n and zeta follow');
lg(fid, cfg, '', '  from Eq. 2.58 without a mass. Disagreements > disagree_pct are flagged.');
lg(fid, cfg, '', 'Isolation onset: first frequency ABOVE the highest |T| peak where |T| < 1 (below resonance |T| < 1 also occurs in');
lg(fid, cfg, '', '  these data, which a base-excited SDOF cannot produce; see the per-file model check).');
lg(fid, cfg, '', 'T smoothing: %d-point %s filter on complex T before peak / half-power picking (T only; D unsmoothed).', cfg.smooth_points_T, cfg.smooth_method_T);
lg(fid, cfg, '', 'Sub-bin refinement (peak_refine = ''sdof''): the picked bins are refined by a local least-squares parabola');
lg(fid, cfg, '', '  in the coordinates where the SDOF curves are symmetric: Im vs x = 1-(w/w_n)^2 (|x| <= 0.35*2*zeta),');
lg(fid, cfg, '', '  Re vs ln|1-(w/w_n)^2| (+/- ln 2). Raw-bin zeta and m_q are also reported (zeta_rawbins, m_modal_rawbins).');

% ---------------- known mass ----------------
lg(fid, cfg, '', '');
if ~isfinite(cfg.m_known_kg)
    lg(fid, cfg, 'WARN', 'm_known_kg is NaN: known-mass k, c, %%diff and all transmissibility k, c are skipped (NaN).');
else
    lg(fid, cfg, '', 'm_known_kg = %g kg is a user-supplied constant (not read from the data or config files).', cfg.m_known_kg);
end

xml_cfg = struct();
results = new_result();
results = results([]);

for t = 1:numel(cfg.tests)
    spec = cfg.tests(t);
    folder = fullfile(cfg.data_root, spec.folder);
    lg(fid, cfg, '', '');
    lg(fid, cfg, '', '--- TEST %s: %s ---', spec.id, folder);
    if exist(folder, 'dir') ~= 7
        lg(fid, cfg, 'WARN', 'folder not found, test %s skipped', spec.id);
        continue
    end

    % acquisition settings (documentation only)
    xf = dir(fullfile(folder, '*config.xml'));
    C = [];
    if ~isempty(xf)
        C = read_quickdaq_config(fullfile(folder, xf(1).name));
        log_acquisition(fid, cfg, C, spec);
    end
    xml_cfg.(spec.id) = C;

    files = list_csv(folder);
    lg(fid, cfg, '', '%d CSV files found', numel(files));
    for i = 1:numel(files)
        R = process_file(files{i}, spec, C, cfg, fid, figdir);
        results(end+1) = R; %#ok<AGROW>
    end
end

% ---------------- cross-check D vs T ----------------
xcheck = cross_check(results, cfg, fid);

% ---------------- outputs ----------------
write_results_csv(fullfile(cfg.results_dir, 'modal_results.csv'), results);
write_xcheck_csv(fullfile(cfg.results_dir, 'crosscheck_D_vs_T.csv'), xcheck);
save(fullfile(cfg.results_dir, 'modal_results.mat'), 'results', 'xcheck', 'cfg', 'xml_cfg', '-v7');
print_table(fid, cfg, results);

% ---------------- warnings summary ----------------
lg(fid, cfg, '', '');
lg(fid, cfg, '', '--- WARNINGS SUMMARY ---');
nw = 0;
for i = 1:numel(results)
    for j = 1:numel(results(i).warnings)
        lg(fid, cfg, '', '%-10s %s', results(i).name, results(i).warnings{j});
        nw = nw + 1;
    end
end
for j = 1:numel(xcheck)
    if ~isempty(xcheck(j).flag)
        lg(fid, cfg, '', '%-10s D vs T: %s', sprintf('%sS%d', xcheck(j).group, xcheck(j).specimen), xcheck(j).flag);
        nw = nw + 1;
    end
end
lg(fid, cfg, '', '%d data/model warnings in total (all handled: values flagged, not dropped).', nw);
[msg, id] = lastwarn;
if ~isempty(msg)
    lg(fid, cfg, 'WARN', 'runtime warning raised during the run: %s (%s)', msg, id);
else
    lg(fid, cfg, '', 'No MATLAB/Octave runtime warnings were raised.');
end

log_limitations(fid, cfg);
lg(fid, cfg, '', '');
lg(fid, cfg, '', 'Outputs: %s', cfg.results_dir);
end

% =========================================================================
function R = process_file(file, spec, C, cfg, fid, figdir)
R = new_result();
[~, name] = fileparts(file);
[parent] = fileparts(file);
[~, group] = fileparts(parent);
R.name = name;
R.file = relpath(file, cfg.data_root);
R.group = group;
R.test = spec.id;
tok = regexp(name, '^TPU(\d+)S(\d+)([DT])$', 'tokens', 'once');
if ~isempty(tok)
    R.specimen = str2double(tok{2});
    if ~strcmpi(tok{3}, spec.id)
        error('tpu_modal_batch:ambiguous', ...
            'AMBIGUOUS: %s is in the "%s" folder (test %s) but its name ends in %s. Fix the file or folder before running.', ...
            R.file, spec.folder, spec.id, tok{3});
    end
else
    R.specimen = NaN;
end

lg(fid, cfg, '', '');
lg(fid, cfg, '', '[%s] %s', name, R.file);
raw = read_quickdaq_csv(file, cfg.phase_units);
R.frf_kind_header = raw.frf_kind;
R.y_units = raw.y_units;
for j = 1:numel(raw.notes), lg(fid, cfg, '', '  read: %s', raw.notes{j}); end
if ~isempty(raw.ignored_cols)
    lg(fid, cfg, '', '  read: ignored non-FRF columns: %s', strjoin(raw.ignored_cols, '; '));
end

% ---- test identity from headers, checked against folder / file name ----
switch spec.input
    case 'force'
        if ~any(strcmp(raw.frf_kind, {'accelerance', 'receptance'}))
            error('tpu_modal_batch:ambiguous', ...
                'AMBIGUOUS: %s is in the force-input folder but its Y units "%s" read as %s.', R.file, raw.y_units, raw.frf_kind);
        end
    case 'base'
        if ~strcmp(raw.frf_kind, 'transmissibility')
            error('tpu_modal_batch:ambiguous', ...
                'AMBIGUOUS: %s is in the transmissibility folder but its Y units "%s" read as %s.', R.file, raw.y_units, raw.frf_kind);
        end
end

if strcmp(spec.input, 'base')
    ns = cfg.smooth_points_T; sm = cfg.smooth_method_T;
else
    ns = cfg.smooth_points; sm = 'box';
end
opts = struct('f_band_Hz', cfg.f_band_Hz, 'min_coherence', cfg.min_coherence, 'smooth_points', ns, 'smooth_method', sm);
P = prepare_frf(raw.f_Hz, raw.H, raw.coh, opts);
for j = 1:numel(P.notes), lg(fid, cfg, '', '  prep: %s', P.notes{j}); end
R.n_band = numel(P.band_f_Hz);
R.n_used = numel(P.f_Hz);
R.n_coh_masked = sum(~P.band_used);

popts = struct('n_modes', cfg.n_modes, 'min_prominence', cfg.min_prominence, ...
    'min_prominence_rel', cfg.min_prominence_rel, 'real_search_frac', cfg.real_search_frac, ...
    'peak_refine', cfg.peak_refine);

if strcmp(spec.input, 'force')
    R = process_force(R, raw, P, popts, C, cfg, fid, figdir);
else
    R = process_base(R, P, popts, cfg, fid, figdir);
end
end

% =========================================================================
function R = process_force(R, raw, P, popts, C, cfg, fid, figdir)
% ---- accelerance -> receptance (explicit, switchable) ----
ftype = cfg.frf_type;
if strcmpi(ftype, 'auto')
    ftype = raw.frf_kind;
    src = sprintf('auto, from Y units "%s"', raw.y_units);
else
    src = 'set by cfg.frf_type';
    if ~strcmpi(ftype, raw.frf_kind)
        R = warn(R, fid, cfg, sprintf('cfg.frf_type = %s overrides header Y units "%s" (%s)', ftype, raw.y_units, raw.frf_kind));
    end
end
wb = 2*pi*P.band_f_Hz;
switch lower(ftype)
    case 'accelerance'
        G = accel_to_receptance(P.H, P.w);
        Gb = accel_to_receptance(P.band_H, wb);
        R.conversion = 'H_rec = H_acc/(-w^2)';
    case 'receptance'
        G = P.H; Gb = P.band_H;
        R.conversion = 'none (input already receptance)';
    otherwise
        error('tpu_modal_batch:frftype', 'frf_type "%s" not supported for the force-input test (use accelerance or receptance).', ftype);
end
R.frf_type_used = lower(ftype);
lg(fid, cfg, '', '  FRF type: %s (%s); conversion applied: %s', R.frf_type_used, src, R.conversion);

% ---- sign convention ----
[s, pinfo] = choose_polarity(P.w, G, cfg);
R.polarity = s;
G = s*G; Gb = s*Gb;
lg(fid, cfg, '', '  polarity: %+d (%s); Im<0 in %.0f%% of points with |H|>0.5 max; Re sign pattern consistent: %d', ...
    s, polarity_src(cfg), 100*pinfo.im_frac_negative, pinfo.re_consistent);
if s < 0
    R = warn(R, fid, cfg, 'sign inverted (x -1) to textbook convention: raw H_acc/(-w^2) had Im > 0 at resonance (hammer/accelerometer polarity opposed?)');
end
if ~pinfo.re_consistent
    R = warn(R, fid, cfg, 'Re(H_rec) sign pattern around the dominant peak does not match the textbook form (conjugate convention or non-driving-point FRF?)');
end

% ---- peak picking (Eqs. 2.58-2.63) ----
[modes, cand] = peak_pick_textbook(P.w, G, popts);
R.cand = cand;
log_candidates(fid, cfg, cand, 'Im(H_rec) minima', 'm/N');
if numel(modes) < cfg.n_modes
    R = warn(R, fid, cfg, sprintf('only %d of n_modes = %d modes passed the prominence threshold', numel(modes), cfg.n_modes));
end

% exponential response window -> added apparent damping (estimate only)
tau = NaN;
if ~isempty(C) && strcmpi(C.WindowTypeResp, 'Exponential') && C.ExpWindowFinalValue > 0 && C.ExpWindowFinalValue < 1
    tau = (C.FFTSize/C.SampleRate)/log(1/C.ExpWindowFinalValue);
end

recs = new_mode();
recs = recs([]);
for q = 1:numel(modes)
    M = modes(q);
    r = new_mode();
    r.mode = q;
    r.method = 'textbook Re/Im peak picking on H_rec (Eq. 2.58-2.63)';
    r.f_n_Hz = M.f_n_Hz; r.w_n = M.w_n;
    r.f_a_Hz = M.f_a_Hz; r.f_b_Hz = M.f_b_Hz;
    r.label_n = M.label_n; r.label_a = M.label_a; r.label_b = M.label_b;
    r.peak_name = M.peak_name; r.Im_peak = M.Im_peak;
    r.zeta = M.zeta;
    r.m_modal = M.m_q; r.k_modal = M.k_q; r.c_modal = M.c_q;
    r.zeta_raw = M.zeta_raw; r.m_modal_raw = M.m_q_raw;
    r.inverted = M.inverted;
    if M.inverted, r.mass_basis = 'none (inverted mode)'; else, r.mass_basis = 'modal (FRF) and known'; end
    r.flags = M.flags;

    ilo = find(P.w >= M.win_lo, 1, 'first'); ihi = find(P.w <= M.win_hi, 1, 'last');
    hp = half_power_damping(P.w, abs(G), ilo, ihi);
    r.alt_method = 'half-power bandwidth of |H_rec|';
    r.f_n_alt_Hz = hp.w_peak/(2*pi);
    r.zeta_alt = hp.zeta;
    r.flags = [r.flags hp.flags];
    [r.f_n_alt_diff_pct, r.zeta_alt_diff_pct] = diffs(r.f_n_Hz, r.f_n_alt_Hz, r.zeta, r.zeta_alt);
    if abs(r.zeta_alt_diff_pct) > cfg.disagree_pct
        r.flags{end+1} = sprintf('zeta real-part %.4f vs half-power %.4f differ by %.1f%% (> %g%%)', r.zeta, r.zeta_alt, r.zeta_alt_diff_pct, cfg.disagree_pct);
    end
    if abs(r.f_n_alt_diff_pct) > cfg.disagree_pct
        r.flags{end+1} = sprintf('w_n Im-min vs |H| peak differ by %.1f%%', r.f_n_alt_diff_pct);
    end

    r = known_mass(r, cfg);
    if isfinite(r.m_known)
        r.m_diff_pct = 100*(r.m_modal - r.m_known)/r.m_known;
    end
    r.zeta_expwin_est = 1/(tau*r.w_n);

    inb = P.band_f_Hz >= r.f_a_Hz & P.band_f_Hz <= r.f_b_Hz;
    r.coh_masked_in_mode = sum(inb & ~P.band_used);
    if r.coh_masked_in_mode > 0
        r.flags{end+1} = sprintf('%d low-coherence point(s) removed inside [f_a, f_b]', r.coh_masked_in_mode);
    end
    recs(q) = r;

    lg(fid, cfg, '', ['  mode %d: point %d f_n = %.3f Hz (w_n = %.2f rad/s), %s = Im(H_rec) = %.4g m/N; points %d/%d: ', ...
        'f_a = %.3f, f_b = %.3f Hz'], q, r.label_n, r.f_n_Hz, r.w_n, r.peak_name, r.Im_peak, r.label_a, r.label_b, r.f_a_Hz, r.f_b_Hz);
    lg(fid, cfg, '', '          zeta_q = %.4f (Eq. 2.58) | half-power %.4f (%+.1f%%) | k_q = %.4g N/m, m_q = %.4g kg, c_q = %.4g N s/m', ...
        r.zeta, r.zeta_alt, r.zeta_alt_diff_pct, r.k_modal, r.m_modal, r.c_modal);
    lg(fid, cfg, '', '          raw bins (no refinement): zeta = %.4f, m_q = %.4g kg', r.zeta_raw, r.m_modal_raw);
    if isfinite(r.m_known)
        lg(fid, cfg, '', '          m_known = %g kg -> k = %.4g N/m, c = %.4g N s/m; m_q vs m_known: %+.1f%%', ...
            r.m_known, r.k_known, r.c_known, r.m_diff_pct);
    end
    for j = 1:numel(r.flags)
        R = warn(R, fid, cfg, sprintf('mode %d: %s', q, r.flags{j}));
    end
end
R.modes = recs;
if ~isempty(modes)
    R.K_q = diag([modes.k_q]); R.M_q = diag([modes.m_q]); R.C_q = diag([modes.c_q]);
end
% passivity: a driving-point receptance has Im(H_rec) <= 0 at every w
% (absorbed power = -w |F|^2 Im(H)/2 >= 0); an inverted mode means the FRF
% is a cross / non-collocated FRF and Eq. 2.62 does not give the physical mass
inv = cand.inverted & cand.passed;
R.passivity_ok = ~any(inv);
if any(inv)
    [~, j] = max(cand.prominence .* inv);
    R = warn(R, fid, cfg, sprintf(['PASSIVITY: Im(H_rec) has a positive (inverted) peak at %.1f Hz, %.0f%% of max|Im|. ', ...
        'A driving-point receptance has Im <= 0 at all frequencies, so this behaves like a cross / non-collocated ', ...
        'FRF: m_q, k_q, c_q are mode-shape-scaled modal values, not the physical mass/stiffness.'], ...
        cand.f_Hz(j), 100*cand.prominence(j)/max(abs(cand.Im))));
end
R.f_Hz = P.band_f_Hz; R.G = Gb; R.used = P.band_used;

if cfg.make_figures
    ttl = sprintf('%s (%s, test D): receptance, %s, polarity %+d', R.name, R.group, R.conversion, R.polarity);
    plot_fig221(P.band_f_Hz, Gb, P.band_used, modes, ttl, fullfile(figdir, [R.name '.png']), 'm/N');
end
end

% =========================================================================
function R = process_base(R, P, popts, cfg, fid, figdir)
R.frf_type_used = 'transmissibility';
R.conversion = 'none (dimensionless accel/accel)';
lg(fid, cfg, '', '  FRF type: transmissibility (from Y units "%s"); conversion applied: none. Mass cannot be identified.', R.y_units);

Tr = transmissibility_analysis(P.w, P.H, popts);
R.T = Tr;
[fs_, is_] = sort(Tr.cand_f_Hz);
keep = Tr.cand_prom(is_) >= 0.2*Tr.threshold;
lg(fid, cfg, '', '  |T| peaks (prominence threshold %.3g): %s', Tr.threshold, ...
    strjoin(arrayfun(@(f, m, p) sprintf('%.2f Hz |T|=%.3f prom=%.3f', f, m, p), fs_(keep)', Tr.cand_mag(is_(keep))', ...
    Tr.cand_prom(is_(keep))', 'UniformOutput', false), '; '));
if ~isempty(Tr.merged_f_Hz)
    lg(fid, cfg, '', '  |T| peaks merged into a stronger peak (inside its half-power band, noise split): %s Hz', ...
        sprintf('%.2f ', Tr.merged_f_Hz));
end
if isempty(Tr.peaks)
    R = warn(R, fid, cfg, 'no |T| peak passed the prominence threshold');
end

% model check: for a viscous SDOF on a moving base |T| >= 1 for all w <= sqrt(2) w_n
% (any zeta). Median |T| of the coherent points below 0.7 f_n is compared with 1.
R.T_low_median = NaN;
model_ok = true;
if ~isempty(Tr.peaks)
    fd = Tr.peaks(Tr.i_dom).f_n_Hz;
    lowsel = P.f_Hz < 0.7*fd;
    if sum(lowsel) >= 5
        R.T_low_median = median(abs(P.H(lowsel)));
        model_ok = R.T_low_median >= 0.9;
        lg(fid, cfg, '', '  model check: median |T| = %.3f for %.0f-%.0f Hz (below 0.7 f_n; base-excited SDOF requires >= 1)', ...
            R.T_low_median, P.f_Hz(find(lowsel, 1)), 0.7*fd);
        if ~model_ok
            R = warn(R, fid, cfg, sprintf(['MODEL: |T| = %.2f (median) below resonance, but a base-excited SDOF has |T| >= 1 for ', ...
                'all w <= sqrt(2) w_n. The data do not follow the base-excitation model there (input-noise bias of the H1 ', ...
                'estimate, or the input accelerometer not measuring the effective base motion?). (T-1)/w^2 cross-check skipped; ', ...
                'isolation onset is searched only above the peak.'], R.T_low_median));
        end
    else
        lg(fid, cfg, '', '  model check: fewer than 5 coherent points below 0.7 f_n, not possible');
    end
end

% cross-check per |T| peak: textbook Re/Im peak picking on G = (T - 1)/w^2,
% in a sub-band around the peak, sign taken at the peak
G = (P.H - 1)./P.w.^2;
Gb = (P.band_H - 1)./(2*pi*P.band_f_Hz).^2;
gopts = popts; gopts.n_modes = 1;
gplot = struct([]); gsign = [];
recs = new_mode();
recs = recs([]);
R.polarity = 1;
for q = 1:numel(Tr.peaks)
    pk = Tr.peaks(q);
    r = new_mode();
    r.mode = q;
    r.method = 'peak |T| (w_n) + half-power bandwidth of |T| (zeta)';
    r.f_n_Hz = pk.f_n_Hz; r.w_n = pk.w_n;
    r.zeta = pk.zeta_hp;
    r.T_peak = pk.T_peak;
    r.mass_basis = 'mass-assumed (m_known)';
    r.flags = pk.flags;
    if q == Tr.i_dom
        r.f_iso_Hz = Tr.f_iso_Hz;
        r.f_iso_over_f_n = Tr.f_iso_Hz/pk.f_n_Hz;
        if ~isfinite(Tr.f_iso_Hz)
            r.flags{end+1} = '|T| never drops below 1 above the peak inside the band';
        end
    end
    if ~isfinite(r.zeta)
        r.flags{end+1} = 'half-power zeta not available (|T| did not fall to peak/sqrt(2) inside the window)';
    end

    r.alt_method = 'textbook Re/Im peak picking on (T-1)/w^2';
    if ~model_ok
        r.alt_method = 'none (|T| < 1 below resonance: base-excitation model not met)';
        r = known_mass(r, cfg);
        recs(q) = r;
        lg(fid, cfg, '', '  mode %d: f_n = %.3f Hz (peak |T| = %.3f), zeta_hp = %.4f', q, r.f_n_Hz, r.T_peak, r.zeta);
        if q == Tr.i_dom
            lg(fid, cfg, '', '          isolation onset |T| < 1 above the peak at %.2f Hz = %.3f x f_n (viscous SDOF: sqrt(2) = 1.414)', r.f_iso_Hz, r.f_iso_over_f_n);
        end
        if isfinite(r.m_known)
            lg(fid, cfg, '', '          mass-assumed (m_known = %g kg): k = %.4g N/m, c = %.4g N s/m', r.m_known, r.k_known, r.c_known);
        end
        for j = 1:numel(r.flags), R = warn(R, fid, cfg, sprintf('mode %d: %s', q, r.flags{j})); end
        continue
    end
    sb = P.w >= pk.w_n*(1 - cfg.real_search_frac) & P.w <= pk.w_n*(1 + cfg.real_search_frac);
    [s, pinfo] = choose_polarity(P.w(sb), G(sb), cfg, pk.w_n);
    if q == Tr.i_dom, R.polarity = s; end
    [gm, gc] = peak_pick_textbook(P.w(sb), s*G(sb), gopts);
    lg(fid, cfg, '', '  mode %d: G = (T-1)/w^2 sign %+d at the |T| peak (%s); Re sign pattern consistent: %d', ...
        q, s, polarity_src(cfg), pinfo.re_consistent);
    log_candidates(fid, cfg, gc, sprintf('mode %d: Im((T-1)/w^2) extrema', q), 's^2');
    if s < 0
        r.flags{end+1} = '(T-1)/w^2 needed a sign inversion at this peak; unexpected for an accel/accel ratio';
    end
    if ~isempty(gm)
        g = gm(1);
        r.f_n_alt_Hz = g.f_n_Hz; r.zeta_alt = g.zeta; r.zeta_raw = g.zeta_raw;
        r.f_a_Hz = g.f_a_Hz; r.f_b_Hz = g.f_b_Hz;
        r.inverted = g.inverted;
        r.mq_over_m_from_T = g.m_q;
        r.flags = [r.flags cellfun(@(t) ['(T-1)/w^2: ' t], g.flags, 'UniformOutput', false)];
        [r.f_n_alt_diff_pct, r.zeta_alt_diff_pct] = diffs(r.f_n_Hz, r.f_n_alt_Hz, r.zeta, r.zeta_alt);
        if abs(r.zeta_alt_diff_pct) > cfg.disagree_pct
            r.flags{end+1} = sprintf('zeta half-power %.4f vs real-part on (T-1)/w^2 %.4f differ by %.1f%% (> %g%%)', r.zeta, r.zeta_alt, r.zeta_alt_diff_pct, cfg.disagree_pct);
        end
        if abs(r.f_n_alt_diff_pct) > cfg.disagree_pct
            r.flags{end+1} = sprintf('w_n peak|T| vs Im-min of (T-1)/w^2 differ by %.1f%%', r.f_n_alt_diff_pct);
        end
        if isempty(gplot), gplot = g; else, gplot(end+1) = g; end %#ok<AGROW>
        gsign(end+1) = s; %#ok<AGROW>
    else
        r.flags{end+1} = 'no resolved textbook bracket in (T-1)/w^2 near this peak (cross-check not available)';
    end
    r = known_mass(r, cfg);
    if isfinite(pk.hp_w1) && isfinite(pk.hp_w2)
        inb = P.band_f_Hz >= pk.hp_w1/(2*pi) & P.band_f_Hz <= pk.hp_w2/(2*pi);
        r.coh_masked_in_mode = sum(inb & ~P.band_used);
        if r.coh_masked_in_mode > 0
            r.flags{end+1} = sprintf('%d low-coherence point(s) removed inside the half-power band', r.coh_masked_in_mode);
        end
    end
    recs(q) = r;

    lg(fid, cfg, '', '  mode %d: f_n = %.3f Hz (peak |T| = %.3f), zeta_hp = %.4f | (T-1)/w^2: f_n = %.3f Hz, zeta = %.4f (%+.1f%%)', ...
        q, r.f_n_Hz, r.T_peak, r.zeta, r.f_n_alt_Hz, r.zeta_alt, r.zeta_alt_diff_pct);
    if q == Tr.i_dom
        lg(fid, cfg, '', '          isolation onset |T| < 1 above the peak at %.2f Hz = %.3f x f_n (viscous SDOF: sqrt(2) = 1.414)', r.f_iso_Hz, r.f_iso_over_f_n);
    end
    if isfinite(r.m_known)
        lg(fid, cfg, '', '          mass-assumed (m_known = %g kg): k = %.4g N/m, c = %.4g N s/m', r.m_known, r.k_known, r.c_known);
    end
    for j = 1:numel(r.flags)
        R = warn(R, fid, cfg, sprintf('mode %d: %s', q, r.flags{j}));
    end
end
R.modes = recs;
Gb = R.polarity*Gb;
R.f_Hz = P.band_f_Hz; R.G = Gb; R.used = P.band_used;

if cfg.make_figures
    % express every cross-check mode in the plotted sign and number the points 1..2n
    ng = numel(gplot);
    for q = 1:ng
        f_ = gsign(q)*R.polarity;
        gplot(q).Im_peak = f_*gplot(q).Im_peak; gplot(q).Re_a = f_*gplot(q).Re_a; gplot(q).Re_b = f_*gplot(q).Re_b;
        gplot(q).label_n = q; gplot(q).label_a = ng + 2*q - 1; gplot(q).label_b = ng + 2*q;
        gplot(q).peak_name = char('A' + q - 1);
    end
    Tinfo.f_Hz = P.band_f_Hz; Tinfo.magT = abs(P.band_H); Tinfo.used = P.band_used;
    Tinfo.peaks = Tr.peaks; Tinfo.f_iso_Hz = Tr.f_iso_Hz;
    ttl = sprintf('%s (%s, test T): |T| and G = (T-1)/w^2, sign %+d', R.name, R.group, R.polarity);
    plot_fig221(P.band_f_Hz, Gb, P.band_used, gplot, ttl, fullfile(figdir, [R.name '.png']), 's^2', Tinfo);
end
end

% =========================================================================
function X = cross_check(results, cfg, fid)
lg(fid, cfg, '', '');
lg(fid, cfg, '', '--- CROSS-CHECK: D (real-part zeta) vs T (half-power zeta), same group and specimen ---');
X = struct('group', {}, 'specimen', {}, 'T_mode', {}, 'D_mode', {}, 'f_n_D_Hz', {}, 'f_n_T_Hz', {}, ...
    'f_diff_pct', {}, 'zeta_D', {}, 'zeta_T_hp', {}, 'zeta_diff_pct', {}, 'zeta_T_realpart', {}, 'flag', {});
iT = find(strcmp({results.test}, 'T'));
for i = iT
    T = results(i);
    iD = find(strcmp({results.test}, 'D') & strcmp({results.group}, T.group) & [results.specimen] == T.specimen);
    if isempty(iD) || isempty(T.modes)
        lg(fid, cfg, 'WARN', '%s: no matching D result to compare', T.name);
        continue
    end
    D = results(iD(1));
    if isempty(D.modes), continue, end
    for q = 1:numel(T.modes)
        [~, j] = min(abs([D.modes.f_n_Hz] - T.modes(q).f_n_Hz));
        x.group = T.group; x.specimen = T.specimen; x.T_mode = q; x.D_mode = j;
        x.f_n_D_Hz = D.modes(j).f_n_Hz; x.f_n_T_Hz = T.modes(q).f_n_Hz;
        x.f_diff_pct = 100*(x.f_n_T_Hz - x.f_n_D_Hz)/x.f_n_D_Hz;
        x.zeta_D = D.modes(j).zeta; x.zeta_T_hp = T.modes(q).zeta;
        x.zeta_diff_pct = 100*(x.zeta_T_hp - x.zeta_D)/x.zeta_D;
        x.zeta_T_realpart = T.modes(q).zeta_alt;
        fl = {};
        if abs(x.f_diff_pct) > cfg.disagree_pct, fl{end+1} = sprintf('f_n differs %.1f%%', x.f_diff_pct); end %#ok<AGROW>
        if abs(x.zeta_diff_pct) > cfg.disagree_pct, fl{end+1} = sprintf('zeta differs %.1f%%', x.zeta_diff_pct); end %#ok<AGROW>
        if ~(isfinite(x.f_diff_pct) && isfinite(x.zeta_diff_pct)), fl{end+1} = 'comparison not possible (NaN)'; end %#ok<AGROW>
        x.flag = strjoin(fl, '; ');
        X(end+1) = x; %#ok<AGROW>
        lg(fid, cfg, '', '%-5s S%-2d T mode %d vs D mode %d: f_n %.2f / %.2f Hz (%+.1f%%), zeta D %.4f / T %.4f (%+.1f%%), T real-part %.4f %s', ...
            x.group, x.specimen, q, j, x.f_n_D_Hz, x.f_n_T_Hz, x.f_diff_pct, x.zeta_D, x.zeta_T_hp, x.zeta_diff_pct, ...
            x.zeta_T_realpart, ifelse_str(isempty(x.flag), '', ['<-- FLAG: ' x.flag]));
    end
end
end

% =========================================================================
function r = known_mass(r, cfg)
r.m_known = cfg.m_known_kg;
if isfinite(r.m_known)
    r.k_known = r.m_known*r.w_n^2;
    r.c_known = 2*r.zeta*r.m_known*r.w_n;
end
end

function [df, dz] = diffs(f, f_alt, z, z_alt)
df = 100*(f_alt - f)/f;
dz = 100*(z_alt - z)/z;
end

function [s, info] = choose_polarity(w, G, cfg, w_ref)
if nargin < 4, w_ref = []; end
[s, info] = detect_frf_polarity(w, G, w_ref);
if ~ischar(cfg.polarity)
    s = sign(cfg.polarity);
    if s == 0, s = 1; end
end
end

function t = polarity_src(cfg)
if ischar(cfg.polarity), t = 'auto: Im < 0 at the dominant resonance'; else, t = 'set by cfg.polarity'; end
end

function log_candidates(fid, cfg, cand, what, units)
lg(fid, cfg, '', '  %s candidates (prominence threshold %.3g %s; listed if > 20%% of it):', what, cand.threshold, units);
[~, o] = sort(cand.f_Hz);
for j = o(:)'
    if cand.passed(j) || cand.prominence(j) > 0.2*cand.threshold
        if cand.inverted(j), t = 'Im max (inverted)'; else, t = 'Im min'; end
        lg(fid, cfg, '', '      %8.2f Hz  %-17s %10.4g  prominence %10.4g  %s', cand.f_Hz(j), t, cand.Im(j), cand.prominence(j), cand.status{j});
    end
end
end

function log_acquisition(fid, cfg, C, spec)
lg(fid, cfg, '', 'acquisition (%s): FFT size %d, fs = %.4f Hz -> df = %.6f Hz, %d averages, windows: response %s / reference %s', ...
    relpath(C.file, cfg.data_root), C.FFTSize, C.SampleRate, C.SampleRate/C.FFTSize, C.NumberOfAverages, C.WindowTypeResp, C.WindowTypeRef);
for k = 1:numel(C.channels)
    ch = C.channels(k);
    if ch.enabled
        lg(fid, cfg, '', '  channel "%s": unit label %s, sensitivity %g mV/unit', ch.name, ch.unit, ch.mV_per_unit);
    end
end
if strcmp(spec.input, 'force')
    names = {C.channels([C.channels.enabled]).name};
    if any(~cellfun(@isempty, regexpi(names, 'floor')))
        lg(fid, cfg, 'WARN', ['the response channel is named "Floor Accelerometer". The textbook modal mass assumes a ', ...
            'DRIVING-POINT FRF (accelerometer on the struck mass). If the accelerometer was on the floor, m_modal is not the mass. Confirm.']);
    end
    if strcmpi(C.WindowTypeResp, 'Exponential') && C.ExpWindowFinalValue > 0 && C.ExpWindowFinalValue < 1
        tau = (C.FFTSize/C.SampleRate)/log(1/C.ExpWindowFinalValue);
        lg(fid, cfg, '', ['  exponential response window (final value %g): if it decays from 1 at record start to %g at record end ', ...
            '(QuickDAQ definition not verified), tau = %.3f s and the window adds about dzeta = 1/(tau w_n) to the measured zeta. ', ...
            'Reported per mode as zeta_expwin_est; NOT subtracted.'], C.ExpWindowFinalValue, C.ExpWindowFinalValue, tau);
    end
end
end

function log_limitations(fid, cfg)
L = {
 'Peak picking assumes a viscous SDOF model near each mode and proportional damping (textbook p. 42). The rTPU mount is viscoelastic, so k and c are EFFECTIVE values at w_n, not constants across frequency.'
 'The textbook method assumes well-separated modes; overlapping Re brackets are flagged per mode.'
 'Mass-loaded mounts can show rocking and other modes. Modes are chosen by prominence, not by order; all picked modes are reported. The first peak is not assumed to be the vertical mode.'
 'Points with coherence below min_coherence are removed (QuickDAQ coherence from the averaged strikes / records); counts are logged per file and per mode.'
 'm_modal (Eq. 2.62) requires a calibrated DRIVING-POINT FRF. It scales directly with the hammer and accelerometer sensitivities and with the accelerometer location; see the acquisition notes above.'
 'An overall sign inversion of the FRF does not change w_n, zeta, |A| or m_modal; it only decides which extremum the textbook points refer to.'
 'Frequency resolution is df = 0.39 Hz; picked points are refined by a local quadratic fit, but zeta from a Re bracket spanning few bins remains approximate.'
 'For a viscous SDOF the |T| peak lies slightly below w_n and the half-power width of |T| is approximate; they are used as a cross-check per the task definition.'
 'Transmissibility k and c are mass-assumed (k = m_known w_n^2, c = 2 zeta m_known w_n); T contains no force information.'
 };
lg(fid, cfg, '', '');
lg(fid, cfg, '', '--- KNOWN LIMITATIONS ---');
for i = 1:numel(L), lg(fid, cfg, '', '- %s', L{i}); end
end

% =========================================================================
function print_table(fid, cfg, results)
lg(fid, cfg, '', '');
lg(fid, cfg, '', '--- RESULTS (per specimen, test, mode; * = inverted mode, mass not defined) ---');
lg(fid, cfg, '', '%-10s %-4s %2s %9s %8s %10s %8s %8s %11s %11s %11s %11s', 'specimen', 'test', 'q', 'f_n[Hz]', 'zeta', ...
    'm_modal', 'm_known', '%diff', 'k_mod[N/m]', 'c_mod[Ns/m]', 'k_kn[N/m]', 'c_kn[Ns/m]');
for i = 1:numel(results)
    R = results(i);
    for q = 1:numel(R.modes)
        r = R.modes(q);
        tst = R.test; if r.inverted, tst = [tst '*']; end
        lg(fid, cfg, '', '%-10s %-4s %2d %9.3f %8.4f %10.4g %8.4g %8.2f %11.4g %11.4g %11.4g %11.4g', R.name, tst, q, r.f_n_Hz, ...
            r.zeta, r.m_modal, r.m_known, r.m_diff_pct, r.k_modal, r.c_modal, r.k_known, r.c_known);
    end
    if isempty(R.modes)
        lg(fid, cfg, '', '%-10s %-4s  -  no modes identified', R.name, R.test);
    end
end
end

function write_results_csv(fname, results)
f = fopen(fname, 'w');
hdr = {'group', 'specimen', 'test', 'file', 'frf_type_used', 'conversion', 'polarity', 'mode', 'inverted', 'method', ...
    'f_n_Hz', 'w_n_rad_s', 'zeta', 'f_a_Hz', 'f_b_Hz', 'textbook_points', 'Im_peak_name', 'Im_peak_m_per_N', ...
    'm_modal_kg', 'k_modal_N_per_m', 'c_modal_Ns_per_m', 'm_known_kg', 'm_diff_pct', 'k_known_N_per_m', 'c_known_Ns_per_m', ...
    'mass_basis', 'alt_method', 'f_n_alt_Hz', 'zeta_alt', 'f_n_alt_diff_pct', 'zeta_alt_diff_pct', 'zeta_rawbins', 'm_modal_rawbins_kg', ...
    'T_peak', 'f_iso_Hz', 'f_iso_over_f_n', 'zeta_expwin_est', 'n_band_points', 'n_coh_masked', 'coh_masked_in_mode', 'flags'};
fprintf(f, '%s\n', strjoin(hdr, ','));
for i = 1:numel(results)
    R = results(i);
    for q = 1:numel(R.modes)
        r = R.modes(q);
        if strcmp(R.test, 'D'), Ip = r.Im_peak; else, Ip = NaN; end
        if isfinite(r.label_n), pts = sprintf('%d;%d;%d', r.label_n, r.label_a, r.label_b); else, pts = ''; end
        fprintf(f, '%s,%d,%s,%s,%s,%s,%d,%d,%d,%s,', q_(R.group), R.specimen, R.test, q_(R.file), R.frf_type_used, q_(R.conversion), ...
            R.polarity, q, r.inverted, q_(r.method));
        fprintf(f, '%.6g,%.6g,%.6g,%.6g,%.6g,%s,%s,%.6g,', r.f_n_Hz, r.w_n, r.zeta, r.f_a_Hz, r.f_b_Hz, pts, r.peak_name, Ip);
        fprintf(f, '%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%s,', r.m_modal, r.k_modal, r.c_modal, r.m_known, r.m_diff_pct, r.k_known, r.c_known, q_(r.mass_basis));
        fprintf(f, '%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,', q_(r.alt_method), r.f_n_alt_Hz, r.zeta_alt, r.f_n_alt_diff_pct, r.zeta_alt_diff_pct, ...
            r.zeta_raw, r.m_modal_raw);
        fprintf(f, '%.6g,%.6g,%.6g,%.6g,%d,%d,%d,%s\n', r.T_peak, r.f_iso_Hz, r.f_iso_over_f_n, r.zeta_expwin_est, R.n_band, R.n_coh_masked, ...
            r.coh_masked_in_mode, q_(strjoin(r.flags, ' | ')));
    end
end
fclose(f);
end

function write_xcheck_csv(fname, X)
f = fopen(fname, 'w');
fprintf(f, 'group,specimen,T_mode,D_mode,f_n_D_Hz,f_n_T_Hz,f_diff_pct,zeta_D_realpart,zeta_T_halfpower,zeta_diff_pct,zeta_T_realpart,flag\n');
for i = 1:numel(X)
    x = X(i);
    fprintf(f, '%s,%d,%d,%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%s\n', q_(x.group), x.specimen, x.T_mode, x.D_mode, ...
        x.f_n_D_Hz, x.f_n_T_Hz, x.f_diff_pct, x.zeta_D, x.zeta_T_hp, x.zeta_diff_pct, x.zeta_T_realpart, q_(x.flag));
end
fclose(f);
end

function s = q_(s)
s = ['"' strrep(s, '"', '""') '"'];
end

% =========================================================================
function R = new_result()
R.name = ''; R.file = ''; R.group = ''; R.specimen = NaN; R.test = '';
R.frf_kind_header = ''; R.y_units = ''; R.frf_type_used = ''; R.conversion = ''; R.polarity = 1;
R.n_band = 0; R.n_used = 0; R.n_coh_masked = 0;
R.modes = new_mode(); R.modes = R.modes([]);
R.cand = []; R.T = [];
R.K_q = []; R.M_q = []; R.C_q = [];
R.f_Hz = []; R.G = []; R.used = [];
R.warnings = {};
R.passivity_ok = true;
R.T_low_median = NaN;
end

function r = new_mode()
r.mode = NaN; r.method = ''; r.f_n_Hz = NaN; r.w_n = NaN; r.zeta = NaN;
r.f_a_Hz = NaN; r.f_b_Hz = NaN; r.label_n = NaN; r.label_a = NaN; r.label_b = NaN;
r.peak_name = ''; r.Im_peak = NaN;
r.m_modal = NaN; r.k_modal = NaN; r.c_modal = NaN;
r.m_known = NaN; r.m_diff_pct = NaN; r.k_known = NaN; r.c_known = NaN;
r.mass_basis = '';
r.alt_method = ''; r.f_n_alt_Hz = NaN; r.zeta_alt = NaN; r.f_n_alt_diff_pct = NaN; r.zeta_alt_diff_pct = NaN;
r.T_peak = NaN; r.f_iso_Hz = NaN; r.f_iso_over_f_n = NaN; r.mq_over_m_from_T = NaN;
r.zeta_raw = NaN; r.m_modal_raw = NaN; r.inverted = false;
r.zeta_expwin_est = NaN; r.coh_masked_in_mode = 0;
r.flags = {};
end

function cfg = apply_defaults(cfg)
here = fileparts(mfilename('fullpath'));
d = struct('m_known_kg', NaN, 'frf_type', 'auto', 'n_modes', 2, 'min_prominence', [], ...
    'min_prominence_rel', 0.05, 'polarity', 'auto', 'f_band_Hz', [30 500], 'min_coherence', 0.75, ...
    'smooth_points', 1, 'smooth_points_T', 11, 'smooth_method_T', 'sgolay', 'peak_refine', 'sdof', 'real_search_frac', 0.5, 'disagree_pct', 15, ...
    'phase_units', 'auto', 'data_root', here, 'results_dir', '', 'make_figures', true, 'verbose', true);
fn = fieldnames(d);
for i = 1:numel(fn)
    if ~isfield(cfg, fn{i}), cfg.(fn{i}) = d.(fn{i}); end
end
if isempty(cfg.results_dir), cfg.results_dir = fullfile(cfg.data_root, 'results'); end
if ~isfield(cfg, 'tests') || isempty(cfg.tests)
    cfg.tests = struct('id', {'D', 'T'}, ...
        'folder', {'Disturbance Rejection Test', 'Transmissibility test'}, ...
        'input', {'force', 'base'});
end
end

function files = list_csv(folder)
% recursive *.csv listing (dir(**) is not available in Octave)
files = {};
d = dir(folder);
for i = 1:numel(d)
    if d(i).isdir
        if any(strcmp(d(i).name, {'.', '..'})), continue, end
        files = [files list_csv(fullfile(folder, d(i).name))]; %#ok<AGROW>
    elseif numel(d(i).name) > 4 && strcmpi(d(i).name(end-3:end), '.csv')
        files{end+1} = fullfile(folder, d(i).name); %#ok<AGROW>
    end
end
files = sort(files);
end

function p = relpath(p, root)
if strncmp(p, root, numel(root)), p = p(numel(root)+2:end); end
end

function R = warn(R, fid, cfg, msg)
R.warnings{end+1} = msg;
lg(fid, cfg, 'WARN', '  %s', msg);
end

function s = ifelse_str(c, a, b)
if c, s = a; else, s = b; end
end

function lg(fid, cfg, level, fmt, varargin)
msg = sprintf(fmt, varargin{:});
if ~isempty(level), msg = ['[' level '] ' msg]; end
fprintf(fid, '%s\n', msg);
if cfg.verbose, fprintf('%s\n', msg); end
end
