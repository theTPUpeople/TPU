function R = process_modal_file(file, cfg)
%PROCESS_MODAL_FILE Peak-picking modal fit of one FRF export (Schmitz & Smith, Sec. 2.5.1).
%   R = PROCESS_MODAL_FILE(FILE, CFG) runs the same steps on any QuickDAQ CSV
%   export, disturbance rejection (accelerance) or transmissibility:
%     1. read frequency, FRF magnitude, FRF phase [deg] and coherence; the
%        columns are found from the "Measurement Type" header line
%     2. phase [deg] -> [rad]; Re = A_mag cos(phase), Im = A_mag sin(phase)
%     3. accelerance only: H_rec = H/(-omega^2), because Eqs. 2.58-2.63 assume
%        displacement/force (a transmissibility has no force and is kept as is)
%     4. sign convention of Fig. 2.21 (Im < 0 at resonance), CFG.polarity
%     5. quadratic Savitzky-Golay smoothing over CFG.pick_smooth_points bins
%     6. points 1-6 of Fig. 2.21 inside CFG.f_band_Hz, CFG.pick_n_modes modes
%     7. Eqs. 2.58-2.63 and the modal matrices K_q, M_q, C_q (p. 43)
%     8. Re and Im plots with the points, and a one-page PDF
%        "<file name> processed.pdf" in CFG.report_dir
%   R holds every value; the fields are named after the textbook symbols.
%
%   CFG fields used: f_band_Hz, polarity, frf_type, min_coherence,
%   min_prominence_rel, pick_n_modes, pick_smooth_points, pick_show_plots,
%   report_dir.

[~, name] = fileparts(file);
notes = {};

% ---------------- 1. read the export ----------------
L = regexp(fileread(file), '\r?\n', 'split');
r_mt  = find(strncmp(L, 'Measurement Type', 16), 1);
r_yu  = find(strncmp(L, 'Y Axis Units', 12), 1);
r_lab = find(strncmp(L, 'Frequency,', 10), 1);
if isempty(r_mt) || isempty(r_yu) || isempty(r_lab)
    error('process_modal_file:header', '%s: QuickDAQ header lines not found.', file);
end
mt  = strtrim(regexp(L{r_mt}, ',', 'split'));
lab = strtrim(regexp(L{r_lab}, ',', 'split'));
yu  = strtrim(regexp(L{r_yu}, ',', 'split'));
col_frf = find(strcmpi(mt, 'FRF'), 1);                 % FRF magnitude; the phase follows it
col_coh = find(strcmpi(mt, 'Coherence'), 1);
if isempty(col_frf) || numel(lab) < col_frf + 1 || ~strcmpi(lab{col_frf}, 'Magnitude') || ~strcmpi(lab{col_frf + 1}, 'Phase')
    error('process_modal_file:layout', '%s: FRF columns are not "Magnitude, Phase" as expected.', file);
end
y_units = yu{col_frf};
data = dlmread(file, ',', r_lab, 0);                   % numeric block below the header
freq      = data(:, 1);                                % frequency [Hz]
A_mag     = data(:, col_frf);                          % FRF magnitude (the "A" of Re = A cos)
phase_deg = data(:, col_frf + 1);                      % FRF phase [deg]
if isempty(col_coh), coh = nan(size(freq)); else, coh = data(:, col_coh); end

% ---------------- 2. real and imaginary components ----------------
phase_rad = phase_deg*pi/180;                          % phase [rad]
Re = A_mag.*cos(phase_rad);                            % Re = A cos(phase)
Im = A_mag.*sin(phase_rad);                            % Im = A sin(phase)

% ---------------- 3. FRF type and conversion ----------------
switch lower(regexprep(y_units, '\s', ''))
    case 'm/s^2/n',     frf_type = 'accelerance';
    case 'm/n',         frf_type = 'receptance';
    case 'm/s^2/m/s^2', frf_type = 'transmissibility';
    otherwise
        error('process_modal_file:units', '%s: Y axis units "%s" not recognised; set cfg.frf_type.', file, y_units);
end
if ~strcmp(frf_type, 'transmissibility') && isfield(cfg, 'frf_type') && any(strcmpi(cfg.frf_type, {'accelerance', 'receptance'}))
    frf_type = lower(cfg.frf_type);                    % user override for the force-input test
end
keep = freq > 0;                                       % omega = 0 cannot be divided by
freq = freq(keep); Re = Re(keep); Im = Im(keep); coh = coh(keep);
omega = 2*pi*freq;                                     % [rad/s]
switch frf_type
    case 'accelerance'
        H = (Re + 1i*Im)./(-omega.^2);                 % receptance H_rec = H_acc/(-omega^2)  [m/N]
        Re = real(H); Im = imag(H);
        conversion = 'accelerance [m/s^2/N] -> receptance H/(-omega^2) [m/N]';
        u = struct('frf', 'm/N', 'k', 'N/m', 'm', 'kg', 'c', 'N s/m');
    case 'receptance'
        conversion = 'none (receptance [m/N])';
        u = struct('frf', 'm/N', 'k', 'N/m', 'm', 'kg', 'c', 'N s/m');
    case 'transmissibility'
        conversion = 'none (transmissibility, accel/accel, dimensionless)';
        u = struct('frf', '-', 'k', '-', 'm', 's^2', 'c', 's');
        notes{end+1} = ['Transmissibility has no force: A, B are dimensionless and k_q, m_q, c_q come out in ' ...
            '[-], [s^2], [s], not N/m, kg, N s/m. For one mode on a moving base Im(T) = -1/(2 zeta) at ' ...
            'resonance, so Eq. 2.60 gives k_q1 close to 1.'];
end

% ---------------- 4. sign convention of Fig. 2.21 ----------------
inb = find(freq >= cfg.f_band_Hz(1) & freq <= cfg.f_band_Hz(2));
if ischar(cfg.polarity)
    polarity = detect_frf_polarity(omega(inb), Re(inb) + 1i*Im(inb));
else
    polarity = sign(cfg.polarity);
end
Re = polarity*Re; Im = polarity*Im;
if polarity < 0
    notes{end+1} = 'Sign inverted (x -1) to the Fig. 2.21 convention: the FRF had Im > 0 at resonance.';
end

% ---------------- 5. smoothing (quadratic Savitzky-Golay) ----------------
m = floor(cfg.pick_smooth_points/2);
sg = (3*(3*m^2 + 3*m - 1) - 15*(-m:m)'.^2)/((2*m + 1)*(4*m^2 + 4*m - 3));
Re_s = conv(Re, sg, 'same'); Re_s([1:m, end-m+1:end]) = Re([1:m, end-m+1:end]);   % edge bins unsmoothed
Im_s = conv(Im, sg, 'same'); Im_s([1:m, end-m+1:end]) = Im([1:m, end-m+1:end]);

% ---------------- 6. points 1-6 of Fig. 2.21 ----------------
% point 1, point 2   minima of the imaginary part -> omega_n1, omega_n2 and the peak values A, B
% points 3, 4        max / min of the real part around mode 1 -> omega_3, omega_4
% points 5, 6        max / min of the real part around mode 2 -> omega_5, omega_6
[~, k] = min(Im_s(inb));  i_n = inb(k);               % deepest minimum of the imaginary part
if cfg.pick_n_modes >= 2
    % second mode: most prominent other minimum of the imaginary part outside the first mode's
    % real-part bracket (a minimum inside that bracket belongs to the same mode)
    seg = inb(inb <= i_n); [~, k] = max(Re_s(seg)); i_lo = seg(k);
    seg = inb(inb >= i_n); [~, k] = min(Re_s(seg)); i_hi = seg(k);
    [ip, prom] = find_peaks_prominence(-Im_s(inb));
    cand = inb(ip);
    ok = Im_s(cand) < 0 & prom >= cfg.min_prominence_rel*max(abs(Im_s(inb))) & (cand < i_lo | cand > i_hi);
    cand = cand(ok); prom = prom(ok);
    if ~isempty(cand)
        [~, k] = max(prom);
        i_n = sort([i_n; cand(k)]);                    % mode 1 = lower frequency
    else
        notes{end+1} = 'No second mode found outside the first mode''s bracket: mode 2 values are NaN.';
    end
end
has_mode2 = numel(i_n) == 2;
if has_mode2, i_mid = round(mean(i_n)); else, i_mid = inb(end); end   % boundary between the modes
seg = inb(inb <= i_n(1));                [~, k] = max(Re_s(seg)); i_3 = seg(k);    % point 3
seg = inb(inb >= i_n(1) & inb <= i_mid); [~, k] = min(Re_s(seg)); i_4 = seg(k);    % point 4
omega_n1 = omega(i_n(1));  A = Im_s(i_n(1));          % point 1
omega_3 = omega(i_3);      omega_4 = omega(i_4);
if has_mode2
    seg = inb(inb >= i_mid & inb <= i_n(2)); [~, k] = max(Re_s(seg)); i_5 = seg(k); % point 5
    seg = inb(inb >= i_n(2));                [~, k] = min(Re_s(seg)); i_6 = seg(k); % point 6
    omega_n2 = omega(i_n(2));  B = Im_s(i_n(2));      % point 2
    omega_5 = omega(i_5);      omega_6 = omega(i_6);
    idx = [i_n(1) i_n(2) i_3 i_4 i_5 i_6];
else
    omega_n2 = NaN; B = NaN; omega_5 = NaN; omega_6 = NaN;
    idx = [i_n(1) NaN i_3 i_4 NaN NaN];
    if cfg.pick_n_modes < 2
        notes{end+1} = 'One mode fitted (cfg.pick_n_modes = 1): Eqs. 2.59 and 2.61 and the mode 2 parts of 2.62-2.63 are NaN.';
    end
end

% ---------------- 7. Eqs. 2.58-2.63 ----------------
% Eq. 2.58: omega_4 - omega_3 = omega_n1(1 + zeta_q1) - omega_n1(1 - zeta_q1) = 2 zeta_q1 omega_n1
zeta_q1 = (omega_4 - omega_3)/(2*omega_n1);
% Eq. 2.59
zeta_q2 = (omega_6 - omega_5)/(2*omega_n2);
% Eq. 2.60: A = -1/(2 k_q1 zeta_q1), so k_q1 = -1/(2 zeta_q1 A)
k_q1 = -1/(2*zeta_q1*A);
% Eq. 2.61
k_q2 = -1/(2*zeta_q2*B);
% Eq. 2.62: omega_n1 = sqrt(k_q1/m_q1), so m_q1 = k_q1/omega_n1^2 and m_q2 = k_q2/omega_n2^2
m_q1 = k_q1/omega_n1^2;
m_q2 = k_q2/omega_n2^2;
% Eq. 2.63: zeta_q1 = c_q1/(2 sqrt(k_q1 m_q1)), so c_q1 = 2 zeta_q1 sqrt(k_q1 m_q1) and c_q2 = 2 zeta_q2 sqrt(k_q2 m_q2)
c_q1 = 2*zeta_q1*sqrt(k_q1*m_q1);
c_q2 = 2*zeta_q2*sqrt(k_q2*m_q2);
% modal matrices [K_q], [M_q], [C_q] (p. 43), one diagonal entry per fitted mode
n_q = 1 + has_mode2;
k_q_list = [k_q1 k_q2]; m_q_list = [m_q1 m_q2]; c_q_list = [c_q1 c_q2];
K_q = diag(k_q_list(1:n_q));
M_q = diag(m_q_list(1:n_q));
C_q = diag(c_q_list(1:n_q));

% data-quality notes
coh_pts = nan(1, 6);
for j = 1:6
    if isfinite(idx(j)), coh_pts(j) = coh(idx(j)); end
end
low = find(coh_pts < cfg.min_coherence);
if ~isempty(low)
    notes{end+1} = sprintf('Coherence below %.2f at point(s) %s (%s): those values are uncertain.', cfg.min_coherence, ...
        strjoin(arrayfun(@num2str, low, 'UniformOutput', false), ', '), ...
        strjoin(arrayfun(@(v) sprintf('%.2f', v), coh_pts(low), 'UniformOutput', false), ', '));
end
d_lo = omega_n1 - omega_3; d_hi = omega_4 - omega_n1;
if max(d_lo, d_hi) > 2*min(d_lo, d_hi)
    notes{end+1} = sprintf(['Points 3 and 4 are not symmetric about point 1 (%.2f Hz below, %.2f Hz above); for one ' ...
        'viscous mode they lie near omega_n1(1 -/+ zeta_q1).'], d_lo/(2*pi), d_hi/(2*pi));
end

% ---------------- results ----------------
R.file = file; R.name = name; R.frf_type = frf_type; R.y_units = y_units; R.conversion = conversion;
R.units = u; R.polarity = polarity; R.f_band_Hz = cfg.f_band_Hz; R.smooth_points = 2*m + 1;
R.n_modes = n_q; R.n_points = numel(freq); R.f_range_Hz = [freq(1) freq(end)];
R.omega_n1 = omega_n1; R.omega_n2 = omega_n2; R.omega_3 = omega_3; R.omega_4 = omega_4;
R.omega_5 = omega_5; R.omega_6 = omega_6; R.A = A; R.B = B;
R.zeta_q1 = zeta_q1; R.zeta_q2 = zeta_q2; R.k_q1 = k_q1; R.k_q2 = k_q2;
R.m_q1 = m_q1; R.m_q2 = m_q2; R.c_q1 = c_q1; R.c_q2 = c_q2;
R.K_q = K_q; R.M_q = M_q; R.C_q = C_q;
R.point_index = idx; R.point_coherence = coh_pts;
R.point_Re = nan(1, 6); R.point_Im = nan(1, 6);
for j = 1:6
    if isfinite(idx(j)), R.point_Re(j) = Re_s(idx(j)); R.point_Im(j) = Im_s(idx(j)); end
end
R.notes = notes;

% console summary
fprintf('\n%s  (%s; %s; sign %+d; %d-point smoothing; %g-%g Hz)\n', name, frf_type, conversion, polarity, 2*m + 1, cfg.f_band_Hz);
fprintf('  point 1: omega_n1 = %9.2f rad/s (%8.3f Hz)   A  = %11.4g [%s]\n', omega_n1, omega_n1/(2*pi), A, u.frf);
fprintf('  point 3: omega_3  = %9.2f rad/s (%8.3f Hz)   Re = %11.4g [%s]\n', omega_3, omega_3/(2*pi), Re_s(i_3), u.frf);
fprintf('  point 4: omega_4  = %9.2f rad/s (%8.3f Hz)   Re = %11.4g [%s]\n', omega_4, omega_4/(2*pi), Re_s(i_4), u.frf);
if has_mode2
    fprintf('  point 2: omega_n2 = %9.2f rad/s (%8.3f Hz)   B  = %11.4g [%s]\n', omega_n2, omega_n2/(2*pi), B, u.frf);
    fprintf('  point 5: omega_5  = %9.2f rad/s (%8.3f Hz)   Re = %11.4g [%s]\n', omega_5, omega_5/(2*pi), Re_s(i_5), u.frf);
    fprintf('  point 6: omega_6  = %9.2f rad/s (%8.3f Hz)   Re = %11.4g [%s]\n', omega_6, omega_6/(2*pi), Re_s(i_6), u.frf);
end
fprintf('  Eq. 2.58 zeta_q1 = %.4f   Eq. 2.60 k_q1 = %.4g [%s]   Eq. 2.62 m_q1 = %.4g [%s]   Eq. 2.63 c_q1 = %.4g [%s]\n', ...
    zeta_q1, k_q1, u.k, m_q1, u.m, c_q1, u.c);
fprintf('  Eq. 2.59 zeta_q2 = %.4f   Eq. 2.61 k_q2 = %.4g [%s]   Eq. 2.62 m_q2 = %.4g [%s]   Eq. 2.63 c_q2 = %.4g [%s]\n', ...
    zeta_q2, k_q2, u.k, m_q2, u.m, c_q2, u.c);
for j = 1:numel(notes), fprintf('  note: %s\n', notes{j}); end

% ---------------- 8. plots and PDF report ----------------
C.freq = freq; C.Re = Re_s; C.Im = Im_s; C.idx = idx; C.has_mode2 = has_mode2; C.u = u;
C.smooth = 2*m + 1; C.name = name;
if isfield(cfg, 'pick_show_plots') && cfg.pick_show_plots
    for part = 1:2
        fh = figure('Color', 'w', 'Name', sprintf('%s %s', name, ifelse(part == 1, 'real part', 'imaginary part')));
        draw_part(axes('Parent', fh), C, part, 12);
    end
end
if exist(cfg.report_dir, 'dir') ~= 7, mkdir(cfg.report_dir); end
R.pdf_file = fullfile(cfg.report_dir, [name ' processed.pdf']);
write_report(R, C);
fprintf('  report: %s\n', R.pdf_file);
end

% =========================================================================
function draw_part(ax, C, part, fs)
% Real (part 1) or imaginary (part 2) component with the Fig. 2.21 points.
line_col = [42 120 214]/255; pt_col = [235 104 52]/255; grey = [0.55 0.55 0.55];
f = C.freq;
if part == 1, y = C.Re; else, y = C.Im; end
plot(ax, f, y, '-', 'Color', line_col, 'LineWidth', 1.5);
hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
plot(ax, [f(1) f(end)], [0 0], ':', 'Color', grey);
set(ax, 'FontSize', fs);
xlim(ax, [f(1) f(end)]);
xlabel(ax, 'Frequency [Hz]');
lbl = {'Real', 'Imaginary'};
ylabel(ax, sprintf('%s part [%s]', lbl{part}, C.u.frf));
title(ax, sprintf('%s: %s part (smoothed: %d-point Savitzky-Golay)', C.name, lower(lbl{part}), C.smooth), 'Interpreter', 'none');
dx = 0.02*(f(end) - f(1));
if part == 1, pts = [3 4 5 6]; else, pts = [1 2]; end
names = {'A', 'B'};
for p = pts
    i = C.idx(p);
    if ~isfinite(i), continue, end
    if part == 2
        plot(ax, [f(1) f(i)], [y(i) y(i)], '--', 'Color', grey);
        txt = sprintf('%d: f_{n%d} = %.2f Hz,  %s = %.4g', p, p, f(i), names{p}, y(i));
    else
        txt = sprintf('%d: %.2f Hz,  Re = %.4g', p, f(i), y(i));
    end
    plot(ax, f(i), y(i), 'o', 'Color', pt_col, 'MarkerFaceColor', pt_col, 'MarkerSize', 7);
    text(f(i) + dx, y(i), txt, 'Parent', ax, 'FontSize', fs - 1, 'VerticalAlignment', 'middle');
end
end

% =========================================================================
function write_report(R, C)
% One-page PDF (8.5 x 14 in): test information, picked points, modal parameters,
% notes and the two plots.
W = 8.5; Hh = 14;
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0 0 W Hh]);
set(fig, 'PaperUnits', 'inches', 'PaperSize', [W Hh], 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 W Hh]);
ax = axes('Position', [0 0 1 1], 'Visible', 'off'); xlim(ax, [0 1]); ylim(ax, [0 1]); hold(ax, 'on');
x0 = 0.06; fs = 9; dy = 0.0125;
u = R.units;
T = @(x, y, s, varargin) text(x, y, s, 'Parent', ax, 'FontSize', fs, 'Interpreter', 'none', 'VerticalAlignment', 'middle', varargin{:});
rule = @(y) plot(ax, [x0 1 - x0], [y y], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);

text(x0, 0.978, [R.name ' processed'], 'Parent', ax, 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');
T(x0, 0.962, 'Peak-picking modal fit after Schmitz & Smith, Machining Dynamics, 2nd ed., Sec. 2.5.1 (Fig. 2.21, Eqs. 2.58-2.63)', ...
    'Color', [0.32 0.32 0.31]);

% 1. test
y = 0.940;
text(x0, y, '1. Test', 'Parent', ax, 'FontSize', 11, 'FontWeight', 'bold'); y = y - 1.3*dy;
info = {
    'File', R.file
    'Experiment / FRF type', sprintf('%s (Y axis units in file: %s)', R.frf_type, R.y_units)
    'Conversion applied', R.conversion
    'Sign convention', sprintf('%+d (Fig. 2.21: Im < 0 at resonance)', R.polarity)
    'Data', sprintf('%d points, %.3f-%.3f Hz', R.n_points, R.f_range_Hz)
    'Search band', sprintf('%g-%g Hz', R.f_band_Hz)
    'Smoothing', sprintf('%d-point quadratic Savitzky-Golay (curves and picking)', R.smooth_points)
    'Modes fitted', sprintf('%d', R.n_modes)
    'Processed', datestr(now, 'yyyy-mm-dd HH:MM')};
for j = 1:size(info, 1)
    T(x0, y, info{j, 1}, 'FontWeight', 'bold'); T(x0 + 0.20, y, info{j, 2}); y = y - dy;
end

% 2. picked points
y = y - 0.6*dy;
text(x0, y, '2. Picked points (Fig. 2.21)', 'Parent', ax, 'FontSize', 11, 'FontWeight', 'bold'); y = y - 1.3*dy;
cx = x0 + [0 0.05 0.30 0.43 0.55 0.67 0.82];
hdr = {'Point', 'Meaning', 'Variable', 'f [Hz]', 'omega [rad/s]', sprintf('Value [%s]', u.frf), 'Coherence'};
for c = 1:numel(hdr), T(cx(c), y, hdr{c}, 'FontWeight', 'bold'); end
rule(y - dy/2); y = y - dy;
pdesc = {'Im minimum, mode 1', 'Im minimum, mode 2', 'Re maximum, mode 1', 'Re minimum, mode 1', 'Re maximum, mode 2', 'Re minimum, mode 2'};
pvar = {'omega_n1, A', 'omega_n2, B', 'omega_3', 'omega_4', 'omega_5', 'omega_6'};
pom = [R.omega_n1 R.omega_n2 R.omega_3 R.omega_4 R.omega_5 R.omega_6];
pval = [R.point_Im(1:2) R.point_Re(3:6)];
for p = [1 3 4 2 5 6]
    if ~isfinite(pom(p)) && p ~= 2, continue, end
    T(cx(1), y, sprintf('%d', p)); T(cx(2), y, pdesc{p}); T(cx(3), y, pvar{p});
    T(cx(4), y, num(pom(p)/(2*pi), '%.3f')); T(cx(5), y, num(pom(p), '%.2f'));
    T(cx(6), y, num(pval(p), '%.4g')); T(cx(7), y, num(R.point_coherence(p), '%.3f'));
    y = y - dy;
end

% 3. modal parameters
y = y - 0.6*dy;
text(x0, y, '3. Modal parameters (Eqs. 2.58-2.63)', 'Parent', ax, 'FontSize', 11, 'FontWeight', 'bold'); y = y - 1.3*dy;
cx = x0 + [0 0.12 0.37 0.58 0.74];
hdr = {'Equation', 'Quantity', 'Variables (mode 1 / mode 2)', 'Mode 1', 'Mode 2'};
for c = 1:numel(hdr), T(cx(c), y, hdr{c}, 'FontWeight', 'bold'); end
rule(y - dy/2); y = y - dy;
rows = {
    '2.58 / 2.59', 'damping ratio [-]',            'zeta_q1 / zeta_q2', R.zeta_q1, R.zeta_q2, '%.4f'
    '2.60 / 2.61', sprintf('stiffness [%s]', u.k), 'k_q1 / k_q2',       R.k_q1,    R.k_q2,    '%.4g'
    '2.62',        sprintf('mass [%s]', u.m),      'm_q1 / m_q2',       R.m_q1,    R.m_q2,    '%.4g'
    '2.63',        sprintf('damping [%s]', u.c),   'c_q1 / c_q2',       R.c_q1,    R.c_q2,    '%.4g'};
for j = 1:size(rows, 1)
    T(cx(1), y, rows{j, 1}); T(cx(2), y, rows{j, 2}); T(cx(3), y, rows{j, 3});
    T(cx(4), y, num(rows{j, 4}, rows{j, 6})); T(cx(5), y, num(rows{j, 5}, rows{j, 6}));
    y = y - dy;
end
y = y - 0.4*dy;
T(x0, y, sprintf('Modal matrices (p. 43):  [K_q] = diag(%s) [%s]   [M_q] = diag(%s) [%s]   [C_q] = diag(%s) [%s]', ...
    vec(diag(R.K_q)), u.k, vec(diag(R.M_q)), u.m, vec(diag(R.C_q)), u.c));
y = y - dy;

% 4. notes
y = y - 0.6*dy;
text(x0, y, '4. Notes', 'Parent', ax, 'FontSize', 11, 'FontWeight', 'bold'); y = y - 1.3*dy;
if isempty(R.notes), T(x0, y, '-  none'); y = y - dy; end
for j = 1:numel(R.notes)
    lines = wrap(R.notes{j}, 130);
    for l = 1:numel(lines)
        if l == 1, T(x0, y, ['-  ' lines{l}]); else, T(x0, y, ['   ' lines{l}]); end
        y = y - dy;
    end
end

% 5. plots
y = y - 0.6*dy;
text(x0, y, '5. Real and imaginary parts with the picked points', 'Parent', ax, 'FontSize', 11, 'FontWeight', 'bold');
top = y - 0.03; h = (top - 0.05 - 0.05)/2;
a1 = axes('Position', [0.11 top - h 0.84 h]);       draw_part(a1, C, 1, 9);
a2 = axes('Position', [0.11 0.045 0.84 h]);         draw_part(a2, C, 2, 9);

print(fig, R.pdf_file, '-dpdf');
close(fig);
end

function s = num(v, fmt)
if ~isfinite(v), s = 'NaN'; else, s = sprintf(fmt, v); end
end

function s = vec(v)
s = strjoin(arrayfun(@(x) sprintf('%.4g', x), v(:)', 'UniformOutput', false), ', ');
end

function out = wrap(s, n)
% split a sentence into lines of at most n characters at spaces
out = {};
while numel(s) > n
    k = find(s(1:n) == ' ', 1, 'last');
    if isempty(k), k = n; end
    out{end+1} = s(1:k-1); %#ok<AGROW>
    s = s(k+1:end);
end
out{end+1} = s;
end

function v = ifelse(c, a, b)
if c, v = a; else, v = b; end
end
