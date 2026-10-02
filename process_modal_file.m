function R = process_modal_file(file, cfg)
%PROCESS_MODAL_FILE Peak-picking modal fit of one SDOF FRF export (Schmitz & Smith, Sec. 2.5.1).
%   R = PROCESS_MODAL_FILE(FILE, CFG) runs the same steps on any QuickDAQ CSV
%   export, disturbance rejection (accelerance) or transmissibility:
%     1. read frequency, FRF magnitude, FRF phase [deg] and coherence; the
%        columns are found from the "Measurement Type" header line
%     2. phase [deg] -> [rad]; Re = A_mag cos(phase), Im = A_mag sin(phase)
%     3. accelerance only: H_rec = H/(-omega^2), because Eqs. 2.58-2.63 assume
%        displacement/force (a transmissibility has no force and is kept as is)
%     4. sign convention of Fig. 2.21 (Im < 0 at resonance), CFG.polarity
%     5. quadratic Savitzky-Golay smoothing over CFG.pick_smooth_points bins
%     6. natural frequency from the raw data: maximum of the FRF magnitude as
%        read from the file (no smoothing, no conversion)
%     7. SDOF points of Fig. 2.21 inside CFG.f_band_Hz: point 1 = minimum of the
%        imaginary part (omega_n1, A), point 3 = biggest maximum of the real
%        part (omega_3), point 4 = biggest minimum of the real part (omega_4)
%     8. Eqs. 2.58, 2.60, 2.62, 2.63 (Eqs. 2.59 and 2.61 are the same equations
%        for a second mode and do not apply to SDOF data)
%     9. plots (raw magnitude, real and imaginary parts) and a one-page PDF
%        "<file name> processed.pdf" in CFG.report_dir
%   R holds every value; the fields are named after the textbook symbols.
%
%   CFG fields used: f_band_Hz, polarity, frf_type, min_coherence,
%   pick_smooth_points, pick_show_plots, report_dir.

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
freq = freq(keep); A_mag = A_mag(keep); Re = Re(keep); Im = Im(keep); coh = coh(keep);
omega = 2*pi*freq;                                     % [rad/s]
switch frf_type
    case 'accelerance'
        H = (Re + 1i*Im)./(-omega.^2);                 % receptance H_rec = H_acc/(-omega^2)  [m/N]
        Re = real(H); Im = imag(H);
        conversion = 'accelerance [m/s^2/N] -> receptance H/(-omega^2) [m/N]';
        u = struct('frf', 'm/N', 'k', 'N/m', 'm', 'kg', 'c', 'N s/m', 'raw', 'm/s^2/N');
    case 'receptance'
        conversion = 'none (receptance [m/N])';
        u = struct('frf', 'm/N', 'k', 'N/m', 'm', 'kg', 'c', 'N s/m', 'raw', 'm/N');
    case 'transmissibility'
        conversion = 'none (transmissibility, accel/accel, dimensionless)';
        u = struct('frf', '-', 'k', '-', 'm', 's^2', 'c', 's', 'raw', '-');
        notes{end+1} = ['Transmissibility has no force: A is dimensionless and k_q1, m_q1, c_q1 come out in ' ...
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

% ---------------- 6. natural frequency from the raw data ----------------
% maximum of the FRF magnitude as read from the file (no smoothing, no conversion), inside the band
[A_mag_max, k] = max(A_mag(inb));  i_raw = inb(k);
f_n_raw = freq(i_raw);                                 % natural frequency from the raw data [Hz]
omega_n_raw = 2*pi*f_n_raw;                            % [rad/s]

% ---------------- 7. SDOF points of Fig. 2.21 ----------------
[A, k] = min(Im_s(inb));     i_1 = inb(k);            % point 1: minimum of the imaginary part and its value A
omega_n1 = omega(i_1);                                 % natural frequency omega_n1 [rad/s]
[Re_3, k] = max(Re_s(inb));  i_3 = inb(k);            % point 3: biggest maximum of the real part
omega_3 = omega(i_3);
[Re_4, k] = min(Re_s(inb));  i_4 = inb(k);            % point 4: biggest minimum of the real part
omega_4 = omega(i_4);

% ---------------- 8. Eqs. 2.58, 2.60, 2.62, 2.63 (one mode) ----------------
% Eq. 2.58: omega_4 - omega_3 = omega_n1(1 + zeta_q1) - omega_n1(1 - zeta_q1) = 2 zeta_q1 omega_n1
zeta_q1 = (omega_4 - omega_3)/(2*omega_n1);
% Eq. 2.60: A = -1/(2 k_q1 zeta_q1), so k_q1 = -1/(2 zeta_q1 A)
k_q1 = -1/(2*zeta_q1*A);
% Eq. 2.62: omega_n1 = sqrt(k_q1/m_q1), so m_q1 = k_q1/omega_n1^2
m_q1 = k_q1/omega_n1^2;
% Eq. 2.63: zeta_q1 = c_q1/(2 sqrt(k_q1 m_q1)), so c_q1 = 2 zeta_q1 sqrt(k_q1 m_q1)
c_q1 = 2*zeta_q1*sqrt(k_q1*m_q1);
% (Eqs. 2.59 and 2.61 repeat 2.58 and 2.60 for a second mode with omega_5, omega_6, B; not used for SDOF data)

% data-quality notes
idx = [i_raw i_1 i_3 i_4];                             % raw maximum, points 1, 3, 4
coh_pts = coh(idx)';
lab_pts = {'raw maximum', '1', '3', '4'};
low = find(coh_pts < cfg.min_coherence);
if ~isempty(low)
    notes{end+1} = sprintf('Coherence below %.2f at %s (%s): those values are uncertain.', cfg.min_coherence, ...
        strjoin(lab_pts(low), ', '), strjoin(arrayfun(@(v) sprintf('%.2f', v), coh_pts(low), 'UniformOutput', false), ', '));
end
if omega_3 >= omega_4
    notes{end+1} = 'The biggest Re maximum lies above the biggest Re minimum, unlike Fig. 2.21: zeta_q1 is not valid.';
else
    d_lo = omega_n1 - omega_3; d_hi = omega_4 - omega_n1;
    if d_lo <= 0 || d_hi <= 0 || max(d_lo, d_hi) > 2*min(d_lo, d_hi)
        notes{end+1} = sprintf(['Points 3 and 4 are not symmetric about point 1 (%.2f Hz below, %.2f Hz above); for one ' ...
            'viscous mode they lie near omega_n1(1 -/+ zeta_q1).'], d_lo/(2*pi), d_hi/(2*pi));
    end
end

% ---------------- results ----------------
R.file = file; R.name = name; R.frf_type = frf_type; R.y_units = y_units; R.conversion = conversion;
R.units = u; R.polarity = polarity; R.f_band_Hz = cfg.f_band_Hz; R.smooth_points = 2*m + 1;
R.n_points = numel(freq); R.f_range_Hz = [freq(1) freq(end)];
R.f_n_raw = f_n_raw; R.omega_n_raw = omega_n_raw; R.A_mag_max = A_mag_max;
R.omega_n1 = omega_n1; R.A = A; R.omega_3 = omega_3; R.Re_3 = Re_3; R.omega_4 = omega_4; R.Re_4 = Re_4;
R.zeta_q1 = zeta_q1; R.k_q1 = k_q1; R.m_q1 = m_q1; R.c_q1 = c_q1;
R.point_index = idx; R.point_coherence = coh_pts;
R.notes = notes;

% console summary
fprintf('\n%s  (%s; %s; sign %+d; %d-point smoothing; %g-%g Hz)\n', name, frf_type, conversion, polarity, 2*m + 1, cfg.f_band_Hz);
fprintf('  raw maximum: omega_n_raw = %9.2f rad/s (%8.3f Hz)   |FRF| = %11.4g [%s]\n', omega_n_raw, f_n_raw, A_mag_max, u.raw);
fprintf('  point 1:     omega_n1    = %9.2f rad/s (%8.3f Hz)   A     = %11.4g [%s]\n', omega_n1, omega_n1/(2*pi), A, u.frf);
fprintf('  point 3:     omega_3     = %9.2f rad/s (%8.3f Hz)   Re    = %11.4g [%s]\n', omega_3, omega_3/(2*pi), Re_3, u.frf);
fprintf('  point 4:     omega_4     = %9.2f rad/s (%8.3f Hz)   Re    = %11.4g [%s]\n', omega_4, omega_4/(2*pi), Re_4, u.frf);
fprintf('  Eq. 2.58 zeta_q1 = %.4f   Eq. 2.60 k_q1 = %.4g [%s]   Eq. 2.62 m_q1 = %.4g [%s]   Eq. 2.63 c_q1 = %.4g [%s]\n', ...
    zeta_q1, k_q1, u.k, m_q1, u.m, c_q1, u.c);
for j = 1:numel(notes), fprintf('  note: %s\n', notes{j}); end

% ---------------- 9. plots and PDF report ----------------
C.freq = freq; C.A_mag = A_mag; C.Re = Re_s; C.Im = Im_s; C.idx = idx; C.u = u;
C.smooth = 2*m + 1; C.name = name;
if isfield(cfg, 'pick_show_plots') && cfg.pick_show_plots
    for part = 0:2
        fh = figure('Color', 'w', 'Name', sprintf('%s %s', name, part_name(part)));
        draw_part(axes('Parent', fh), C, part, 12);
    end
end
if exist(cfg.report_dir, 'dir') ~= 7, mkdir(cfg.report_dir); end
R.pdf_file = fullfile(cfg.report_dir, [name ' processed.pdf']);
write_report(R, C);
fprintf('  report: %s\n', R.pdf_file);
end

% =========================================================================
function s = part_name(part)
names = {'raw FRF magnitude', 'real part', 'imaginary part'};
s = names{part + 1};
end

function draw_part(ax, C, part, fs)
% part 0: raw FRF magnitude with its maximum; part 1: real part with points 3 and 4;
% part 2: imaginary part with point 1 and A.
line_col = [42 120 214]/255; pt_col = [235 104 52]/255; grey = [0.55 0.55 0.55];
f = C.freq;
switch part
    case 0, y = C.A_mag; ylab = sprintf('|FRF| as measured [%s]', C.u.raw);
        ttl = sprintf('%s: raw FRF magnitude (no smoothing)', C.name);
    case 1, y = C.Re; ylab = sprintf('Real part [%s]', C.u.frf);
        ttl = sprintf('%s: real part (smoothed: %d-point Savitzky-Golay)', C.name, C.smooth);
    case 2, y = C.Im; ylab = sprintf('Imaginary part [%s]', C.u.frf);
        ttl = sprintf('%s: imaginary part (smoothed: %d-point Savitzky-Golay)', C.name, C.smooth);
end
plot(ax, f, y, '-', 'Color', line_col, 'LineWidth', 1.5);
hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
if part > 0, plot(ax, [f(1) f(end)], [0 0], ':', 'Color', grey); end
set(ax, 'FontSize', fs);
xlim(ax, [f(1) f(end)]);
xlabel(ax, 'Frequency [Hz]');
ylabel(ax, ylab);
title(ax, ttl, 'Interpreter', 'none');
dx = 0.02*(f(end) - f(1));
switch part
    case 0, pts = 1;          txt = {sprintf('natural frequency (raw maximum): %.2f Hz,  |FRF| = %.4g', f(C.idx(1)), y(C.idx(1)))};
    case 1, pts = [3 4];      txt = {sprintf('3: %.2f Hz,  Re = %.4g', f(C.idx(3)), y(C.idx(3))), ...
                                     sprintf('4: %.2f Hz,  Re = %.4g', f(C.idx(4)), y(C.idx(4)))};
    case 2, pts = 2;          txt = {sprintf('1: f_{n1} = %.2f Hz,  A = %.4g', f(C.idx(2)), y(C.idx(2)))};
end
for j = 1:numel(pts)
    i = C.idx(pts(j));
    if part == 2, plot(ax, [f(1) f(i)], [y(i) y(i)], '--', 'Color', grey); end
    plot(ax, f(i), y(i), 'o', 'Color', pt_col, 'MarkerFaceColor', pt_col, 'MarkerSize', 7);
    text(f(i) + dx, y(i), txt{j}, 'Parent', ax, 'FontSize', fs - 1, 'VerticalAlignment', 'middle');
end
end

% =========================================================================
function write_report(R, C)
% One-page PDF (8.5 x 14 in): test information, natural frequency and picked
% points, modal parameters, notes and the three plots.
W = 8.5; Hh = 14;
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [0 0 W Hh]);
set(fig, 'PaperUnits', 'inches', 'PaperSize', [W Hh], 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 W Hh]);
ax = axes('Position', [0 0 1 1], 'Visible', 'off'); xlim(ax, [0 1]); ylim(ax, [0 1]); hold(ax, 'on');
x0 = 0.06; fs = 9; dy = 0.0125;
u = R.units;
T = @(x, y, s, varargin) text(x, y, s, 'Parent', ax, 'FontSize', fs, 'Interpreter', 'none', 'VerticalAlignment', 'middle', varargin{:});
rule = @(y) plot(ax, [x0 1 - x0], [y y], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
H1 = @(y, s) text(x0, y, s, 'Parent', ax, 'FontSize', 11, 'FontWeight', 'bold');

text(x0, 0.978, [R.name ' processed'], 'Parent', ax, 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');
T(x0, 0.962, 'SDOF peak-picking modal fit after Schmitz & Smith, Machining Dynamics, 2nd ed., Sec. 2.5.1 (Fig. 2.21, Eqs. 2.58, 2.60, 2.62, 2.63)', ...
    'Color', [0.32 0.32 0.31]);

% 1. test
y = 0.940;
H1(y, '1. Test'); y = y - 1.3*dy;
info = {
    'File', R.file
    'Experiment / FRF type', sprintf('%s (Y axis units in file: %s)', R.frf_type, R.y_units)
    'Conversion applied', R.conversion
    'Sign convention', sprintf('%+d (Fig. 2.21: Im < 0 at resonance)', R.polarity)
    'Data', sprintf('%d points, %.3f-%.3f Hz', R.n_points, R.f_range_Hz)
    'Search band', sprintf('%g-%g Hz', R.f_band_Hz)
    'Smoothing', sprintf('%d-point quadratic Savitzky-Golay (real and imaginary parts; raw maximum unsmoothed)', R.smooth_points)
    'Processed', datestr(now, 'yyyy-mm-dd HH:MM')};
for j = 1:size(info, 1)
    T(x0, y, info{j, 1}, 'FontWeight', 'bold'); T(x0 + 0.20, y, info{j, 2}); y = y - dy;
end

% 2. natural frequency and picked points
y = y - 0.6*dy;
H1(y, '2. Natural frequency and picked points (Fig. 2.21, one mode)'); y = y - 1.3*dy;
cx = x0 + [0 0.10 0.36 0.50 0.60 0.71 0.85];
hdr = {'Point', 'Meaning', 'Variable', 'f [Hz]', 'omega [rad/s]', 'Value', 'Coherence'};
for c = 1:numel(hdr), T(cx(c), y, hdr{c}, 'FontWeight', 'bold'); end
rule(y - dy/2); y = y - dy;
rows = {
    'raw max', 'maximum of the raw FRF magnitude', 'omega_n_raw', R.omega_n_raw, sprintf('|FRF| = %.4g [%s]', R.A_mag_max, u.raw)
    '1',       'minimum of the imaginary part',    'omega_n1, A', R.omega_n1,   sprintf('A = %.4g [%s]', R.A, u.frf)
    '3',       'biggest maximum of the real part', 'omega_3',     R.omega_3,    sprintf('Re = %.4g [%s]', R.Re_3, u.frf)
    '4',       'biggest minimum of the real part', 'omega_4',     R.omega_4,    sprintf('Re = %.4g [%s]', R.Re_4, u.frf)};
for j = 1:size(rows, 1)
    T(cx(1), y, rows{j, 1}); T(cx(2), y, rows{j, 2}); T(cx(3), y, rows{j, 3});
    T(cx(4), y, sprintf('%.3f', rows{j, 4}/(2*pi))); T(cx(5), y, sprintf('%.2f', rows{j, 4}));
    T(cx(6), y, rows{j, 5}); T(cx(7), y, sprintf('%.3f', R.point_coherence(j)));
    y = y - dy;
end

% 3. modal parameters
y = y - 0.6*dy;
H1(y, '3. Modal parameters'); y = y - 1.3*dy;
cx = x0 + [0 0.08 0.22 0.66 0.80];
hdr = {'Equation', 'Quantity', 'Variable', 'Value', 'Units'};
for c = 1:numel(hdr), T(cx(c), y, hdr{c}, 'FontWeight', 'bold'); end
rule(y - dy/2); y = y - dy;
rows = {
    '2.58', 'damping ratio', 'zeta_q1 = (omega_4 - omega_3)/(2 omega_n1)', sprintf('%.4f', R.zeta_q1), '-'
    '2.60', 'stiffness',     'k_q1 = -1/(2 zeta_q1 A)',                     sprintf('%.4g', R.k_q1),    u.k
    '2.62', 'mass',          'm_q1 = k_q1/omega_n1^2',                      sprintf('%.4g', R.m_q1),    u.m
    '2.63', 'damping',       'c_q1 = 2 zeta_q1 sqrt(k_q1 m_q1)',            sprintf('%.4g', R.c_q1),    u.c};
for j = 1:size(rows, 1)
    for c = 1:5, T(cx(c), y, rows{j, c}); end
    y = y - dy;
end

% 4. notes
y = y - 0.6*dy;
H1(y, '4. Notes'); y = y - 1.3*dy;
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
H1(y, '5. Raw FRF magnitude, real and imaginary parts');
top = y - 0.03; gap = 0.068; h = (top - 0.045 - 2*gap)/3;
for part = 0:2
    a = axes('Position', [0.11 top - (part + 1)*h - part*gap 0.84 h]);
    draw_part(a, C, part, 9);
end

print(fig, R.pdf_file, '-dpdf');
close(fig);
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
