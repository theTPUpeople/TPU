clc
clear all
close all

%% ========================= CONFIGURATION (edit only here) =========================
% Peak-picking modal fit after Schmitz & Smith, Machining Dynamics 2nd ed.,
% Sec. 2.5.1 (Fig. 2.21, Eqs. 2.58-2.63). 'batch' runs the QuickDAQ CSV
% pipeline over all TPU specimens (tpu_modal_batch.m); 'legacy' runs the
% original single-FRF workflow below on a .mat file with Freq, Real, Imag.
cfg.m_known_kg   = NaN;   % TODO: preload mass in kg (plates + accelerometer). Static preload target 1.1 kN gives about 112 kg; verify actual.
cfg.frf_type     = 'auto'; % 'auto' infers from file/headers; override if wrong ('accelerance' | 'receptance', force-input test only)
cfg.n_modes      = 2;
cfg.min_prominence = [];  % tune on real data. Absolute prominence of the Im(H_rec) minima [m/N]; [] -> min_prominence_rel*max|Im|
cfg.min_prominence_rel = 0.05; % used when min_prominence = [] (also for |T| peaks)
cfg.polarity     = 'auto'; % 'auto' | +1 | -1 : sign giving Im(H_rec) < 0 at resonance (textbook Fig. 2.21 convention)
cfg.f_band_Hz    = [30 500]; % analysis band [Hz]. From 30 Hz up >= 95 % of D-test bins have coherence >= 0.75 (58 % at 15-25 Hz);
                            % below 30 Hz 98-100 % of T-test bins fail it. /w^2 amplifies low-frequency noise.
cfg.min_coherence = 0.75; % points with lower coherence are skipped (threshold used in the T spreadsheets' "Coherence valid?" column)
cfg.smooth_points = 1;    % D test: centred moving average length on the complex FRF (1 = off, like the original windowSize = 1);
                          % measured |H| scatter 0.9 % median / 1.75 % max, inside the validated 1-2 % range
cfg.smooth_points_T = 11;      % T test only: |T| scatter is 3-16 % per bin (median 6 %). Unsmoothed, half-power zeta is biased
cfg.smooth_method_T = 'sgolay'; % -14..-19 % and missing in 55-86 % of Monte-Carlo runs; quadratic Savitzky-Golay over 11 bins gave
                               % zeta RMS 6-9 % at 6 % scatter with 2.5-6 % noise-free bias (box-7: 9-17 % bias)
cfg.peak_refine  = 'sdof'; % sub-bin refinement of points 1-6: 'sdof' (parabola in SDOF-symmetric coordinates) | 'none' (raw bins) | h (parabola over +/- h bins)
cfg.real_search_frac = 0.5; % Re max/min searched within w_n*(1 -/+ frac), clipped at neighbouring modes
cfg.disagree_pct = 15;    % cross-check tolerance [%] (methods and D vs T)
cfg.phase_units  = 'deg';  % phase angles in the exports are in degrees ('auto' | 'deg' | 'rad')
cfg.data_root    = fileparts(mfilename('fullpath'));
if isempty(cfg.data_root), cfg.data_root = pwd; end
cfg.results_dir  = fullfile(cfg.data_root, 'results');
% one specimen of each test, read at the top of this script (CSV exports, open in Excel)
cfg.disturbance_file      = fullfile(cfg.data_root, 'Disturbance Rejection Test', '0%', 'TPU0S1D.csv');
cfg.transmissibility_file = fullfile(cfg.data_root, 'Transmissibility test', '0%', 'TPU0S1T.csv');
% SDOF peak picking and modal parameters (process_modal_file.m): the same steps and settings for
% every file and both experiment types (disturbance rejection and transmissibility)
cfg.process_files      = {cfg.transmissibility_file};   % list of exports to process, any mix of D and T files
cfg.pick_smooth_points = 11;    % quadratic Savitzky-Golay window [bins, odd] for the curves and the picking; 1 = raw data
cfg.pick_show_plots    = true;  % open the raw-magnitude, real and imaginary figures of each file
cfg.report_dir         = fullfile(cfg.results_dir, 'processed');   % "<file name> processed.pdf" is written here
cfg.make_figures = true;
cfg.run_mode     = 'batch'; % 'batch' | 'legacy'
% legacy mode only (original hard-coded values kept as defaults)
cfg.legacy_source     = 'csv'; % 'csv': FRF extracted at the top of this script | 'mat': original .mat with Freq, Real, Imag
cfg.legacy_test       = 'D';   % with 'csv': 'D' disturbance rejection | 'T' transmissibility
cfg.legacy_mat_file   = 'H:\Boeing Project\Data\Dynamometer Transfer Functions\Makino Mounted without Aluminum Plate\dyna_tf_makino_no_AL_plate_Z.mat';
cfg.legacy_peak_index = []; % [] -> prompt for peaks as before; e.g. [1 2] runs non-interactively
cfg.legacy_fig_dir    = ''; % original: 'H:\Boeing Project\Journal Papers\ASPE Journal\Cutting Force Coefficient Paper\Inverse Filter Code and Figures'
%% ====================================================================================
addpath(cfg.data_root);

%% ------------- FRF magnitude and phase from the CSV (Excel) exports -------------
% QuickDAQ exports: 9 header lines, numbers from line 10 on. Column layout:
%   Disturbance rejection: 1 Frequency [Hz] | 2 FRF magnitude [m/s^2/N] | 3 FRF phase [deg] | 4 coherence
%   Transmissibility:      1 Frequency [Hz] | 2 coherence | 3 FRF magnitude [-]     | 4 FRF phase [deg]
% The "Measurement Type" header line (line 5) is checked so a different layout stops here.
hdr_D = regexp(fileread(cfg.disturbance_file), '\r?\n', 'split');
hdr_T = regexp(fileread(cfg.transmissibility_file), '\r?\n', 'split');
lay_D = 'Measurement Type,FRF,,Coherence';
lay_T = 'Measurement Type,Coherence,FRF';
if ~strncmp(hdr_D{5}, lay_D, numel(lay_D)) || ~strncmp(hdr_T{5}, lay_T, numel(lay_T))
    error('FRF:layout', 'Column layout differs from the QuickDAQ export expected here; check line 5 of the files.');
end

% Disturbance rejection (force input, acceleration output)
data_D      = dlmread(cfg.disturbance_file, ',', 9, 0);   % numeric block, 9 header lines skipped
freq_D      = data_D(:, 1);                 % frequency [Hz]
A_D         = data_D(:, 2);                 % FRF magnitude A [m/s^2/N]
phase_D_deg = data_D(:, 3);                 % FRF phase [deg]
coh_D       = data_D(:, 4);                 % coherence
phase_D_rad = phase_D_deg*pi/180;           % phase [rad]
Re_D        = A_D.*cos(phase_D_rad);        % Re = A cos(phase)
Im_D        = A_D.*sin(phase_D_rad);        % Im = A sin(phase)

% Transmissibility (base acceleration input, mass acceleration output)
data_T      = dlmread(cfg.transmissibility_file, ',', 9, 0);
freq_T      = data_T(:, 1);                 % frequency [Hz]
coh_T       = data_T(:, 2);                 % coherence
A_T         = data_T(:, 3);                 % FRF magnitude A [-]
phase_T_deg = data_T(:, 4);                 % FRF phase [deg]
phase_T_rad = phase_T_deg*pi/180;           % phase [rad]
Re_T        = A_T.*cos(phase_T_rad);        % Re = A cos(phase)
Im_T        = A_T.*sin(phase_T_rad);        % Im = A sin(phase)

fprintf('Read %d points from %s and %d points from %s (phase deg -> rad, Re = A cos, Im = A sin)\n', ...
    numel(freq_D), cfg.disturbance_file, numel(freq_T), cfg.transmissibility_file);

%% ------------- Peak picking and modal parameters for every file in cfg.process_files -------------
% process_modal_file.m applies the same SDOF steps to each file: read -> Re = A cos(phase), Im = A sin(phase)
% -> (accelerance only) receptance H/(-omega^2) -> Fig. 2.21 sign -> smoothing
% -> natural frequency from the raw data (maximum of the measured FRF magnitude)
% -> point 1 (Im minimum: omega_n1, A), point 3 (biggest Re maximum), point 4 (biggest Re minimum)
% -> Eqs. 2.58, 2.60, 2.62, 2.63 -> plots -> "<file name> processed.pdf" in cfg.report_dir.
% The results (named after the textbook symbols) are collected in the struct array "processed".
clear processed
for i_file = 1:numel(cfg.process_files)
    processed(i_file) = process_modal_file(cfg.process_files{i_file}, cfg); %#ok<SAGROW>
end

if strcmpi(cfg.run_mode, 'batch')
    [results, xcheck] = tpu_modal_batch(cfg);
    return
end

%% ---------------- legacy mode: original single-FRF workflow ----------------
% Load FRF data ( Freq, Real, Imag)
if strcmpi(cfg.legacy_source, 'mat')
    load(cfg.legacy_mat_file)
else
    % FRF extracted from the CSV exports at the top of this script
    if strcmpi(cfg.legacy_test, 'T')
        Freq = freq_T; Real = Re_T; Imag = Im_T;
        warning('FRF:legacyT', ['Transmissibility is dimensionless (no force): k_q, m_q, c_q from ', ...
            'Eqs. 2.60-2.63 are not physical for T; only fn and zeta are meaningful.']);
    else
        keep = freq_D > 0;                                    % w = 0 cannot be divided by
        Freq = freq_D(keep); Real = Re_D(keep); Imag = Im_D(keep);
        is_acc = strcmpi(cfg.frf_type, 'accelerance') || (strcmpi(cfg.frf_type, 'auto') && ...
            ~isempty(strfind(strrep(hdr_D{8}, ' ', ''), 'm/s^2/N')));
        if is_acc
            % accelerance -> receptance: Eqs. 2.58-2.63 assume displacement/force
            H_rec = (Real + 1i*Imag)./(-(2*pi*Freq).^2);
            Real = real(H_rec); Imag = imag(H_rec);
            fprintf('legacy: H_rec = H_acc/(-w^2) applied\n');
        end
    end
    % sign convention of textbook Fig. 2.21 (Im < 0 at resonance), checked inside the analysis band
    if ischar(cfg.polarity)
        inb = Freq >= cfg.f_band_Hz(1) & Freq <= cfg.f_band_Hz(2);
        sgn = detect_frf_polarity(2*pi*Freq(inb), Real(inb) + 1i*Imag(inb));
    else
        sgn = sign(cfg.polarity);
    end
    Real = sgn*Real; Imag = sgn*Imag;
    fprintf('legacy: polarity %+d applied\n', sgn);
end

% row vectors throughout (column data broadcast to N x N further down)
freq = Freq(:).';
omega = freq*2*pi;
real_X_F_measured = Real(:).';
imag_X_F_measured = Imag(:).';
clearvars Real1 Imag1 Freq1

figure(1)

subplot(211)
plot(freq,real_X_F_measured, 'LineWidth', 1.5)
set(gca,'FontSize',14)
title('Dynamometer X-Direction FRF')
xlabel('Frequency [Hz]')
ylabel('Real[N/N]')
xlim([0 max(freq)])
subplot(212)
plot(freq,imag_X_F_measured, 'LineWidth', 1.5)
set(gca,'FontSize',14)
xlabel('Frequency [Hz]')
ylabel('Imaginary[N/N]')
xlim([0 max(freq)])
%%
%**************************************************************************
% Apply a moving average filter represented by the following difference
% equation: 

% y(n) = 1/windowSize * (X(n) + X(n-1) + ... + X(n-(windowSize-1))
windowSize = 1;
b = (1/windowSize)*ones(1,windowSize);
a = 1;                                      % denominator coefficient
 
% centred moving average: filter() is causal and shifts peaks by
% (windowSize-1)/2 bins when windowSize > 1 (identical for windowSize = 1)
real_filtered = conv(real_X_F_measured, b, 'same');
imag_filtered = conv(imag_X_F_measured, b, 'same');

%**************************************************************************
% Remove all FRF data below 'minumum FRF frequency' due to integration
% noise
minFreq =0;                          % minimum FRF frequency
maxFreq = 9000;                         % maximum FRF frequency
if ~strcmpi(cfg.legacy_source, 'mat'), minFreq = cfg.f_band_Hz(1); maxFreq = cfg.f_band_Hz(2); end
Freq_index = find(freq >= minFreq & freq <= maxFreq);
freq = freq(Freq_index);
omega = omega(Freq_index);
real_X_F_measured = real_X_F_measured(Freq_index);
imag_X_F_measured = imag_X_F_measured(Freq_index);
real_filtered = real_filtered(Freq_index);
imag_filtered = imag_filtered(Freq_index);

mag_X_F_measured = sqrt(real_filtered.^2 + imag_filtered.^2);
phase_X_F_measured = atan2(imag_filtered,real_filtered)*(180/pi);   % [deg] (was /2/pi = cycles, plotted against degrees in Fig. 8)

figure(2)
subplot(211)
plot(freq,real_X_F_measured,'k',freq,real_filtered, 'LineWidth', 1.5)
set(gca,'FontSize',14)
title('Dynamometer X-Direction Force-to-Force FRF')
% legend('Raw Data','Filtered Data')
xlabel('Frequency [Hz]')
ylabel('Real[N/N]')
xlim([minFreq maxFreq])
ylim([-15 15])
subplot(212)
plot(freq,imag_X_F_measured,'k',freq,imag_filtered, 'LineWidth', 1.5)
set(gca,'FontSize',14)
% legend('Raw Data','Filtered Data')
xlabel('Frequency [Hz]')
ylabel('Imaginary[N/N]')
xlim([minFreq maxFreq])
ylim([-30 1])
%%
%**************************************************************************
% Find local minima of the Imaginary part of the FRF to identify natural
% frequencies
if max(imag_filtered) > abs(min(imag_filtered))
    warning('FRF:polarity', ['Im(FRF) is mostly positive: the sign is inverted relative to the textbook ', ...
        'convention (Im < 0 at resonance). Multiply Real and Imag by -1 before peak picking.']);
end
[imag_pks,imag_locs] = findpeaks(-imag_filtered,'MinPeakDistance',5,'MinPeakHeight',0.01*max(-imag_filtered));
imag_pks = -imag_pks;                                                                       
natural_freqs = freq(imag_locs);                                                          % create a vector of natural frequencies [Hz]

figure(3)
plot(freq, imag_filtered, natural_freqs, imag_pks, 'o', 'LineWidth', 1.5)
set(gca,'FontSize',14)
text(natural_freqs+10,imag_pks-2*10^-9,num2str((1:numel(imag_pks))'))
xlabel('Frequency [Hz]')
ylabel('Imaginary[N/N]')
xlim([minFreq maxFreq])

%**************************************************************************
% Prompt the user to enter the number of DOFs to be modeled.
if isempty(cfg.legacy_peak_index)
    prompt = 'Enter the number of peaks from Figure 3 to be used.';
    num_DOF = input(prompt);

    for cnt = 1:num_DOF
       prompt = 'Enter the number of the peak to be kept.';
       index(cnt) = input(prompt);
    end
else
    index = cfg.legacy_peak_index;
    num_DOF = numel(index);
end

natural_freqs = natural_freqs(index);            % create a vector of natural frequencies [Hz]
imag_pks = imag_pks(index);
 
figure(4)
plot(freq, imag_filtered, natural_freqs, imag_pks, 'o', 'LineWidth', 1.5)
set(gca,'FontSize',14)
xlabel('Frequency [Hz]')
ylabel('Imaginary[m/N]')
xlim([minFreq maxFreq])

%%
%**************************************************************************

num_modes = length(imag_pks);                                                               % counts number of vibration modes [int]
real_X_F = real_filtered;
imag_X_F = imag_filtered;
Q_R_total = zeros(1,length(freq));

for cnt1 = 1:num_modes
    
    [val,pos] = min(imag_pks);                                                                % returns the peak value and index of the most negative peak in the imaginary part of the FRF
    fn(cnt1) = natural_freqs(pos);                                                              % assigns the frequency of the peak (natural frequency)
    A(cnt1) = val;                                                                              % assigns the amplitude of the peak
    
    % (each if/elseif pair below was missing the tie case, leaving f_lower /
    % f_upper undefined or stale; the elseif's are now else's. A single mode
    % indexed natural_freqs(pos + 1) out of range; it now uses the full band.)
    if num_modes == 1
        f_lower = min(freq);
        f_upper = max(freq);

    elseif pos == 1
        if abs((natural_freqs(pos + 1) - natural_freqs(pos))*0.40) < abs(natural_freqs(pos) - min(freq)) 
            f_lower = natural_freqs(pos) - ((natural_freqs(pos + 1) - natural_freqs(pos))*0.40);
            f_upper = natural_freqs(pos) + ((natural_freqs(pos + 1) - natural_freqs(pos))*0.40);
            
        else
            f_lower = min(freq);
            f_upper = natural_freqs(pos) + (natural_freqs(pos) - min(freq));
            
        end
        
    elseif pos == num_modes
        if abs((natural_freqs(pos) - natural_freqs(pos-1))*0.40) < abs(max(freq) - natural_freqs(pos))
            f_lower = natural_freqs(pos) - ((natural_freqs(pos) - natural_freqs(pos - 1)) * 0.40);
            f_upper = natural_freqs(pos) + ((natural_freqs(pos) - natural_freqs(pos - 1)) * 0.40);
            
        else
            f_lower = natural_freqs(pos) - (max(freq) - natural_freqs(pos));
            f_upper = max(freq);
            
        end
        
    else
        if abs((natural_freqs(pos) - natural_freqs(pos - 1)) * 0.40) < abs((natural_freqs(pos + 1) - natural_freqs(pos))*0.40)
            f_lower = natural_freqs(pos) - ((natural_freqs(pos) - natural_freqs(pos - 1)) * 0.40);
            f_upper = natural_freqs(pos) + ((natural_freqs(pos) - natural_freqs(pos - 1))*0.40);
            
        else
            f_lower = natural_freqs(pos) - ((natural_freqs(pos + 1) - natural_freqs(pos))*0.40);
            f_upper = natural_freqs(pos) + ((natural_freqs(pos + 1) - natural_freqs(pos))*0.40);
            
        end
            
    end
    
    f_index = find(freq >= f_lower & freq <= f_upper);                                        % defines an index pointing to a valid frequency range for the selected vibration mode
    f_range = freq(f_index);                                                                  % defines valid frequency range
    real_range = real_filtered(f_index);                                                      % assigns the real part of the FRF in the valid frequency range
    [f1_pks,f1_locs] = max(real_range);
    f1(cnt1) = f_range(f1_locs);
    [f2_pks,f2_locs] = max(-real_range);
    f2(cnt1) = f_range(f2_locs);
    
    if f1(cnt1) > f2(cnt1)
        f1(cnt1) = f_range(f2_locs);
        f2(cnt1) = f_range(f1_locs);
    end
    
    % Textbook notation (Sec. 2.5.1, Fig. 2.21): for mode 1 omega_n = omega_n1, omega_low = omega_3,
    % omega_high = omega_4 and A = A; for mode 2 omega_n = omega_n2, omega_low = omega_5, omega_high = omega_6, A = B.
    omega_n(cnt1) = fn(cnt1)*(2*pi);            % natural frequency (point 1 / 2) [rad/s]
    omega_low(cnt1) = f1(cnt1)*(2*pi);          % Re maximum below omega_n (point 3 / 5) [rad/s]
    omega_high(cnt1) = f2(cnt1)*(2*pi);         % Re minimum above omega_n (point 4 / 6) [rad/s]
    
    zeta_q(cnt1) = (omega_high(cnt1)-omega_low(cnt1))/(2*omega_n(cnt1));   % Eq. 2.58 / 2.59: modal damping ratio [unitless]
    k_q(cnt1) = -1/(2*zeta_q(cnt1)*A(cnt1));                    % Eq. 2.60 / 2.61: modal stiffness [N/m]
    m_q(cnt1) = k_q(cnt1)/omega_n(cnt1)^2;                      % Eq. 2.62: modal mass [kg]
    c_q(cnt1) = 2*zeta_q(cnt1)*sqrt(k_q(cnt1)*m_q(cnt1));       % Eq. 2.63: modal damping coefficient [(N-s)/m]
    
    r = (omega/omega_n(cnt1));
    Q_R(cnt1,:) = (1/k_q(cnt1))*(((1-r.^2)-1i*(2*zeta_q(cnt1)*r))./((1-r.^2).^2+(2*zeta_q(cnt1)*r).^2));
    real_Q_R(cnt1,:) = real(Q_R(cnt1,:));
    imag_Q_R(cnt1,:) = imag(Q_R(cnt1,:));
    
    real_X_F = real_X_F - real_Q_R(cnt1,:);
    imag_X_F = imag_X_F - imag_Q_R(cnt1,:);
    
    imag_pks(pos) = 0;
            
    Q_R_total = Q_R_total + Q_R(cnt1,:);
    
    figure(5)
    plot(f_range,real_range)
    set(gca,'FontSize',14)
    xlabel('Frequency [Hz]')
    ylabel('Real[m/N]')
        
    figure(6)
    subplot(211)
    plot(freq, real_Q_R(cnt1,:), ':', freq,real_filtered, 'LineWidth', 1.5)
    set(gca,'FontSize',14)
    xlabel('Frequency [Hz]')
    ylabel('Real[N/N]')
    hold on
    
    subplot(212)
    plot(freq, imag_Q_R(cnt1,:), ':',freq,imag_filtered, 'LineWidth', 1.5)
    set(gca,'FontSize',14)
    xlabel('Frequency [Hz]')
    ylabel('Imaginary[N/N]')
    hold on
end
% Order the modes by frequency so that index q follows the textbook numbering (mode 1 = lowest
% omega_n); the loop above handled them from the deepest imaginary peak down.
[~, q_order] = sort(fn);
fn = fn(q_order); f1 = f1(q_order); f2 = f2(q_order); A = A(q_order);
omega_n = omega_n(q_order); omega_low = omega_low(q_order); omega_high = omega_high(q_order);
zeta_q = zeta_q(q_order); k_q = k_q(q_order); m_q = m_q(q_order); c_q = c_q(q_order);

x0 = [omega_n; k_q; zeta_q];
x0 = x0(:)';


mag_Q_R_total = sqrt(real(Q_R_total).^2 + imag(Q_R_total).^2);
phase_Q_R_total = atan2(imag(Q_R_total),real(Q_R_total))*(180/pi);

figure(7)
subplot(211)
plot(freq, real_filtered, freq, real(Q_R_total), ':', 'LineWidth', 1.5)
set(gca,'FontSize',14)
title('Fit from Peak Picking')
legend('Measured FRF','Modal Fit')
xlabel('Frequency [Hz]')
ylabel('Real[N/N]')
subplot(212)
plot(freq,imag_filtered,freq,imag(Q_R_total), ':', 'LineWidth', 1.5)
set(gca,'FontSize',14)
legend('Measured FRF','Modal Fit')
xlabel('Frequency [Hz]')
ylabel('Imaginary[m/N]')

figure(8)
subplot(211)
plot(freq,mag_X_F_measured,freq,mag_Q_R_total)
set(gca,'FontSize',14)
title('Original Fit')
legend('Measured FRF','Modal Fit')
xlabel('Frequency [Hz]')
ylabel('Magnitude [m/N]')
subplot(212)
plot(freq,phase_X_F_measured,freq,phase_Q_R_total)
set(gca,'FontSize',14)
legend('Measured FRF','Modal Fit')
xlabel('Frequency [Hz]')
ylabel('Phase [deg]')
%%
[x] = modalfit(freq,real_filtered,imag_filtered,x0)          % modalfit.m sits next to this script

Q_R_total = zeros(1,length(freq));

for cnt2 = 1:num_modes
    r = (freq*2*pi)/x(1, (1 + (cnt2 - 1)*3));
    Q_R(cnt2,:) = (1/x(1, (2 + (cnt2 - 1)*3)))*(((1-r.^2)-1i*(2*x(1, (3 + (cnt2 - 1)*3))*r))./((1-r.^2).^2+(2*x(1, (3 + (cnt2 - 1)*3))*r).^2));
    Q_R_total = Q_R_total + Q_R(cnt2,:);
end

X_F_Mag = sqrt(real_filtered.^2 + imag_filtered.^2);
X_F_Phase = atan2(imag_filtered, real_filtered)*(180/pi);

Q_R_Mag = sqrt(real(Q_R_total).^2 + imag(Q_R_total).^2);
Q_R_Phase = atan2(imag(Q_R_total), real(Q_R_total))*(180/pi);
Q_R_Real = real(Q_R_total);
Q_R_Imag = imag(Q_R_total); 

for cnt = 1:num_modes
    omega_n(cnt) = x(1 + (cnt - 1)*3);
    fn(cnt) = omega_n(cnt)/(2*pi);
    k_q(cnt) = x(2 + (cnt - 1)*3);
    zeta_q(cnt) = x(3 + (cnt - 1)*3);
    m_q(cnt) = k_q(cnt)/omega_n(cnt)^2;                         % Eq. 2.62 with the optimised omega_n, k_q
    c_q(cnt) = 2*zeta_q(cnt)*sqrt(k_q(cnt)*m_q(cnt));           % Eq. 2.63 with the optimised zeta_q, k_q, m_q
end
K_q = diag(k_q);                                                % modal matrices [K_q], [M_q], [C_q] (p. 43)
M_q = diag(m_q);
C_q = diag(c_q);

figure(9)
subplot(211)
plot(freq, real_filtered, freq, Q_R_Real, ':', 'LineWidth', 1.5)
set(gca,'FontSize', 14)
title('Optimized Fit')
legend('Measured','Fit')
xlabel('Frequency (Hz)')
ylabel('Real (N/N)')
% xlim([1900 2100])
subplot(212)
plot(freq, imag_filtered, freq, Q_R_Imag, ':', 'LineWidth', 1.5)
set(gca,'FontSize', 14)
xlabel('Frequency (Hz)')
ylabel('Imaginary (N/N)')
% xlim([1900 2100])

if ~isempty(cfg.legacy_fig_dir), cd(cfg.legacy_fig_dir); end

figure(10)
subplot(211)
plot(freq, X_F_Mag, freq,mag_Q_R_total, freq, Q_R_Mag)
set(gca,'FontSize', 14)
legend('Measured FRF', 'Peak Picking Fit', 'Optimized Fit')
xlabel('Frequency (Hz)')
ylabel('Magnitude (m/N)')
% magnifyOnFigure(fig)
subplot(212)
plot(freq, X_F_Phase, freq,phase_Q_R_total, freq, Q_R_Phase)
set(gca,'FontSize', 14)
xlabel('Frequency (Hz)')
ylabel(['Phase (' char(176) ')'])

