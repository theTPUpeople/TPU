clc
clear all 
close all

% Load FRF data ( Freq, Real, Imag)
cd('H:\Boeing Project\Data\Dynamometer Transfer Functions\Makino Mounted without Aluminum Plate')
load('dyna_tf_makino_no_AL_plate_Z.mat')

freq = Freq;
omega = freq*2*pi;
real_X_F_measured = Real;
imag_X_F_measured = Imag;
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
 
real_filtered = filter(b,a,real_X_F_measured);
imag_filtered = filter(b,a,imag_X_F_measured);

%**************************************************************************
% Remove all FRF data below 'minumum FRF frequency' due to integration
% noise
minFreq =0;                          % minimum FRF frequency
maxFreq = 9000;                         % maximum FRF frequency
Freq_index = find(freq >= minFreq & freq <= maxFreq);
freq = freq(Freq_index);
omega = omega(Freq_index);
real_X_F_measured = real_X_F_measured(Freq_index);
imag_X_F_measured = imag_X_F_measured(Freq_index);
real_filtered = real_filtered(Freq_index);
imag_filtered = imag_filtered(Freq_index);

mag_X_F_measured = sqrt(real_filtered.^2 + imag_filtered.^2);
phase_X_F_measured = atan2(imag_filtered,real_filtered)/2/pi;

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
prompt = 'Enter the number of peaks from Figure 3 to be used.';
num_DOF = input(prompt);

for cnt = 1:num_DOF
   prompt = 'Enter the number of the peak to be kept.';
   index(cnt) = input(prompt);
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
    
    if pos == 1
        if abs((natural_freqs(pos + 1) - natural_freqs(pos))*0.40) < abs(natural_freqs(pos) - min(freq)) 
            f_lower = natural_freqs(pos) - ((natural_freqs(pos + 1) - natural_freqs(pos))*0.40);
            f_upper = natural_freqs(pos) + ((natural_freqs(pos + 1) - natural_freqs(pos))*0.40);
            
        elseif abs(natural_freqs(pos) - min(freq)) < abs((natural_freqs(pos + 1) - natural_freqs(pos))*0.40)
            f_lower = min(freq);
            f_upper = natural_freqs(pos) + (natural_freqs(pos) - min(freq));
            
        end
        
    elseif pos == num_modes
        if abs((natural_freqs(pos) - natural_freqs(pos-1))*0.40) < abs(max(freq) - natural_freqs(pos))
            f_lower = natural_freqs(pos) - ((natural_freqs(pos) - natural_freqs(pos - 1)) * 0.40);
            f_upper = natural_freqs(pos) + ((natural_freqs(pos) - natural_freqs(pos - 1)) * 0.40);
            
        elseif abs(max(freq) - natural_freqs(pos)) < abs((natural_freqs(pos) - natural_freqs(pos-1))*0.40)
            f_lower = natural_freqs(pos) - (max(freq) - natural_freqs(pos));
            f_upper = max(freq);
            
        end
        
    else
        if abs((natural_freqs(pos) - natural_freqs(pos - 1)) * 0.40) < abs((natural_freqs(pos + 1) - natural_freqs(pos))*0.40)
            f_lower = natural_freqs(pos) - ((natural_freqs(pos) - natural_freqs(pos - 1)) * 0.40);
            f_upper = natural_freqs(pos) + ((natural_freqs(pos) - natural_freqs(pos - 1))*0.40);
            
        elseif abs((natural_freqs(pos + 1) - natural_freqs(pos))*0.40) < abs((natural_freqs(pos) - natural_freqs(pos - 1))*0.40)
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
    
    wn(cnt1) = fn(cnt1)*(2*pi);                 % convert to [rad/s]
    w1(cnt1) = f1(cnt1)*(2*pi);                 % convert to [rad/s]
    w2(cnt1) = f2(cnt1)*(2*pi);                 % convert to [rad/s]
    
    zeta_q(cnt1) = (w2(cnt1)-w1(cnt1))/(2*wn(cnt1));            % modal damping ratio [unitless]
    k_q(cnt1) = -1/(2*zeta_q(cnt1)*A(cnt1));                    % modal stiffness [N/m]
    m_q(cnt1) = k_q(cnt1)/wn(cnt1)^2;                           % modal mass [kg]
    c_q(cnt1) = 2*zeta_q(cnt1)*sqrt(k_q(cnt1)*m_q(cnt1));       % modal damping coefficient [(N-s)/m]
    
    r = (omega/wn(cnt1));
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
x0 = [wn; k_q; zeta_q];
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
cd('H:\Boeing Project\Matlab Toolbox\FRF Modal Parameter Fitting')
[x] = modalfit(freq,real_filtered,imag_filtered,x0)

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
    fn(cnt) = x(1 + (cnt - 1)*3)/(2*pi);
    k_q(cnt) = x(2 + (cnt - 1)*3);
    zeta_q(cnt) = x(3 + (cnt - 1)*3);
end

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

cd('H:\Boeing Project\Journal Papers\ASPE Journal\Cutting Force Coefficient Paper\Inverse Filter Code and Figures');

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
ylabel('Phase (char(176))')

