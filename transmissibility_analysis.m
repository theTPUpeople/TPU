function Tr = transmissibility_analysis(w, T, opts)
%TRANSMISSIBILITY_ANALYSIS Natural frequency, damping and isolation onset from |T|.
%   TR = TRANSMISSIBILITY_ANALYSIS(W, T, OPTS) for a base-excitation
%   transmissibility T(w) = output accel / base accel (dimensionless):
%     - peaks of |T| (up to OPTS.n_modes, by prominence >= OPTS.min_prominence_rel*max|T|)
%       give w_n (w at peak |T|); a peak lying inside the half-power band of a
%       more prominent peak is merged into it (TR.merged_f_Hz)
%     - zeta from the half-power bandwidth of |T| around each peak
%     - isolation onset: first frequency above the highest peak where |T|
%       drops below 1 (linear interpolation). For a viscous SDOF this is
%       exactly sqrt(2)*w_n for any zeta.
%   No force is measured, so mass cannot be identified from T.
%   Note: for a viscous SDOF the |T| peak lies slightly below w_n
%   (r_peak^2 = (sqrt(1+8 zeta^2) - 1)/(4 zeta^2)); e.g. 0.5 % low at zeta = 0.07.

w = w(:); mag = abs(T(:));
n = numel(w);
[ip, prom] = find_peaks_prominence(mag);
thr = opts.min_prominence_rel*max(mag);
pass = prom >= thr;
[~, order] = sort(prom, 'descend');
order = order(pass(order));
% a peak inside the half-power band of a more prominent, already accepted
% peak cannot be resolved by peak picking (noise split of one resonance):
% it is merged into that peak
sel = zeros(0, 1);
Tr.merged_f_Hz = zeros(0, 1);
bands = zeros(0, 2);
for c = order(:)'
    if numel(sel) >= opts.n_modes, break, end
    i = ip(c);
    if any(w(i) >= bands(:, 1) & w(i) <= bands(:, 2))
        Tr.merged_f_Hz(end+1, 1) = w(i)/(2*pi);
        continue
    end
    win = find(w >= w(i)*(1 - opts.real_search_frac) & w <= w(i)*(1 + opts.real_search_frac));
    hp0 = half_power_damping(w, mag, win(1), win(end));
    b1 = hp0.w1; b2 = hp0.w2;
    if ~isfinite(b1), b1 = w(win(1)); end
    if ~isfinite(b2), b2 = w(win(end)); end
    bands(end+1, :) = [b1 b2]; %#ok<AGROW>
    sel(end+1, 1) = i; %#ok<AGROW>
end
sel = sort(sel);

Tr.cand_f_Hz = w(ip)/(2*pi);
Tr.cand_mag = mag(ip);
Tr.cand_prom = prom;
Tr.threshold = thr;
Tr.peaks = struct([]);
for q = 1:numel(sel)
    i = sel(q);
    lo = w(i)*(1 - opts.real_search_frac);
    hi = w(i)*(1 + opts.real_search_frac);
    if q > 1,           lo = max(lo, 0.5*(w(sel(q-1)) + w(i))); end
    if q < numel(sel),  hi = min(hi, 0.5*(w(i) + w(sel(q+1)))); end
    win = find(w >= lo & w <= hi);
    hp = half_power_damping(w, mag, win(1), win(end));
    pk.q = q;
    pk.i_peak = hp.i_peak;
    pk.w_n = hp.w_peak;
    pk.f_n_Hz = hp.w_peak/(2*pi);
    pk.T_peak = hp.mag_peak;
    pk.zeta_hp = hp.zeta;
    pk.hp_w1 = hp.w1;
    pk.hp_w2 = hp.w2;
    pk.flags = hp.flags;
    if pk.T_peak <= 1
        pk.flags{end+1} = '|T| peak <= 1: no resonant amplification';
    end
    if isempty(Tr.peaks), Tr.peaks = pk; else, Tr.peaks(q) = pk; end
end

% isolation onset above the dominant (largest) peak
Tr.f_iso_Hz = NaN;
Tr.i_dom = [];
if ~isempty(Tr.peaks)
    [~, d] = max([Tr.peaks.T_peak]);
    Tr.i_dom = d;
    i0 = Tr.peaks(d).i_peak;
    j = find(mag(i0:n) < 1, 1, 'first');
    if ~isempty(j)
        j = i0 + j - 1;
        Tr.f_iso_Hz = (w(j-1) + (mag(j-1) - 1)*(w(j) - w(j-1))/(mag(j-1) - mag(j)))/(2*pi);
    end
end
end
