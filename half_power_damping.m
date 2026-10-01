function hp = half_power_damping(w, mag, ilo, ihi)
%HALF_POWER_DAMPING Half-power (-3 dB) bandwidth damping estimate.
%   HP = HALF_POWER_DAMPING(W, MAG, ILO, IHI) takes the largest MAG within
%   indices ILO:IHI as the peak, finds where MAG falls below peak/sqrt(2) on
%   each side (linear interpolation between bins, search limited to
%   ILO:IHI) and returns
%       zeta = (w2 - w1) / (2 w_peak)
%   HP fields: zeta, w1, w2, w_peak, i_peak, mag_peak, ok, flags.
%   ok = false (zeta = NaN) if either crossing is not found inside ILO:IHI.

w = w(:); mag = mag(:);
hp.flags = {};
[mp, k] = max(mag(ilo:ihi));
i = ilo + k - 1;
hp.i_peak = i; hp.w_peak = w(i); hp.mag_peak = mp;
lvl = mp/sqrt(2);

j = i;
while j > ilo && mag(j) > lvl, j = j - 1; end
if mag(j) > lvl
    hp.w1 = NaN; hp.flags{end+1} = 'lower half-power point not found in window';
else
    hp.w1 = w(j) + (lvl - mag(j))*(w(j+1) - w(j))/(mag(j+1) - mag(j));
end
j = i;
while j < ihi && mag(j) > lvl, j = j + 1; end
if mag(j) > lvl
    hp.w2 = NaN; hp.flags{end+1} = 'upper half-power point not found in window';
else
    hp.w2 = w(j-1) + (mag(j-1) - lvl)*(w(j) - w(j-1))/(mag(j-1) - mag(j));
end
hp.zeta = (hp.w2 - hp.w1)/(2*hp.w_peak);
hp.ok = isfinite(hp.zeta) && hp.zeta > 0;
if ~hp.ok, hp.zeta = NaN; end
end
