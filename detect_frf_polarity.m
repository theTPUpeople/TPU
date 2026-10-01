function [s, info] = detect_frf_polarity(w, G, w_ref)
%DETECT_FRF_POLARITY Sign that puts an FRF in the textbook convention.
%   [S, INFO] = DETECT_FRF_POLARITY(W, G) returns S = +1 or -1 such that S*G
%   has Im(S*G) < 0 at the largest |G| (the dominant resonance), as in
%   textbook Fig. 2.21 (receptance with x = X e^{i w t}: Re > 0 below
%   resonance, Im < 0 at resonance).
%
%   [S, INFO] = DETECT_FRF_POLARITY(W, G, W_REF) evaluates the sign at the bin
%   nearest W_REF [rad/s] instead (e.g. at the |T| peak, when |G| is dominated
%   by low-frequency error elsewhere).
%
%   An overall sign inversion (hammer and accelerometer positive directions
%   opposed) flips Re and Im together. A conjugate (e^{-i w t}) convention
%   flips Im only. INFO.re_consistent reports whether Re(S*G) is > 0 just
%   below and < 0 just above the dominant peak; if Im is consistent but Re is
%   not, the data do not follow the textbook form and the caller should warn.
%   An overall sign does not change w_n, zeta or |Im| at resonance.

G = G(:); w = w(:);
if nargin < 3 || isempty(w_ref)
    [~, i] = max(abs(G));
else
    [~, i] = min(abs(w - w_ref));
end
s = -sign(sum(imag(G(max(1, i-3):min(numel(G), i+3)))));   % +/- 3 bins: robust to noise
if s == 0, s = 1; end
Gs = s*G;
wp = w(i);
near = abs(G) > 0.5*abs(G(i));
info.w_peak = wp;
info.im_frac_negative = mean(imag(Gs(near)) < 0);
below = w >= 0.85*wp & w <= 0.97*wp;
above = w >= 1.03*wp & w <= 1.15*wp;
if any(below) && any(above)
    info.re_consistent = mean(real(Gs(below))) > 0 && mean(real(Gs(above))) < 0;
else
    info.re_consistent = true;   % not enough points to test
end
end
