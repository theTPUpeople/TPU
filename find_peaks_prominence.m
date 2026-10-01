function [idx, prom] = find_peaks_prominence(y)
%FIND_PEAKS_PROMINENCE Local maxima of a vector and their prominence.
%   [IDX, PROM] = FIND_PEAKS_PROMINENCE(Y) returns indices of local maxima of
%   Y and their topographic prominence (height above the higher of the two
%   lowest points between the peak and the nearest higher sample on each
%   side, or the signal end). Same definition as findpeaks' 'Prominence',
%   but without the Signal Processing Toolbox. To find minima, pass -Y.

y = y(:);
n = numel(y);
if n < 3, idx = zeros(0, 1); prom = zeros(0, 1); return, end
is_pk = [false; y(2:end-1) > y(1:end-2) & y(2:end-1) >= y(3:end); false];
idx = find(is_pk);
prom = zeros(size(idx));
for j = 1:numel(idx)
    i = idx(j);
    L = find(y(1:i-1) > y(i), 1, 'last');
    if isempty(L), L = 1; end
    R = find(y(i+1:end) > y(i), 1, 'first');
    if isempty(R), R = n; else, R = i + R; end
    prom(j) = y(i) - max(min(y(L:i)), min(y(i:R)));
end
end
