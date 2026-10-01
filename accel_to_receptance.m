function H_rec = accel_to_receptance(H_acc, w)
%ACCEL_TO_RECEPTANCE Convert accelerance to receptance.
%   H_REC = ACCEL_TO_RECEPTANCE(H_ACC, W) returns
%       H_rec(w) = H_acc(w) / (-w^2)
%   with W in rad/s. For harmonic motion x = X e^{i w t}, a = -w^2 x, so an
%   accelerance [(m/s^2)/N] divided by -w^2 is a receptance [m/N], which is the
%   displacement/force FRF assumed by the textbook peak-picking equations
%   (Schmitz & Smith, Machining Dynamics 2nd ed., Eqs. 2.58-2.63).
%   W must be > 0 (drop the f = 0 bin first).

w = reshape(w, size(H_acc));
if any(w <= 0)
    error('accel_to_receptance:w0', 'w must be > 0; drop the f = 0 point before converting.');
end
H_rec = H_acc ./ (-(w.^2));
end
