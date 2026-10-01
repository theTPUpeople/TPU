function [x] = modalfit(freq,real_filtered,imag_filtered,x0)
%MODALFIT: This function optimizes the input modal parameters, solved via
%peak picking, for a given frequency response function (FRF).
%   x0 = [wn_1 k_1 zeta_1 wn_2 k_2 zeta_2 ...]  (wn in rad/s, k in N/m)
%   freq in Hz; real_filtered / imag_filtered = receptance [m/N] in the
%   textbook sign convention (Im < 0 at resonance).

% Force row vectors: with column inputs, Q_R_total (row) - X_F (column)
% silently broadcast to an N x N matrix.
freq = freq(:).';
real_filtered = real_filtered(:).';
imag_filtered = imag_filtered(:).';
x0 = x0(:).';

if exist('optimoptions', 'file') || exist('optimoptions', 'builtin')
    options = optimoptions('lsqnonlin','Display','iter','FunValCheck','on','TolFun',1e-15);
else
    % GNU Octave (optim package) has lsqnonlin but no optimoptions
    options = optimset('Display','iter','FunValCheck','on','TolFun',1e-15);
end
d = size(x0);
numrows = d(1);
numcols = d(2);
num_modes = numcols/3;
lb = zeros(numrows, numcols);
ub = ones(numrows, numcols);
for cnt1 = 1:num_modes
    lb(1, (1 + (cnt1 - 1)*3)) = x0(1, (1 + (cnt1 - 1)*3)) - 1;
    ub(1, (1 + (cnt1 - 1)*3)) = x0(1, (1 + (cnt1 - 1)*3)) + 1;
    ub(1, (2 + (cnt1 - 1)*3)) = 1e10;
    ub(1, (3 + (cnt1 - 1)*3)) = 0.2;
end

x = lsqnonlin(@obj_fct,x0,lb,ub,options);
x = reshape(x, size(x0));           % Octave's lsqnonlin returns a column
    
    function [y] = obj_fct(x)
        x = x(:).';                     % Octave passes x as a column; indexing below is x(1, k)
        w = freq*2*pi;                  % frequency vector, [Hz]
        Q_R_total = zeros(1,length(freq));
        
        for cnt2 = 1:num_modes
            r = w/x(1, (1 + (cnt2 - 1)*3));
            Q_R(cnt2,:) = (1/x(1, (2 + (cnt2 - 1)*3)))*(((1-r.^2)-1i*(2*x(1, (3 + (cnt2 - 1)*3))*r))./((1-r.^2).^2+(2*x(1, (3 + (cnt2 - 1)*3))*r).^2));
            Q_R_total = Q_R_total + Q_R(cnt2,:);
        end
        
        X_F_Mag = sqrt(real_filtered.^2 + imag_filtered.^2);
        Q_R_Mag = sqrt(real(Q_R_total).^2 + imag(Q_R_total).^2);
        Q_R_Real = real(Q_R_total);
        Q_R_Imag = imag(Q_R_total);
        X_F = [X_F_Mag; real_filtered; imag_filtered];
        Q_R = [Q_R_Mag; Q_R_Real; Q_R_Imag];
        y = Q_R - X_F;
    end
end

