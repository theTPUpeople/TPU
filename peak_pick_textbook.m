function [modes, cand] = peak_pick_textbook(w, G, opts)
%PEAK_PICK_TEXTBOOK Peak-picking modal fit, Schmitz & Smith Sec. 2.5.1.
%   [MODES, CAND] = PEAK_PICK_TEXTBOOK(W, G, OPTS) identifies up to
%   OPTS.n_modes modes from a direct receptance G(w) [m/N] given in the
%   textbook sign convention (Re > 0 below resonance, Im < 0 at resonance;
%   Fig. 2.21). W is in rad/s, sorted ascending.
%
%   Textbook notation (2 modes, Fig. 2.21):
%     points 1, 2  minima of Im(G)            -> w_n1, w_n2
%     points 3, 4  max / min of Re(G), mode 1 -> w3, w4
%     points 5, 6  max / min of Re(G), mode 2 -> w5, w6
%     A, B         (negative) Im(G) at points 1, 2
%   Eq. 2.58/2.59  zeta_q = (w4 - w3) / (2 w_n1)       [(w6 - w5)/(2 w_n2)]
%   Eq. 2.60/2.61  k_q    = -1 / (2 zeta_q A)          [-1/(2 zeta_q2 B)]
%   Eq. 2.62       m_q    = k_q / w_n^2
%   Eq. 2.63       c_q    = 2 zeta_q sqrt(k_q m_q)
%   For n modes the labels generalise to: Im point q, Re points n+2q-1, n+2q.
%   Here w_a = lower Re point (w3, w5, ...), w_b = upper (w4, w6, ...).
%
%   OPTS fields:
%     n_modes            max number of modes to keep (most prominent Im minima)
%     min_prominence     absolute prominence threshold on -Im(G) (same units
%                        as G); [] -> min_prominence_rel * max(-Im(G))
%     min_prominence_rel relative threshold used when min_prominence is []
%     real_search_frac   Re extrema searched within w_n*(1 -/+ frac), clipped
%                        at the midpoint to neighbouring picked modes
%     peak_refine        how the picked bins are refined to sub-bin accuracy:
%       'sdof'  (default) local least-squares parabola in the coordinates in
%               which the SDOF curves are symmetric about their extremum:
%               Im vs x = 1-(w/w_n)^2 over |x| <= 0.35*a, and Re vs
%               u = ln|1-(w/w_n)^2| over |u - u_bin| <= ln 2, where
%               a = 2*zeta from the raw bins. For an SDOF,
%               Re = 1/(2 a k cosh(u - ln a)) and Im = -a/(k (x^2 + a^2)), so
%               the vertex is not biased by the peak asymmetry and a constant
%               residual from other modes does not move it.
%       'none'  raw bins (frequency resolution limits zeta)
%       h       (integer) parabola in w over +/- h bins
%
%     allow_inverted     true (default): Im maxima are candidates too. An
%                        inverted mode (Im > 0 at resonance) cannot occur in a
%                        passive driving-point receptance; it appears in a cross
%                        / non-collocated FRF (textbook Fig. 2.20). Its bracket
%                        is mirrored (Re min below, Re max above); w_n and zeta
%                        are reported, k_q, m_q, c_q are NaN.
%
%   Mode acceptance: candidates above the prominence threshold are tried in
%   order of prominence; one is kept only if, with it added, every kept mode
%   still has a resolved Re bracket (the Re extremes lie strictly inside the
%   mode's window and on opposite sides of w_n). Rejected candidates and the
%   reason are listed in CAND.status (band-edge artefacts, closely spaced
%   modes that peak picking cannot separate).
%
%   MODES is a struct array sorted by frequency (raw-bin values are kept in
%   the *_raw fields). CAND lists every Im extremum found with its prominence
%   and status.

w = w(:); G = G(:);
re = real(G); im = imag(G);
if ~isfield(opts, 'peak_refine'), opts.peak_refine = 'sdof'; end
if ~isfield(opts, 'allow_inverted'), opts.allow_inverted = true; end

% candidates: minima of Im (normal modes) and, optionally, maxima of Im
% (inverted modes, only possible in a cross / non-collocated FRF)
[ip, prom] = find_peaks_prominence(-im);
keep = im(ip) < 0; ip = ip(keep); prom = prom(keep);
sg = ones(size(ip));
if opts.allow_inverted
    [ipx, promx] = find_peaks_prominence(im);
    keep = im(ipx) > 0;
    ip = [ip; ipx(keep)]; prom = [prom; promx(keep)]; sg = [sg; -ones(sum(keep), 1)];
end
if isempty(opts.min_prominence)
    thr = opts.min_prominence_rel * max(abs(im));
else
    thr = opts.min_prominence;
end
pass = prom >= thr;
status = repmat({'below threshold'}, size(ip));

% accept candidates in order of prominence; a candidate is kept only if,
% with it added, every accepted mode still has a resolved Re bracket
% (Re extremes strictly inside its window, on opposite sides of w_n)
[~, order] = sort(prom, 'descend');
acc = zeros(0, 1);
notes = {};
for c = order(:)'
    if ~pass(c), continue, end
    if numel(acc) >= opts.n_modes
        status{c} = 'passed, not used (n_modes reached)';
        continue
    end
    trial = [acc; c];
    [~, o] = sort(ip(trial)); trial = trial(o);     % by frequency (ip mixes minima and maxima)
    [~, ok, why] = build_modes(w, re, im, ip(trial), sg(trial), opts, false);
    if all(ok)
        acc = trial;
        status{c} = 'USED';
    else
        k = find(trial == c);
        if ~ok(k)
            status{c} = ['rejected: ' why{k}];
        else
            other = trial(find(~ok, 1));
            status{c} = sprintf('rejected: unresolvable next to the mode at %.2f Hz (closely spaced)', w(ip(other))/(2*pi));
            notes{end+1} = [other, c]; %#ok<AGROW>
        end
    end
end

cand.w = w(ip);
cand.f_Hz = w(ip)/(2*pi);
cand.Im = im(ip);
cand.inverted = sg < 0;
cand.prominence = prom;
cand.threshold = thr;
cand.passed = pass;
cand.used = ismember((1:numel(ip))', acc);
cand.status = status;

modes = build_modes(w, re, im, ip(acc), sg(acc), opts, true);
for j = 1:numel(notes)
    q = find(acc == notes{j}(1));
    if ~isempty(q)
        modes(q).flags{end+1} = sprintf('a weaker Im extremum at %.2f Hz could not be separated (closely spaced modes?)', ...
            w(ip(notes{j}(2)))/(2*pi));
    end
end

% well-separated-mode check (textbook assumption)
for q = 1:numel(modes)-1
    if modes(q).w_b >= modes(q+1).w_a
        msg = sprintf('modes %d and %d overlap: Re brackets [%.1f %.1f] and [%.1f %.1f] Hz intersect', ...
            q, q+1, modes(q).f_a_Hz, modes(q).f_b_Hz, modes(q+1).f_a_Hz, modes(q+1).f_b_Hz);
        modes(q).flags{end+1} = msg;
        modes(q+1).flags{end+1} = msg;
    end
end
end

% -------------------------------------------------------------------------
function [modes, ok, why] = build_modes(w, re, im, idx, sg, opts, full)
% Modes at Im-extremum indices idx (sorted by frequency), sg = +1 normal /
% -1 inverted. ok(q) = Re bracket resolved. With full = false only the
% bracket check is done (used while selecting candidates).
nm = numel(idx);
n = numel(w);
ok = true(nm, 1);
why = repmat({''}, nm, 1);
modes = struct([]);
for q = 1:nm
    i_n = idx(q);
    s = sg(q);
    % search window: w_n*(1 -/+ frac), clipped at midpoints to neighbours
    lo = w(i_n)*(1 - opts.real_search_frac);
    hi = w(i_n)*(1 + opts.real_search_frac);
    if q > 1,  lo = max(lo, 0.5*(w(idx(q-1)) + w(i_n))); end
    if q < nm, hi = min(hi, 0.5*(w(i_n) + w(idx(q+1)))); end
    win = find(w >= lo & w <= hi);
    if numel(win) < 5
        ok(q) = false; why{q} = 'search window has fewer than 5 points';
        if full, error('peak_pick_textbook:window', 'accepted mode with an empty search window'); end
        continue
    end
    ilo = win(1); ihi = win(end);
    % points 3/4 (5/6), raw bins: max of s*Re below w_n, min of s*Re above
    % (an inverted mode mirrors the textbook shape, so work on s*G)
    rs = s*re; is_ = s*im;
    ia_set = (ilo:i_n)';
    ib_set = (i_n:ihi)';
    [~, k] = max(rs(ia_set)); i_a = ia_set(k);
    [~, k] = min(rs(ib_set)); i_b = ib_set(k);
    if ~(ilo < i_a && i_a < i_n)
        ok(q) = false;
        why{q} = 'no Re maximum below the Im extremum inside its window (bracket not resolved)';
    elseif ~(i_n < i_b && i_b < ihi)
        ok(q) = false;
        why{q} = 'no Re minimum above the Im extremum inside its window (bracket not resolved)';
    end
    if ~full, continue, end

    M = struct();
    M.q = q;
    M.label_n = q;
    M.label_a = nm + 2*q - 1;
    M.label_b = nm + 2*q;
    M.peak_name = char('A' + q - 1);
    M.inverted = s < 0;
    M.flags = {};
    M.win_lo = w(ilo); M.win_hi = w(ihi);
    if rs(i_a) <= 0
        M.flags{end+1} = sprintf('Re at point %d does not have the textbook sign (offset from other modes?)', M.label_a);
    end
    if rs(i_b) >= 0
        M.flags{end+1} = sprintf('Re at point %d does not have the textbook sign (offset from other modes?)', M.label_b);
    end
    M.i_n = i_n; M.i_a = i_a; M.i_b = i_b;
    M.w_n_raw = w(i_n); M.w_a_raw = w(i_a); M.w_b_raw = w(i_b);
    M.Im_peak_raw = im(i_n);
    M.zeta_raw = (w(i_b) - w(i_a))/(2*w(i_n));
    M.m_q_raw = -1/(2*M.zeta_raw*im(i_n))/w(i_n)^2;

    % sub-bin refinement of points 1-6 (on s*G, then mapped back)
    w_n = w(i_n); y_n = is_(i_n);
    w_a = w(i_a); y_a = rs(i_a);
    w_b = w(i_b); y_b = rs(i_b);
    rf = opts.peak_refine;
    if ischar(rf) && strcmpi(rf, 'sdof')
        a = 2*M.zeta_raw;
        [w_n, y_n, ok1] = refine_im_sdof(w, is_, i_n, ilo, ihi, a);
        [w_a, y_a, ok2] = refine_re_sdof(w, rs, i_a, ilo, ihi, w_n, +1);
        [w_b, y_b, ok3] = refine_re_sdof(w, rs, i_b, ilo, ihi, w_n, -1);
        lab = [M.label_n M.label_a M.label_b];
        for b = find(~[ok1 ok2 ok3])
            M.flags{end+1} = sprintf('point %d: SDOF-symmetric refinement not possible, raw bin used', lab(b));
        end
    elseif isnumeric(rf) && rf > 0
        [w_n, y_n] = refine_quad(w, is_, i_n, rf, ilo, ihi, 'min');
        [w_a, y_a] = refine_quad(w, rs, i_a, rf, ilo, i_n, 'max');
        [w_b, y_b] = refine_quad(w, rs, i_b, rf, i_n, ihi, 'min');
    end
    M.w_n = w_n; M.Im_peak = s*y_n;
    M.w_a = w_a; M.Re_a = s*y_a;
    M.w_b = w_b; M.Re_b = s*y_b;

    M.f_n_Hz = M.w_n/(2*pi);
    M.f_a_Hz = M.w_a/(2*pi);
    M.f_b_Hz = M.w_b/(2*pi);

    % Eq. 2.58 / 2.59
    M.zeta = (M.w_b - M.w_a)/(2*M.w_n);
    if ~(M.zeta > 0), M.zeta = NaN; end
    if ~M.inverted
        % Eq. 2.60 / 2.61 (A < 0)
        M.k_q = -1/(2*M.zeta*M.Im_peak);
        % Eq. 2.62
        M.m_q = M.k_q/M.w_n^2;
        % Eq. 2.63
        M.c_q = 2*M.zeta*sqrt(M.k_q*M.m_q);
    else
        % Im > 0 at resonance: impossible in a passive driving-point
        % receptance; Eq. 2.60 would give k_q < 0
        M.k_q = NaN; M.m_q = NaN; M.c_q = NaN; M.m_q_raw = NaN;
        M.flags{end+1} = sprintf(['inverted mode (Im = %+.3g > 0 at resonance): only possible in a cross / ', ...
            'non-collocated FRF; w_n and zeta kept, k_q m_q c_q not defined'], M.Im_peak);
    end

    if isempty(modes), modes = M; else, modes(q) = M; end
end
end

% -------------------------------------------------------------------------
function [wv, yv, ok] = refine_im_sdof(w, im, i, ilo, ihi, a)
% Im = -a/(k (x^2 + a^2)) is even in x = 1-(w/w_n)^2: parabola in x over |x| <= 0.35 a.
wv = w(i); yv = im(i); ok = false;
if ~(a > 0), return, end
j = (ilo:ihi)';
t = (1 - (w(j)/w(i)).^2)/a;
j = j(abs(t) <= 0.35); t = (1 - (w(j)/w(i)).^2)/a;
if numel(j) < 3, return, end
p = [t.^2 t ones(size(t))] \ im(j);
if p(1) <= 0, return, end
tv = -p(2)/(2*p(1));
if tv < min(t) || tv > max(t), return, end
wv = w(i)*sqrt(1 - a*tv);
yv = p(1)*tv^2 + p(2)*tv + p(3);
ok = true;
end

function [wv, yv, ok] = refine_re_sdof(w, re, i, ilo, ihi, wn, side)
% side = +1: Re max below w_n (x = 1-(w/w_n)^2 > 0); side = -1: Re min above w_n.
% Re is even in u = ln|x| about ln a: parabola in u over |u - u_i| <= ln 2.
wv = w(i); yv = re(i); ok = false;
j = (ilo:ihi)';
x = side*(1 - (w(j)/wn).^2);
xi = side*(1 - (w(i)/wn)^2);
if xi <= 0, return, end
keep = x > 0;
j = j(keep); u = log(x(keep)) - log(xi);
keep = abs(u) <= log(2);
j = j(keep); u = u(keep);
if numel(j) < 4, return, end
p = [u.^2 u ones(size(u))] \ re(j);
if side*p(1) >= 0, return, end           % max needs p(1) < 0, min needs p(1) > 0
uv = -p(2)/(2*p(1));
if uv < min(u) || uv > max(u), return, end
wv = wn*sqrt(1 - side*xi*exp(uv));
yv = p(1)*uv^2 + p(2)*uv + p(3);
ok = true;
end

function [x0, y0] = refine_quad(x, y, i, h, ilo, ihi, kind)
% Local least-squares parabola in w over bins i-h..i+h (clipped to [ilo ihi]).
x0 = x(i); y0 = y(i);
j = (max(ilo, i-h):min(ihi, i+h))';
if numel(j) < 3, return, end
s = x(j(end)) - x(j(1));
if s <= 0, return, end
u = (x(j) - x(i))/s;
p = [u.^2, u, ones(size(u))] \ y(j);
if (strcmp(kind, 'min') && p(1) <= 0) || (strcmp(kind, 'max') && p(1) >= 0), return, end
uv = -p(2)/(2*p(1));
if uv < min(u) || uv > max(u), return, end
x0 = x(i) + uv*s;
y0 = p(1)*uv^2 + p(2)*uv + p(3);
end
