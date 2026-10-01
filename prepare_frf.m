function P = prepare_frf(f_Hz, H, coh, opts)
%PREPARE_FRF Clean a measured FRF before peak picking.
%   P = PREPARE_FRF(F_HZ, H, COH, OPTS) sorts by frequency, drops NaN and
%   duplicate-frequency points, drops f <= 0, keeps OPTS.f_band_Hz = [lo hi],
%   optionally smooths the complex FRF over OPTS.smooth_points (odd, 1 = off)
%   with OPTS.smooth_method = 'box' (centred moving average, default) or
%   'sgolay' (quadratic Savitzky-Golay: preserves peak curvature, so it
%   broadens a resonance far less for similar noise reduction), and masks
%   points with coherence < OPTS.min_coherence.
%
%   P.f_Hz, P.w, P.H, P.coh   points used for identification (w = 2*pi*f [rad/s])
%   P.band_f_Hz, P.band_H, P.band_used   all in-band points and the usage mask
%                                        (for plotting the masked points)
%   P.notes                   cell array of what was done (for the run log)

f_Hz = f_Hz(:); H = H(:);
has_coh = ~isempty(coh);
if has_coh, coh = coh(:); else, coh = nan(size(f_Hz)); end
notes = {};

ok = isfinite(f_Hz) & isfinite(real(H)) & isfinite(imag(H));
if has_coh, ok = ok & isfinite(coh); end
if any(~ok), notes{end+1} = sprintf('dropped %d NaN/Inf points', sum(~ok)); end
f_Hz = f_Hz(ok); H = H(ok); coh = coh(ok);

if any(diff(f_Hz) < 0)
    notes{end+1} = 'frequency vector was not sorted; sorted ascending';
end
[f_Hz, is] = sort(f_Hz); H = H(is); coh = coh(is);

dup = [false; diff(f_Hz) == 0];
if any(dup)
    notes{end+1} = sprintf('dropped %d duplicate-frequency points (kept first)', sum(dup));
    f_Hz = f_Hz(~dup); H = H(~dup); coh = coh(~dup);
end

nz = f_Hz <= 0;
if any(nz)
    notes{end+1} = sprintf('dropped %d point(s) at f <= 0 Hz', sum(nz));
    f_Hz = f_Hz(~nz); H = H(~nz); coh = coh(~nz);
end

band = f_Hz >= opts.f_band_Hz(1) & f_Hz <= opts.f_band_Hz(2);
f_Hz = f_Hz(band); H = H(band); coh = coh(band);
notes{end+1} = sprintf('analysis band %.2f-%.2f Hz: %d points', ...
    opts.f_band_Hz(1), opts.f_band_Hz(2), numel(f_Hz));
if numel(f_Hz) < 10
    error('prepare_frf:band', 'fewer than 10 points in the analysis band.');
end

ns = opts.smooth_points;
if ~isfield(opts, 'smooth_method'), opts.smooth_method = 'box'; end
if ns > 1
    if mod(ns, 2) == 0, ns = ns + 1; end
    m = (ns - 1)/2;
    switch lower(opts.smooth_method)
        case 'box'
            % centred (zero-phase) moving average of the complex FRF; the original
            % script used a causal filter() which shifts peaks by (ns-1)/2 bins
            k = ones(ns, 1)/ns;
            H = conv(H, k, 'same') ./ conv(ones(size(H)), k, 'same');   % edge bins renormalised
            notes{end+1} = sprintf('centred moving average over %d points applied to complex FRF', ns);
        case 'sgolay'
            % quadratic Savitzky-Golay smoothing weights (closed form, sum to 1,
            % reproduce any quadratic exactly); the m edge bins are left unsmoothed
            j = (-m:m)';
            k = (3*(3*m^2 + 3*m - 1) - 15*j.^2)/((2*m + 1)*(4*m^2 + 4*m - 3));
            Hs = conv(H, k, 'same');
            Hs([1:m, end-m+1:end]) = H([1:m, end-m+1:end]);
            H = Hs;
            notes{end+1} = sprintf('quadratic Savitzky-Golay smoothing over %d points applied to complex FRF', ns);
        otherwise
            error('prepare_frf:smooth', 'smooth_method must be box or sgolay');
    end
end

used = true(size(f_Hz));
if has_coh && opts.min_coherence > 0
    used = coh >= opts.min_coherence;
    if any(~used)
        notes{end+1} = sprintf('masked %d of %d in-band points with coherence < %.2f', ...
            sum(~used), numel(used), opts.min_coherence);
    end
elseif ~has_coh
    notes{end+1} = 'no coherence column: no coherence screening possible';
end

P.band_f_Hz = f_Hz;
P.band_H = H;
P.band_coh = coh;
P.band_used = used;
P.f_Hz = f_Hz(used);
P.w = 2*pi*P.f_Hz;
P.H = H(used);
P.coh = coh(used);
P.has_coh = has_coh;
P.notes = notes;
end
