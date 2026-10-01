function plot_fig221(f_Hz, G, used, modes, ttl, outfile, units, Tinfo)
%PLOT_FIG221 Stacked Re/Im plot with peak-picking points, after textbook Fig. 2.21.
%   PLOT_FIG221(F_HZ, G, USED, MODES, TTL, OUTFILE, UNITS) plots Re(G) and
%   Im(G) versus frequency [Hz], marks points 1..2n (Im minima and Re
%   max/min per mode), draws w_n / w_a / w_b as vertical dotted lines and the
%   A, B, ... peak values as horizontal dashed lines, and saves a PNG.
%   USED is the coherence mask (unused points are drawn as grey dots).
%   UNITS is the y-axis unit string.
%   PLOT_FIG221(..., TINFO) adds a |T| panel on top (transmissibility):
%   TINFO fields f_Hz, magT, used, peaks (from TRANSMISSIBILITY_ANALYSIS),
%   f_iso_Hz.

if nargin < 8, Tinfo = []; end
has_T = ~isempty(Tinfo);
nrow = 2 + has_T;
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 300*nrow]);
f_Hz = f_Hz(:); G = G(:); used = logical(used(:));
xl = [min(f_Hz) max(f_Hz)];
grey = [0.6 0.6 0.6];
blue = [0 0.2 0.8];

if has_T
    ax0 = subplot(nrow, 1, 1); hold on; box on
    semilogy(Tinfo.f_Hz, Tinfo.magT, '-', 'Color', blue, 'LineWidth', 1);
    if any(~Tinfo.used)
        semilogy(Tinfo.f_Hz(~Tinfo.used), Tinfo.magT(~Tinfo.used), '.', 'Color', grey, 'MarkerSize', 6);
    end
    set(ax0, 'YScale', 'log');
    plot(xl, [1 1], 'k:');
    for q = 1:numel(Tinfo.peaks)
        pk = Tinfo.peaks(q);
        plot(pk.f_n_Hz, pk.T_peak, 'ro', 'MarkerSize', 7, 'LineWidth', 1.5);
        text(pk.f_n_Hz, pk.T_peak*1.25, sprintf('f_n=%.1f Hz', pk.f_n_Hz), 'HorizontalAlignment', 'center', 'FontSize', 9);
        if isfinite(pk.hp_w1) && isfinite(pk.hp_w2)
            lv = pk.T_peak/sqrt(2);
            plot([pk.hp_w1 pk.hp_w2]/(2*pi), [lv lv], 'm-', 'LineWidth', 1.5);
        end
    end
    if isfinite(Tinfo.f_iso_Hz)
        plot(Tinfo.f_iso_Hz, 1, 'ks', 'MarkerSize', 7, 'LineWidth', 1.5);
        text(Tinfo.f_iso_Hz, 0.7, sprintf('|T|<1 from %.1f Hz', Tinfo.f_iso_Hz), 'FontSize', 9);
    end
    xlim(xl); ylabel('|T| [-]');
    title(ttl, 'Interpreter', 'none');
end

% ---- Real part ----
subplot(nrow, 1, 1 + has_T); hold on; box on
plot(f_Hz(used), real(G(used)), '-', 'Color', blue, 'LineWidth', 1);
if any(~used), plot(f_Hz(~used), real(G(~used)), '.', 'Color', grey, 'MarkerSize', 6); end
plot(xl, [0 0], 'k:');
yr = real(G(used)); ylr = [min(yr) max(yr)]; ylr = ylr + [-0.15 0.15]*max(diff(ylr), eps);
for q = 1:numel(modes)
    M = modes(q);
    plot([M.f_a_Hz M.f_a_Hz], ylr, 'r:');
    plot([M.f_b_Hz M.f_b_Hz], ylr, 'r:');
    plot([M.f_n_Hz M.f_n_Hz], ylr, 'k--');
    plot(M.f_a_Hz, M.Re_a, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
    plot(M.f_b_Hz, M.Re_b, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
    text(M.f_a_Hz, ylr(1) + 0.05*diff(ylr), sprintf('%d', M.label_a), 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
    text(M.f_b_Hz, ylr(1) + 0.05*diff(ylr), sprintf(' %d', M.label_b), 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
end
xlim(xl); ylim(ylr);
ylabel(['Real [' units ']']);
if ~has_T, title(ttl, 'Interpreter', 'none'); end

% ---- Imaginary part ----
subplot(nrow, 1, 2 + has_T); hold on; box on
plot(f_Hz(used), imag(G(used)), '-', 'Color', blue, 'LineWidth', 1);
if any(~used), plot(f_Hz(~used), imag(G(~used)), '.', 'Color', grey, 'MarkerSize', 6); end
plot(xl, [0 0], 'k:');
yi = imag(G(used)); yli = [min(yi) max(yi)]; yli = yli + [-0.15 0.15]*max(diff(yli), eps);
leg = {};
for q = 1:numel(modes)
    M = modes(q);
    plot([M.f_n_Hz M.f_n_Hz], yli, 'k--');
    plot([xl(1) M.f_n_Hz], [M.Im_peak M.Im_peak], 'k-.');
    plot(M.f_n_Hz, M.Im_peak, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
    text(M.f_n_Hz, yli(2) - 0.08*diff(yli), sprintf('%d', M.label_n), 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(xl(1) + 0.01*diff(xl), M.Im_peak, sprintf('%s', M.peak_name), 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
    leg{end+1} = sprintf('mode %d: f_n=%.2f Hz (\\omega_n=%.1f rad/s), \\zeta=%.4f, f_a=%.2f, f_b=%.2f Hz', ...
        q, M.f_n_Hz, M.w_n, M.zeta, M.f_a_Hz, M.f_b_Hz); %#ok<AGROW>
end
xlim(xl); ylim(yli);
ylabel(['Imag [' units ']']);
xlabel('Frequency [Hz]');
if ~isempty(leg)
    text(xl(1) + 0.35*diff(xl), yli(1) + 0.08*diff(yli)*numel(leg), leg, ...
        'FontSize', 8, 'VerticalAlignment', 'bottom', 'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7]);
end

print(fig, outfile, '-dpng', '-r110');
close(fig);
end
