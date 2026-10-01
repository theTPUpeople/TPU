function D = read_quickdaq_csv(file, phase_units)
%READ_QUICKDAQ_CSV Read a QuickDAQ frequency-domain CSV export (FRF + coherence).
%   D = READ_QUICKDAQ_CSV(FILE) parses the QuickDAQ header block and locates
%   the FRF and coherence columns from the "Measurement Type" row and the
%   column-label row ("Frequency, Magnitude, Phase, ..."). Column positions are
%   NOT assumed: the disturbance-rejection exports are (FRF mag, FRF phase,
%   Coherence) while the transmissibility exports are (Coherence, FRF mag,
%   FRF phase), and some files carry extra spreadsheet columns after a blank.
%
%   D = READ_QUICKDAQ_CSV(FILE, PHASE_UNITS) with PHASE_UNITS = 'deg' | 'rad'
%   overrides the phase-unit inference ('auto', the default, infers degrees
%   when any |phase| > 2*pi and stops with an error if it cannot tell).
%
%   Output fields:
%     f_Hz         frequency column [Hz] (as stored, not yet cleaned)
%     H            complex FRF (column)
%     coh          coherence (column) or [] if the file has none
%     x_units      X axis units of the FRF column (must be 'Hz')
%     y_units      Y axis units of the FRF column, as written in the file
%     frf_kind     'accelerance' | 'receptance' | 'mobility' |
%                  'transmissibility' | 'unknown'   (from y_units)
%     storage      'mag-phase' | 're-im'
%     phase_units  'deg' | 'rad' | ''  (empty for re-im storage)
%     frf_channel, coh_channel   channel names from the header
%     ignored_cols labels of columns that were not part of the FRF/coherence
%     notes        cell array of strings describing what was inferred

if nargin < 2 || isempty(phase_units), phase_units = 'auto'; end

txt = fileread(file);
lines = regexp(txt, '\r?\n', 'split');

% ---- locate header rows by their first cell, not by fixed line numbers ----
r_mtype = 0; r_chan = 0; r_xu = 0; r_yu = 0; r_lab = 0;
for i = 1:min(numel(lines), 50)
    c = split_csv_line(lines{i});
    switch lower(strtrim(c{1}))
        case 'measurement type', r_mtype = i;
        case 'channel name',     r_chan  = i;
        case 'x axis units',     r_xu    = i;
        case 'y axis units',     r_yu    = i;
        case 'frequency',        r_lab   = i; break
    end
end
if any([r_mtype r_chan r_xu r_yu r_lab] == 0)
    error('read_quickdaq_csv:header', ...
        '%s: QuickDAQ header rows (Measurement Type / Channel Name / X,Y Axis Units / Frequency) not all found.', file);
end

mt  = split_csv_line(lines{r_mtype});
ch  = split_csv_line(lines{r_chan});
xu  = split_csv_line(lines{r_xu});
yu  = split_csv_line(lines{r_yu});
lab = split_csv_line(lines{r_lab});
ncol = max([numel(mt) numel(ch) numel(xu) numel(yu) numel(lab)]);
mt = pad(mt, ncol); ch = pad(ch, ncol); xu = pad(xu, ncol); yu = pad(yu, ncol); lab = pad(lab, ncol);

% ---- assign each data column to a measurement block ----
% A block starts where "Measurement Type" is non-empty and continues over the
% following columns whose type cell is empty but whose label is non-empty.
owner = cell(1, ncol); owner(:) = {''};
cur = '';
for j = 2:ncol
    if ~isempty(strtrim(mt{j}))
        cur = strtrim(mt{j});
    elseif isempty(strtrim(lab{j}))
        cur = '';                       % blank label ends the block
    end
    if ~isempty(strtrim(lab{j})), owner{j} = cur; end
end

frf_cols = find(strcmpi(owner, 'FRF'));
coh_cols = find(strcmpi(owner, 'Coherence'));
if isempty(frf_cols)
    error('read_quickdaq_csv:nofrf', '%s: no column with Measurement Type = FRF.', file);
end
ignored = find(~cellfun(@isempty, strtrim(lab)) & cellfun(@isempty, owner));
ignored = ignored(ignored > 1);

frf_labels = lower(strtrim(lab(frf_cols)));
i_mag = frf_cols(strcmp(frf_labels, 'magnitude'));
i_ph  = frf_cols(strcmp(frf_labels, 'phase'));
i_re  = frf_cols(strcmp(frf_labels, 'real'));
i_im  = frf_cols(ismember(frf_labels, {'imaginary', 'imag'}));
if numel(i_mag) == 1 && numel(i_ph) == 1
    storage = 'mag-phase';
elseif numel(i_re) == 1 && numel(i_im) == 1
    storage = 're-im';
else
    error('read_quickdaq_csv:frflayout', ...
        '%s: FRF block labels {%s} are neither Magnitude+Phase nor Real+Imaginary.', ...
        file, strjoin(lab(frf_cols), ', '));
end
i_coh = [];
if ~isempty(coh_cols)
    i_coh = coh_cols(strcmpi(strtrim(lab(coh_cols)), 'magnitude'));
    if numel(i_coh) ~= 1, i_coh = coh_cols(1); end
end

% ---- numeric block ----
data_lines = lines(r_lab+1:end);
data_lines = data_lines(~cellfun(@(s) isempty(strtrim(s)), data_lines));
need = unique([1 i_mag i_ph i_re i_im i_coh]);
M = nan(numel(data_lines), ncol);
for i = 1:numel(data_lines)
    c = split_csv_line(data_lines{i});
    n = min(numel(c), ncol);
    v = str2double(c(1:n));
    M(i, 1:n) = v;
end
M = M(:, 1:ncol);
if all(isnan(M(:, need(end))))
    error('read_quickdaq_csv:nodata', '%s: no numeric data under the header.', file);
end

D.file = file;
D.f_Hz = M(:, 1);
D.x_units = strtrim(xu{frf_cols(1)});
D.y_units = strtrim(yu{frf_cols(1)});
D.frf_channel = strtrim(ch{frf_cols(1)});
D.coh_channel = '';
D.storage = storage;
D.phase_units = '';
D.notes = {};
D.ignored_cols = strtrim(lab(ignored));

if ~strcmpi(D.x_units, 'Hz')
    error('read_quickdaq_csv:xunits', '%s: X axis units are "%s", expected Hz.', file, D.x_units);
end

if strcmp(storage, 'mag-phase')
    mag = M(:, i_mag);
    ph  = M(:, i_ph);
    switch lower(phase_units)
        case 'auto'
            if max(abs(ph)) > 2*pi + 1e-6
                D.phase_units = 'deg';
                D.notes{end+1} = sprintf('phase units inferred as degrees (max |phase| = %.1f > 2*pi)', max(abs(ph)));
            else
                error('read_quickdaq_csv:phaseunits', ...
                    '%s: cannot infer phase units (all |phase| <= 2*pi). Set cfg.phase_units to ''deg'' or ''rad''.', file);
            end
        case {'deg', 'rad'}
            D.phase_units = lower(phase_units);
            D.notes{end+1} = sprintf('phase units set by user: %s', D.phase_units);
        otherwise
            error('read_quickdaq_csv:phaseunits', 'phase_units must be auto, deg or rad.');
    end
    if strcmp(D.phase_units, 'deg'), ph = ph*pi/180; end
    D.H = mag .* exp(1i*ph);
else
    D.H = M(:, i_re) + 1i*M(:, i_im);
end

if ~isempty(i_coh)
    D.coh = M(:, i_coh);
    D.coh_channel = strtrim(ch{i_coh});
else
    D.coh = [];
end

D.frf_kind = classify_units(D.y_units);
D.notes{end+1} = sprintf('FRF columns %s (%s), coherence column %s, Y units "%s" -> %s', ...
    mat2str(frf_cols), storage, mat2str(i_coh), D.y_units, D.frf_kind);
end

% -------------------------------------------------------------------------
function kind = classify_units(u)
s = lower(regexprep(u, '\s', ''));
switch s
    case {'m/s^2/n', '(m/s^2)/n', 'm/s2/n'}
        kind = 'accelerance';
    case {'m/n'}
        kind = 'receptance';
    case {'m/s/n', '(m/s)/n'}
        kind = 'mobility';
    case {'m/s^2/m/s^2', '(m/s^2)/(m/s^2)', 'g/g'}
        kind = 'transmissibility';
    otherwise
        kind = 'unknown';
end
end

function c = split_csv_line(s)
c = regexp(s, ',', 'split');
if isempty(c), c = {''}; end
end

function c = pad(c, n)
if numel(c) < n, c(end+1:n) = {''}; end
end
