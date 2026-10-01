function C = read_quickdaq_config(xmlfile)
%READ_QUICKDAQ_CONFIG Acquisition settings from a QuickDAQ "... config.xml".
%   C = READ_QUICKDAQ_CONFIG(XMLFILE) extracts, by tag name, the enabled
%   input channels (name, unit label, mV/unit sensitivity) and the FFT
%   settings (FFT size, sample rate, averages, windows). Used only to
%   document the acquisition in the run log; no value here is applied to the
%   FRF. Returns [] if the file does not exist.

C = [];
if exist(xmlfile, 'file') ~= 2, return, end
s = fileread(xmlfile);

C.file = xmlfile;
C.FFTSize        = num(tag(s, 'FFTSize'));
C.SampleRate     = num(tag(s, 'SampleRate'));
C.NumberOfAverages = num(tag(s, 'NumberOfAverages'));
C.WindowTypeResp = tag(s, 'WindowTypeResp');
C.WindowTypeRef  = tag(s, 'WindowTypeRef');
C.ExpWindowFinalValue = num(tag(s, 'ExpWindowFinalValue'));
C.ForceWindowEndPercent = num(tag(s, 'ForceWindowEndPercent'));

% channel list: each <ListItem> ... </ListItem> with Name / Enabled / EditUnitLabel / MilliVoltsPerUnit
items = regexp(s, '<ListItem>([\s\S]*?)</ListItem>', 'tokens');
C.channels = struct('name', {}, 'unit', {}, 'mV_per_unit', {}, 'enabled', {});
for i = 1:numel(items)
    b = items{i}{1};
    nm = tag(b, 'Name');
    if isempty(nm), continue, end
    ch.name = nm;
    ch.unit = tag(b, 'EditUnitLabel');
    ch.mV_per_unit = num(tag(b, 'MilliVoltsPerUnit'));
    ch.enabled = strcmpi(tag(b, 'Enabled'), 'true');
    C.channels(end+1) = ch;
end
end

function v = tag(s, t)
m = regexp(s, ['<' t '>([\s\S]*?)</' t '>'], 'tokens', 'once');
if isempty(m), v = ''; else, v = strtrim(m{1}); end
end

function x = num(v)
if isempty(v), x = NaN; else, x = str2double(v); end
end
