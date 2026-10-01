function write_quickdaq_csv(file, f_Hz, H, coh, layout)
%WRITE_QUICKDAQ_CSV Write an FRF in the QuickDAQ CSV export layout (for tests).
%   WRITE_QUICKDAQ_CSV(FILE, F_HZ, H, COH, LAYOUT) writes magnitude and
%   phase [deg] of the complex FRF H plus coherence COH using the same
%   9-line header as the exports in this repo:
%     LAYOUT = 'D'  disturbance rejection: Frequency, FRF mag, FRF phase, Coherence
%                   (Y units m/s^2/N)
%     LAYOUT = 'T'  transmissibility:      Frequency, Coherence, FRF mag, FRF phase
%                   (Y units m/s^2/m/s^2)

d = fileparts(file);
if ~isempty(d) && exist(d, 'dir') ~= 7, mkdir(d); end
f = fopen(file, 'w');
fprintf(f, 'QuickDAQ Data\r\n1/1/0001 12:00:00 AM\r\nNotes:,,\r\nSample Rate:,1600.00938, Hz,\r\n');
mag = abs(H(:)); ph = angle(H(:))*180/pi;
switch upper(layout)
    case 'D'
        fprintf(f, 'Measurement Type,FRF,,Coherence,\r\n');
        fprintf(f, 'Channel Name,SYNTH.Floor Accelerometer-FRF,,SYNTH.Acc-Coherence,\r\n');
        fprintf(f, 'X Axis Units,Hz,,Hz,\r\n');
        fprintf(f, 'Y Axis Units,m/s^2/N,,Hammer-Floor Accelerometer Coherence,\r\n');
        fprintf(f, 'Frequency,Magnitude,Phase,Magnitude,\r\n');
        for i = 1:numel(f_Hz)
            fprintf(f, '%.15g, %.10E,%.10f, %.10f,\r\n', f_Hz(i), mag(i), ph(i), coh(i));
        end
    case 'T'
        fprintf(f, 'Measurement Type,Coherence,FRF,,\r\n');
        fprintf(f, 'Channel Name,SYNTH.Output Acc-Coherence,SYNTH.Output Acc-FRF,,\r\n');
        fprintf(f, 'X Axis Units,Hz,Hz,,\r\n');
        fprintf(f, 'Y Axis Units,Input Acc-Output Acc Coherence,m/s^2/m/s^2,,\r\n');
        fprintf(f, 'Frequency,Magnitude,Magnitude,Phase,\r\n');
        for i = 1:numel(f_Hz)
            fprintf(f, '%.15g, %.10f, %.10E, %.10f,\r\n', f_Hz(i), coh(i), mag(i), ph(i));
        end
    otherwise
        fclose(f);
        error('write_quickdaq_csv:layout', 'layout must be D or T');
end
fclose(f);
end
