%%  Load EMG and LFP data from file
%%
%%
%%  Input:
%%          file        filename of datafile until K, e.g.,
%%                      file='EMG182002834K'
%%          part        part of data to be picked out for analysis in
%%                      seconds, e.g., part=[7 29];
%%  Output:
%%          data        data matrix n x m with n number of channels and m
%%                      number of data points, channels 1-5 LFP, 6-8 EMG
%%          t           time vector in seconds
%%
%% Author:  M. von Papen
%% Date:    28.10.2014

function [data, t]=load_LFP_EMG(file,part)

pathStr = 'C:\Dokumente und Einstellungen\maierf\Eigene Dateien\Biosig-1.97\';
sf=2500; %sampling frequency
nchannels=8; %Number of channels, 1-5 are LFP, 6-8 EMG


for i=[1:8]
    
    fid = fopen([file num2str(i) '.dat'], 'r');  %opens file
    fseek(fid, 602, 'bof');  %sets Pointer to Position 602 (bof)
    data(i,:) = fread(fid, 'ubit12', 4);  %loads data
    fclose(fid);    %closes file

end

if nargin>1
    data = data(:,part(1)*sf:part(2)*sf);
    t=[0:length(data)-1]/sf+part(1);
else
    data = (data-2048)*0.4/4096;  %transforms Bits to Volt
    t=[0:length(data)-1]/sf;
end

data=data';

clear xlab ylab1 ylab2 i lim h