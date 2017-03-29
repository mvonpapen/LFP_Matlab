%% Load data and calculate Wavelet transformation and coherency
function [t, datf, W, CLE, WLE, data, f] = get_data ( patient, site, T, varargin )

args=struct('freq',logspace(-1,2,30),...
            'dLFP',0);
args=parseArgs(varargin,args);

f=args.freq;
nf=length(f);

%% Load patient site data
       load(['/afs/geo/usr/vonpapen/Neuro/DATA/' patient, ...
        '_S' mat2str(site) '.mat'])
% load(['D:\Uni\Neuro\DATA\' patient '_S' mat2str(site) '.mat'])
[nt, nch] = size( data );
nch=nch-1; %nch without STIM channel

%% Parameters
dt=1/2500;

%% Check time vector
t = (0:nt-1)*dt;
ti = t>=T(1) & t<=T(2);
t = t(ti);
data=data(ti,1:nch);
nt = length(t);

%% Variables
if args.dLFP==1
    for i=2:nch-3
        data(:,i)=data(:,i)-data(:,1);
    end
    data=data(:,2:end);
    nch=nch-1; %one channel less if dLFP
end

W=NaN(nf,nt,nch);
P=NaN(nf,nch);
datf=NaN(nt,nch);

for n=1:nch
	
    if n<=nch-3
        %% LFP data
        [datf(:,n), W(:,:,n), coi, P(:,n)] = procdata (data(:,n), 'freq', f);
    else
        %% EMG data
    	[datf(:,n), W(:,:,n), coi, P(:,n)] = procdata (data(:,n), ...
            'filter', [0 60; 90 110], 'rect', 1, 'freq', f);
    end
    
end

%% Coherency
[CLE, WLE] = wave_coh (W(:,:,1:nch-3), W(:,:,nch-2:nch), f, dt);