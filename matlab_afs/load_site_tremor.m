%LOAD_SITE_TREMOR  Load tremor data for a given site
%
%   [t, dat, CHlab, art, loc, pos] = load_site_tremor( site, varargin )
%
%   Loads data from DATA_tremor.mat for site and gives out time series,
%   channel labels, artefact intervals, location and position of
%   measurement. Output time series is already picked tremor interval.
%
%
%   INPUT:
%       site:       site number of measurement
%
%   OPTIONAL INPUT:
%       channel:    cell of channel names or indice to use (default all)
%       noEMG:      Disregard all EMGs (default is 0)
%       commontime: find common interval of all channels (default is 1)
%
%   OUTPUT:
%       t:          time vector if commontime, otherwise cell array
%       dat:        matrix of data if commontime, otherwise cell array
%       pat:        name of patient
%       CHlab:      names of output channels (cell array)
%       art:        artifacts for channels (cell array)
%       loc:        location of measurement (0=unknown, 1=STN, 2=ZI)
%       pos:        position along trajectory w.r.t. target (STN) in mm
%       bad:        logical indicating if channel is bad
%
%
%   Authors:
%       Michael von Papen, Felix Gerick
%
%   Date:
%       06.05.16

function [t, dat, pat, CHlab, art, loc, pos, bad] = load_site_tremor( site, varargin )


%% Arguments
args = struct('channel',[],...
              'noEMG',0, ...
              'commontime',1);
        
args = parseArgs(varargin,args);
channel = args.channel;
dt = 1/2456;


%% Load data
i           = site;
DATA_Tremor = [];
data        = [];
load('/afs/rrz/geo/project/neuro/tremor_data/DATA/DATA_tremor.mat');
pat         = DATA_Tremor(i).patient;

load(['/afs/rrz/geo/project/neuro/tremor_data/DATA/' ...
    pat '_S' mat2str(DATA_Tremor(i).site) '.mat']);
       
if isempty(DATA_Tremor(i).interval) %Skip data set if interval empty (eg Apo)
    error(['no interval given for site S' mat2str(DATA_Tremor(i).site) './n']);
end

% disregard 1st DP, unlabelled and stim channel
data = data(2:end,1:end); %[2:end] to crop 1st data point
t_original = (1:length(data)-1)*dt;
out = find(strcmp(DATA_Tremor(i).channel,'-') | strcmp(DATA_Tremor(i).channel,'stim'));

% disregard EMG if desired
if args.noEMG == 1 
    out = unique([out;find(strncmp(cellfun(@fliplr,DATA_Tremor(i).channel, ...
                        'UniformOutput',0),'er',2))]); %find right hand EMGs
    out = unique([out;find(strncmp(cellfun(@fliplr,DATA_Tremor(i).channel, ...
                        'UniformOutput',0),'il',2))]); %find left hand EMGs
end
in = setdiff(1:length(DATA_Tremor(i).channel),out);

% load specific channel if desired
if ~isempty(channel)
     if iscell(channel)
         [~,in,~] = intersect(DATA_Tremor(i).channel,channel); % find given channel names (sorts alphabetically)
     else
        in = channel; % only use given channel indices
     end
end


%% Write data into output variables
cc = 1;
for nch = in
            
    T = DATA_Tremor(i).interval{nch};
    if isempty(T); 
        continue
    end
    CHlab{cc}   = DATA_Tremor(i).channel{nch};
    Tint        = t_original >= T(1) & t_original <= T(2);
    dat{cc}     = data(Tint,nch);
    t{cc}       = t_original(Tint);

    %% Handle artefacts
    art{cc} = DATA_Tremor(i).artefact{nch};
    %Convert art in units of time to art in units of data points
    if ~isempty(art{cc})
        [i1,~] = size (art{cc});
        for k = 1:i1
            ti = [find(t{cc} >= art{cc}(k,1),1,'first') ...
                find(t{cc} <= art{cc}(k,2),1,'last')];
            art{cc}(k,:) = ti;
        end
    end

    cc = cc+1;

end

% Determine common time vector for LFP and EMG if desired
if args.commontime == 1
    for i = 1:length(t)
        ti(i,:) = minmax(t{i});
    end
    ti = [max(ti(:,1)) min(ti(:,2))];
    T = t{i}( t{i} >= ti(1) & t{i} <= ti(2) );

    k = 1;
    for i = 1:length(dat)
        ni = t{i} >= ti(1) & t{i} <= ti(2);
        tmp(:,k) = dat{i}(ni);
        art{i} = art{i}-ni(1)+1; % correct artefact positions
        k = k+1;
    end
    dat = tmp;
    t = T;
end

loc = DATA_Tremor(site).location;
pos = DATA_Tremor(site).position;
bad = DATA_Tremor(site).isbad(in);
end

