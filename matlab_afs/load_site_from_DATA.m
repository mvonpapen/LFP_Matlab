function [ t, dat, CHlab, art, rigor] = load_site_from_DATA( DATA_num_site, varargin )
%LOAD_SITE_FROM_DATA  Load Data for a given Site number (DATA_num_site) in 
% DATA_last control file.
%   OPTIONAL INPUT:
%           activity:   minimum LFP activity to use (default 0 = all)
%
%           channel:    cell of channel names or indice to use (default
%                       all)
%
%           noEMG:      Disregard all EMGs (default is 0)
%
%           commontime: find common interval of all channels (default is 1)
%
%   OUTPUT:
%           t:          time vector if commontime, otherwise cell array
%
%           dat:        matrix of data if commontime, otherwise cell array
%
%           CHlab:      names of output channels (cell array)
%
%           art:        artifacts for channels (cell array)
%
%           rigor:      rigor value for the site
%
%
% Author: Michael von Papen, Felix Gerick
%
% Date: 13.10.15

i = DATA_num_site;
load('/afs/rrz/geo/project/neuro/DATA/DATA_last.mat');

args = struct('activity',0,...
            'channel',[],...
            'noEMG',0, ...
            'commontime',1);
args = parseArgs(varargin,args);
channel = args.channel;
dt = 1/2500;

load(['/afs/rrz/geo/project/neuro/DATA/' DATA(i).patient ...
           '_S' mat2str(DATA(i).site) '.mat']);
       
if isempty(DATA(i).interval) %Skip data set if interval empty (eg Apo)
    error('no interval given');
end

%% disregard false data
data = data(2:end,1:end-1); %[2:end] to crop 1st data point and disregard STIM
t_original = (1:length(data)-1)*dt;
out = find(strcmp(DATA(i).channel,'-')); %find unlabeled channels
%% disregard EMGs
if args.noEMG == 1 
    out = unique([out;find(strncmp(cellfun(@fliplr,DATA(i).channel, ...
                        'UniformOutput',0),'er',2))]); %find right hand EMGs
    out = unique([out;find(strncmp(cellfun(@fliplr,DATA(i).channel, ...
                        'UniformOutput',0),'il',2))]); %find left hand EMGs
end
in = setdiff(1:length(DATA(i).interval),out);

%% for given channel index or names
if ~isempty(channel)
     if iscell(channel)
         [a,in,a] = intersect(DATA(i).channel,channel); % find given channel names (sorts alphabetically)
         clear a;
     else
        in = channel; % only use given channel indice
     end
end
%%
cc = 1;
for nch = in
            
    %% Use only LFP data that shows activity
    if nch <= length(DATA(i).LFPactivity)
        if DATA(i).LFPactivity(nch)<args.activity
            continue
        end
    end
    T = DATA(i).interval{nch};
    if isempty(T); 
        continue
    end
    CHlab{cc} = DATA(i).channel{nch};
    Tint = t_original >= T(1) & t_original <= T(2);
    dat{cc} = data(Tint,nch);
    t{cc} = t_original(Tint);

    %% Handle artefacts
    art{cc} = DATA(i).artefact{nch};
    %Convert art in units of time to art in units of data points
    if ~isempty(art{cc})
        [i1,i2] = size (art{cc});
        for k = 1:i1
            ti = [find(t{cc} >= art{cc}(k,1),1,'first') ...
                find(t{cc} <= art{cc}(k,2),1,'last')];

            art{cc}(k,:) = ti;
        end
    end

     cc = cc+1;

end

if args.commontime == 1
    %% Determine common time vector for LFP and EMG
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

rigor = DATA(DATA_num_site).rigor;
end

