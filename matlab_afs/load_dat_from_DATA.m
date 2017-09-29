

function [t, dat, CHlab, art, rigor, ndat, patients, LFPact] = load_dat_from_DATA ( varargin )
%% load_dat_from_DATA 
% Function to load all data from control file for given values
%
% OPTIONAL INPUT:
%
%   experiment:     e.g. Taps, Ruhe, ... (default Ruhe)
%
%   patient:        patient name e.g. P02_L
%
%   rigor:          find sites with exact rigor value: 0-100
%
%   apomorphine:    yes=1, no=0
%
%   channel:        give channel names or index
%
%   commontime:     find common interval for all channels in one site
%
%   noEMG:          disregard all EMGs
%
% OUTPUT:
%
%   t:              time cell array (2D) of all sites and each channel, if
%                   commontime=1 for site n: t{n}=[t_min(n) t_max(n)] (1D
%                   cell array)
%
%   dat:            data cell array (2D) of all sites and each channel. if
%                   commontime=1 1D cell array and matrix of site per cell
%                   element
%
%   rigor:          rigor value for each site
%
%   CHlab:          channel list of site (1D cell array)
%
%   art:            artifact list of each site and channel (2D cell array)
%
%   ndat:           indices (in DATA_last.mat file) of sites used 
%

args=struct('experiment',[],...
            'patient',[],...
            'rigor',[],...
            'activity',0,...
            'apomorphine', 1,...
            'channel',[],...
            'commontime', 0, ...
            'noEMG',0); %e.g. channel='FDI'
args=parseArgs(varargin,args);
% channel=args.channel;
load /home/mvonpapen/Neuro/DATA/DATA_last.mat


%% Set output in case of crash
t{1}{1}=[];
dat{1}{1}=[];
rigor=[];
CHlab{1}{1}=[];
art{1}{1}=[];
T2{1}=[];
LFPact{1}=[];

if isempty(args.experiment)
    ndat=structfind(DATA,'experiment','Ruhe');
else
    ndat=structfind(DATA,'experiment',args.experiment);
end
if ~isempty(args.patient)
    ndat=intersect(ndat,structfind(DATA,'patient',args.patient));
end
if ~isempty(args.rigor)
    ndatrigor=[];
    for i=1:length(args.rigor)        
        ndatrigor=unique([ndatrigor,structfind(DATA,'rigor', args.rigor(i))]);
    end
    ndat=intersect(ndat,ndatrigor);
end
ndat = intersect(ndat,structfind(DATA,'Apo',args.apomorphine));

%% patientnames
patients=cell(size(ndat));
for i=1:length(ndat)
    patients{i}=DATA(ndat(i)).patient;
end
%% Fixed parameters
j=0;
dt=1/2500;



for i=ndat;

    %% Write out data
    j=j+1;
    load(['/home/mvonpapen/Neuro/DATA/' DATA(i).patient ...
       '_S' mat2str(DATA(i).site) '.mat'])
    if isempty(DATA(i).interval) %Skip data set if interval empty (eg Apo)
        continue
    end
    cc=1; %channel number variable
    data=data(2:end,1:end-1); %[2:end] to crop 1st data point and disregard STIM
    t_original=(1:length(data)-1)*dt;
    out = find(strcmp(DATA(i).channel,'-'));
    %% disregard EMGs
    if args.noEMG==1 
        out=unique([out;find(strncmp(cellfun(@fliplr,DATA(i).channel, ...
                            'UniformOutput',0),'er',2))]); %find right hand EMGs
        out=unique([out;find(strncmp(cellfun(@fliplr,DATA(i).channel, ...
                            'UniformOutput',0),'il',2))]); %find left hand EMGs
    end
    in = setdiff(1:length(DATA(i).interval),out);

    for nch=in

        %% Use only LFP data that shows activity
        if ismember(DATA(i).channel{nch},{'C','L','A','M','P'})
            LFPact{j}(cc)=DATA(i).LFPactivity(nch);
            if DATA(i).LFPactivity(nch)<args.activity
                continue
            end
        end
        T=DATA(i).interval{nch};
        if isempty(T); 
            continue
        end
        CHlab{j}{cc}=DATA(i).channel{nch};
        Tint=t_original>=T(1) & t_original<=T(2);
        dat{j}{cc}=data(Tint,nch);
        t{j}{cc}=t_original(Tint);

        %% Handle artefacts
        art{j}{cc}=DATA(i).artefact{nch};
        %Convert art in units of time to art in units of data points
        if ~isempty(art{j}{cc})
            [i1,i2] = size (art{j}{cc});
            for k=1:i1
                ti=[find(t{j}{cc}>=art{j}{cc}(k,1),1,'first') ...
                    find(t{j}{cc}<=art{j}{cc}(k,2),1,'last')];
                art{j}{cc}(k,:)=ti;
            end
        end
        
        cc=cc+1;
    end
    if args.commontime==1
        %% Determine common time vector for LFP and EMG
            ti2=zeros(length(t{j}),2);
        for i2=1:length(t{j})
            ti2(i2,:)=[min(t{j}{i2}) max(t{j}{i2})];
        end
        ti2=[max(ti2(:,1)) min(ti2(:,2))];
        T2{j} = t{j}{1}( t{j}{1}>=ti2(1) & t{j}{1}<=ti2(2) );

        k=1;
        for i2=1:length(dat{j})
            ni =  t{j}{i2}>=ti2(1) & t{j}{i2}<=ti2(2);
            tmp{j}(:,k)=dat{j}{i2}(ni);
            art{j}{i2}=art{j}{i2}-ni(1)+1; % correct artefact positions
            k=k+1;
        end
        dat{j}=tmp{j};
        t{j}=T2{j};
        clear ti2
    end

    rigor(j)=DATA(i).rigor;

end
end

