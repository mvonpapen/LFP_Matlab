%% Load data and check for changes in the power of frequency bands
clear all
load '/afs/geo.uni-koeln.de/usr/neuro/DATA/DATA_last.mat'
[t, dat, ndat, rigor, CHlab, art]=...
    load_dat_from_DATA ( DATA, 'value1', 'Ruhe', 'activity', 1);

%% Begin and end of experiments (must be known)
i0   = [1  6 11 17 25 31 36 41 50 57];
iend = [5 10 16 24 30 34 40 49 56 58];
nex=length(i0);

%% Sort data into begin and end of experiment
for i=1:nex
    % Start of experiment
    t0{i}=t{i0(i)};
    dat0{i}=dat{i0(i)};
    CHlab0{i}=CHlab{i0(i)};
    art0{i}=art{i0(i)};
    % End of experiment
    tend{i}=t{iend(i)};
    datend{i}=dat{iend(i)};
    CHlabend{i}=CHlab{iend(i)};
    artend{i}=art{iend(i)};
end
ndat0=ndat(i0);
rigor0=rigor(i0);
ndatend=ndat(iend);
rigorend=rigor(iend);
clear t dat ndat rigor CHlab art


%% Set variables
nf=50;
% WLFP0=cell(nex,1);
% WLFPend=cell(nex,1);
% coiLFP0=cell(nex,1);
% coiLFPend=cell(nex,1);
PLFP0=NaN(nf,4,nex);
PLFPend=NaN(nf,4,nex);
nLFP=zeros(nex,1);
site=zeros(nex,2);
dt=1/2500;
FF=[1 4; 4 7; 7 13; 13 30; 30 70; 60 80; 250 350];
f=logspace(0,log10(350),nf);
LFPch={'C', 'L', 'A', 'M', 'P'};
% LFPch={'EDCre', 'EDCli', 'FDIre', 'FDIli', 'FDLre', 'FDLli'};

%% Calculate waveklet transform
for i=1:nex
    
    site(i,1)=DATA(ndat0(i)).site;
    site(i,2)=DATA(ndatend(i)).site;

%     %% LFP data
%     iLFP=find(ismember(CHlab0{i},LFPch));
%     nLFP(i)=length(iLFP);
%     k=1;
%     for j=iLFP
%         [tmp, WLFP0{i}{k}, coiLFP0{i}{k}, PLFP0(:,k,i)]...
%                 = procdata (dat0{i}{j}, 'freq', f, 'filter', [45 55; 90 110], ...
%                             'art', {art0{i}{j}});
%         [tmp, WLFPend{i}{k}, coiLFPend{i}{k}, PLFPend(:,k,i)]...
%                 = procdata (datend{i}{j}, 'freq', f, 'filter', [45 55; 90 110], ...
%                             'art', {artend{i}{j}});
%         % Normalize PSD to 50/100Hz energy
%         PLFP0(:,k,i)=PLFP0(:,j,i)/geomean(PLFP0([34 40],j,i));
%         PLFPend(:,k,i)=PLFPend(:,j,i)/geomean(PLFPend([34 40],j,i));
%         k=k+1;
%     end
%     clear tmp
    
    %% Bipolar LFP data
    iLFP=find(ismember(CHlab0{i},LFPch));
    nLFP(i)=length(iLFP);
    if nLFP(i)<2
        continue
    end
    
    % Determine common times
    k=1;    
    for j=iLFP
        ti0(k,:)=minmax(t0{i}{j});
        tiend(k,:)=minmax(tend{i}{j});
        k=k+1;
    end
    ti0=[max(ti0(:,1)) min(ti0(:,2))];
    tiend=[max(tiend(:,1)) min(tiend(:,2))];
    
    % Caclculate dLFP
    k=1;
    for j=iLFP
        ni =  t0{i}{j}>=ti0(1) & t0{i}{j}<=ti0(2);
        dLFP0(:,k)=dat0{i}{j}(ni);
        if k>1
            dLFP0(:,k)=dLFP0(:,k)-dLFP0(:,1);
        end
        art0{i}{j}=art0{i}{j}-ni(1)+1;
        ni =  tend{i}{j}>=tiend(1) & tend{i}{j}<=tiend(2);
        dLFPend(:,k)=datend{i}{j}(ni);
        if k>1
            dLFPend(:,k)=dLFPend(:,k)-dLFPend(:,1);
        end
        artend{i}{j}=artend{i}{j}-ni(1)+1;
        k=k+1;
    end
    clear j k
    dLFP0=dLFP0(:,2:end);
    dLFPend=dLFPend(:,2:end);
    
    % Wavelet transform
    [tmp, WLFP0{i}, coiLFP0{i}, PLFP0(:,1:nLFP(i)-1,i)]...
            = procdata (dLFP0, 'freq', f, 'filter', [], 'art', art0{i});
    [tmp, WLFPend{i}, coiLFPend{i}, PLFPend(:,1:nLFP(i)-1,i)]...
            = procdata (dLFPend, 'freq', f, 'filter', [], 'art', artend{i});
    clear tmp dLFP0 dLFPend ti0 tiend
    
end


%% Plot ratio of power ON/OFF in specified freq. band for all LFPs
% P0=NaN(sum(nLFP),length(FF));
% Pend=NaN(sum(nLFP),length(FF));
for i=1:length(FF)
    tmp1=nanmean(PLFP0(f>=FF(i,1) & f<=FF(i,2),:,rigorend>0));
    tmp2=nanmean(PLFPend(f>=FF(i,1) & f<=FF(i,2),:,rigorend>0));
    P0(:,i)=tmp1(~isnan(tmp1) & ~isnan(tmp2));
    Pend(:,i)=tmp2(~isnan(tmp1) & ~isnan(tmp2));
end
figure, boxplot(Pend./P0,'labels',...
    {'1-4Hz', '4-7Hz', '7-13Hz', '13-30Hz', '30-70Hz', '60-80Hz', '250-350Hz'})
ylim([-0.5 3])
clear i tmp1 tmp2 