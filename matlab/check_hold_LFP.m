%% Check Halte ON/OFF von LFP

%% Load data from variable DATA and put it in the pipe!
clear all
load '/afs/geo.uni-koeln.de/usr/neuro/DATA/DATA_last.mat'
[t, dat, ndat, rigor, CHlab, art]=...
    load_dat_from_DATA ( DATA, 'value1', 'Halte', 'activity', 1);

%% Set variables
nf=50;
f=logspace(0,2.7,nf);
dt=1/2500;
ns=length(dat);

%% Define LFP channels
LFPch={'C', 'L', 'A', 'M', 'P'};

%% Wavelet trafo von LFP
for n=1:ns
        
    iLFP=find(ismember(CHlab{n},LFPch));
    nLFP(n)=length(iLFP);
    site(n)=DATA(ndat(n)).site;
    pat{n}=DATA(ndat(n)).patient;
    k=1;
    for i=iLFP
        [LFP{n,i}, WLFP{n,i}, coiLFP{n,i}, PLFP(:,k,n)]...
            = procdata (dat{n}{i}, 'freq', f, 'filter', [], 'art', art{n}{i});
        k=k+1;
    end
    clear i k
    
end

%% Plot ratio of power ON/OFF in specified freq. band
% Begin and end of experiments
i0   = [1 3 5 7  9 11 12 15 17 19 21];
iend = [2 4 6 8 10 23 13 16 18 20 22];

%% Brain frequency bands
fb{1}=find(f>=1 & f<4);
fb{2}=find(f>=4 & f<7);
fb{3}=find(f>=7 & f<13); 
fb{4}=find(f>=13 & f<30);
fb{5}=find(f>=30 & f<70);
fb{6}=find(f>=60 & f<80);
fb{7}=find(f>=250 & f<350);

for m=1:7;
    k=1;
    for i=1:length(i0)
        for j=1:nLFP(i0(i))
            P0(k,m)=nanmean(PLFP(fb{m},j,i0(i)));
            Pend(k,m)=nanmean(PLFP(fb{m},j,iend(i)));
            k=k+1;
        end
    end
end

figure, boxplot(Pend./P0,'labels',...
    {'1-4Hz', '4-7Hz', '7-13Hz', '13-30Hz', '30-70Hz', '60-80Hz', '250-350Hz'})