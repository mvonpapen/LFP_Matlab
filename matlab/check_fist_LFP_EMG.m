%% Check Phase Fist ON/OFF von LFP-EMG

%% Load data from variable DATA and put it in the pipe!
clear all
load '/afs/geo.uni-koeln.de/usr/neuro/DATA/DATA_last.mat'
[t, dat, ndat, rigor, CHlab, art]=...
    load_dat_from_DATA ( DATA, 'value1', 'Faust', 'activity', 1);

%% Begin and end of experiments
i0   = [1 3 5  8  9 12 14];
iend = [2 4 6 17 10 13 15];

%% Set variables
nf=30;
f=logspace(-1,1.5,nf);
dt=1/2500;
ns=length(dat);
t_smo=3;
sc_smo=0;

%% Define LFP channels
LFPch={'C', 'L', 'A', 'M', 'P'};
EMGch={'EDCre', 'EDCli', 'FDIre', 'FDIli', 'FDLre', 'FDLli'};

%% Wavelet trafo von LFP
for n=1:ns
        
    %% Determine common time vector for LFP and EMG
    for i=1:length(t{n})
        ti(i,:)=minmax(t{n}{i});
    end
    ti=[max(ti(:,1)) min(ti(:,2))];
    T{n} = t{n}{i}( t{n}{i}>=ti(1) & t{n}{i}<=ti(2) );
    
    %% LFP data
    iLFP=find(ismember(CHlab{n},LFPch));
    nLFP(n)=length(iLFP);
    k=1;
    for i=iLFP
        ni =  t{n}{i}>=ti(1) & t{n}{i}<=ti(2);
        tmp(:,k)=dat{n}{i}(ni);
        art{n}{i}=art{n}{i}-ni(1)+1;
        k=k+1;
    end
    [LFP{n}, WLFP{n}, coiLFP{n}, PLFP(:,1:nLFP(n),n)]...
            = procdata (tmp, 'freq', f, 'filter', [45 55; 90 110], ...
                        'art', art{n});
    clear tmp ni i k
    
    %% EMG data
    iEMG=find(ismember(CHlab{n},EMGch));
    nEMG(n)=length(iEMG);
    if nEMG(n)>0
        k=1;
        for i=iEMG;
            ni =  t{n}{i}>=ti(1) & t{n}{i}<=ti(2);
            tmp(:,k)=dat{n}{i}(ni);
            art{n}{i}=art{n}{i}-ni(1)+1;
            k=k+1;
        end
        [EMG{n}, WEMG{n}, coiEMG{n}, PEMG(:,1:nEMG(n),n)]...
                =procdata(tmp,'filter', [0 60; 90 110], 'freq', f, ...
                'art', art{n}, 'rect', 1);
        clear tmp ni i k
    end
    
    %% Coherence
    if nLFP(n)>0 && nEMG(n)>0
        [CLE{n}, WLE{n}] = wave_coh(WLFP{n}, WEMG{n}, f, t_smo, sc_smo);
    end
    
end