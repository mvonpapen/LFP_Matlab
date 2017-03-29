%% Calculate all needed variables (coherence, PSD, phases) to check phase locking OFF/ON
clear all

%% Load LFP data
[t, dat, CH, art, rigor, ndat, pat, LFPactivity] = ...
    load_dat_from_DATA ( 'experiment', 'Faust', 'activity', 0, 'commontime', 1);
for i=1:length(dat)
    nLFP(i)=size(dat{i},2);
end
load /afs/.geo.uni-koeln.de/usr/neuro/DATA/DATA_last.mat

%% Parameters
Nf = 6;
f = logspace(0, log10(3), Nf);

%% Determine sites and analysis
% Sites with fist ON OFF and rigor>=30
i0   = [1 3 6 8 10 11 14 18];
iend = [2 4 7 9 20 12 15 19];
Nex = length(i0);

%% Define LFP and EMG channels
LFPch={'C', 'L', 'P', 'A', 'M'};
EMGch={'EDCre', 'EDCli', 'FDIre', 'FDIli', 'FDLre', 'FDLli'};

for i=1:Nex
    DataOFF{i}=dat{i0(i)};
    ArtOFF{i}=art{i0(i)};
    
    DataON{i}=dat{iend(i)};
    ArtON{i}=art{iend(i)};
    
    Channel{i}=CH{i0(i)};
    Patient{i}=DATA(ndat(iend(i))).patient;
    Nch(i)=size(dat{iend(i)},2);
    LFPact{i}=LFPactivity{iend(i)};
    
    iLFP{i}=find(ismember(Channel{i},LFPch));
    nLFP(i)=length(iLFP{i});
    iEMG{i}=find(ismember(Channel{i},EMGch));
    nEMG(i)=length(iEMG{i});
end

clear i i0 iend DATA

   

%% Calculate wavelet transformation
Poff    = NaN(Nf,max(Nch),Nex);
Pon     = NaN(Nf,max(Nch),Nex);
c1d_off = NaN(Nf,max(nLFP),max(nEMG));
c1d_on  = NaN(Nf,max(nLFP),max(nEMG));

for i=1:Nex
    [x,Woff,coi_off{i},Poff(:,1:Nch(i),i)] = procdata(DataOFF{i}, 'freq', f, ...
        'filter', [], 'art', ArtOFF{i});
    [x,Won, coi_on{i}, Pon(:,1:Nch(i),i)]  = procdata(DataON{i},  'freq', f, ...
        'filter', [], 'art', ArtON{i} );
    
    % Calculate coherent spectra of all channels w.r.t. non-STN
    [C_off,Wxy_off{i}] = wave_coh ( Woff(:,:,iLFP{i}), Woff(:,:,iEMG{i}), f );
    [C_on, Wxy_on{i}]  = wave_coh ( Won(:,:,iLFP{i}),  Won(:,:,iEMG{i}), f  );
    
    for j=1:nLFP(i)
        for k=1:nEMG(i)
            c1d_off(:,j,k) = nanmean( coi2nan( f, C_off(:,:,j,k), ...
                min([coi_off{i}(:,j)'; coi_off{i}(:,k+nLFP(i))']) ),2 );
            c1d_on(:,j,k)  = nanmean( coi2nan( f, C_on(:,:,j,k),  ...
                min([coi_on{i}(:,j)';  coi_on{i}(:,k+nLFP(i))' ]) ),2 );
        end
    end
                  
end