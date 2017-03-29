%% This procedure analyzes the spectral changes from OFF to ON separately
%% for the volume conducted and non-volume conducted fluctuations.
%%
%% Here, volume conduction is determined by the phase between two signals
%% being smaller than a certain threshold. Phases are only calculated
%% between electrodes within the STN. If more than two electrodes lie
%% within the STN, then the phases of all electrode combinations must be
%% smaller than the threshold.

clear all

%% Load LFP data
[t, dat, CH, art, rigor, ndat, pat, LFPactivity] = ...
    load_dat_from_DATA ( 'experiment', 'Ruhe', 'activity', 2, 'commontime', 1, 'noEMG', 1);

load /afs/.geo.uni-koeln.de/usr/neuro/DATA/DATA_last.mat

%% Parameters
Nf = 30;
f = logspace(0, log10(90), Nf);
sig = 0.495; %significance threshold for 2D-coherence
phase_thresh = 10;

%% Determine sites and analysis
% Sites with rest ON OFF and rigor>=30
i0   = [2  6 11 18 26 32 36 41 47 65];
iend = [5 10 17 25 31 68 40 46 55 67];
% % Sites with hold ON OFF and nLFP>1
% i0   = [1 3 5  9 11 21];
% iend = [2 4 6 10 23 22];
% % Sites with fist ON OFF and rigor>=30
% i0   = [1 3 6 8 10 11 14 18];
% iend = [2 4 7 9 20 12 15 19];
% i0=6; iend=10;

Nex = length(i0);



for i=1:Nex
    DataOFF{i}=dat{i0(i)};
    ArtOFF{i}=art{i0(i)};
    
    DataON{i}=dat{iend(i)};
    ArtON{i}=art{iend(i)};
    
    Channel{i}=CH{i0(i)};
    Patient{i}=DATA(ndat(i0(i))).patient;
    Patient2{i}=DATA(ndat(iend(i))).patient;
    Nch(i)=size(dat{iend(i)},2);
    LFPact{i}=LFPactivity{iend(i)};
    
    if ~strcmp(Patient{i},Patient2{i})
        error('ERROR: i0 and iend refer to different patients for i=%i!', i)
    end
end

clear i i0 iend DATA dat art CH LFPactivity t ndat pat Patient2

   

%% Calculate wavelet transformation
PVC_OFF         = NaN(Nf,max(Nch),Nex);
PVC_ON          = NaN(Nf,max(Nch),Nex);
PTOT_OFF        = NaN(Nf,max(Nch),Nex);
PTOT_ON         = NaN(Nf,max(Nch),Nex);

for i=1:Nex
        
    if Nch(i)>1
    
        [x,Woff,coi_off] = procdata(DataOFF{i}, 'freq', f, 'filter', [], 'art', ArtOFF{i});
        [x,Won, coi_on]  = procdata(DataON{i},  'freq', f, 'filter', [], 'art', ArtON{i} );
        nt_off = length(DataOFF{i});
        nt_on  = length(DataON{i} );
        
        % Calculate cross-spectra of all channels w.r.t. each other
        Wxy_off = NaN(Nf, nt_off, Nch(i), Nch(i));
        Wxy_on  = NaN(Nf, nt_on , Nch(i), Nch(i));
        for j1=1:Nch(i);
            for j2=setdiff(1:Nch(i),j1);
                Wxy_off(:,:,j1,j2)  = squeeze(Woff(:,:,j1).*conj(Woff(:,:,j2)));
                Wxy_on(:,:,j1,j2)   = squeeze(Won(:,:,j1) .*conj(Won(:,:,j2)) );
            end
        end
        
        % Calculate volume conducted PSD
        [PVC_OFF(:,1:Nch(i),i), PTOT_OFF(:,1:Nch(i),i)] = ...
            psd_vc ( Woff, Wxy_off, phase_thresh, coi_off, f);
        [PVC_ON(:,1:Nch(i),i) , PTOT_ON(:,1:Nch(i),i) ] = ...
            psd_vc ( Won,  Wxy_on,  phase_thresh, coi_on,  f); 
    
    end
    
end
clear x DataOFF DataON ArtOFF ArtON