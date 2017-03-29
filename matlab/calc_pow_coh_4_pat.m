%% Test mono- and bipolar PSD changes and coherence from OFF to ON

function calc_pow_coh_4_pat(patient)

%% Load LFP data
    
[t, dat, Channel, art, rigor, ~, ~, LFPact] = ...
    load_dat_from_DATA ( 'experiment', 'Ruhe', 'activity', 0, 'commontime', 1, ...
    'noEMG', 1, 'Patient', patient);

%% Parameters
f       = [1:100];
Nf      = length(f);
phthres = 15;
w0      = 12;
nsig    = 6;
sig     = sig_coh_thresh(w0, nsig);
tag     = 'v1';
usenan  = 0; % 1 = take mean only over times where signal is present
ds      = 5;
Nex     = length(dat);
zero    = true; % define volume-conduction by 0°--phthres only (not 180°-phthres--180°)


LFPch = {'C', 'L', 'A', 'M', 'P'};
% EMGch = {'EDCre', 'FDLre', 'FDIre', 'EDCli', 'FDLli', 'FDIli'};
for i=1:Nex
    iLFP = find(ismember(Channel{i}, LFPch));
%     iEMG = find(ismember(Channel{i}, EMGch));

    LFP{i}=dat{i}(:,iLFP);
%     EMG{i}=dat{i}(:,iEMG);

    for j=1:length(iLFP)
        ArtLFP{i}{j}=art{i}{iLFP(j)};
    end
%     for j=1:length(iEMG)
%         ArtEMG{i}{j}=art{i}{iEMG(j)};
%     end

    Nch(i) = length(iLFP);
    T(i) = max(t{i}) - min(t{i});

end

clear i dat art t LFPch EMGch


%% Calculate spectral powers: total, coherent, incoherent and volume
%% conducted
nch = max(Nch);

% P_EMG    = NaN(Nf, 3 ,Nex);
P_tot    = NaN(Nf,nch,Nex);
P_vc     = NaN(Nf,nch,Nex);
P_nvc    = NaN(Nf,nch,Nex);
P_inc    = NaN(Nf,nch,Nex);
P_coh    = NaN(Nf,nch,Nex);
P_bip    = NaN(Nf, 6 ,Nex);
cLL_tot  = NaN(Nf,nch,nch,Nex);
cLL_nvc  = NaN(Nf,nch,nch,Nex);
cLL_bip  = NaN(Nf, 6 , 6 ,Nex);
% cLE_tot  = NaN(Nf,nch,3,Nex);
% cLE_bip  = NaN(Nf, 6 ,3,Nex);

STNelec      = NaN(6,Nex);
Nbc          = NaN(Nex,1);
% Nemg         = NaN(Nex,1);
BiPolCh      = cell(Nex,1);

scale   = (w0+sqrt(2+w0^2))/4/pi ./ f;

for i=1:Nex
    fprintf('Processing data set %u/%u.\n', i, Nex);
%     Nemg(i)=size(EMG{i},2);

    % Determine bipolar combinations w/o pure nonSTN channels
    combo = nchoosek(1:Nch(i),2);
    STNelec(1:size(combo,1),i) = sum(LFPact{i}(combo)>0,2);
%     combo = combo( STNelec(:,i)>0, : ); %no bipolar PSD of nonSTN electrodes
    Nbc(i)= size(combo,1);
    BiPolCh{i} = Channel{i}(combo);


    % LFP monopolar
    [~,W1,coi1,P_tot(:,1:Nch(i),i)] = procdata(LFP{i}, 'freq', f, 'w0', w0, ...
        'filter', [], 'art', ArtLFP{i}); %'filter', [44 56; 90 110], 
    [Cs, Wxys, W1s] = wave_cohere ( W1, scale, nsig, ds );
    coi1s = coi1(1:ds:end,:)/nsig;
    [cLL_tot(:,1:Nch(i),1:Nch(i),i), cLL_nvc(:,1:Nch(i),1:Nch(i),i)] ...
            = acoh1d_uptri( f, Cs, coi1s, Wxys, phthres, sig );
    [Pcoh, Pinc, Pvc, Pnvc] = psd_acoh ( f, W1s, Cs, coi1s, sig, Wxys, usenan, ...
        phthres, zero );
    P_coh(:,1:Nch(i),i) = squeeze(nanmean(Pcoh,3));
    P_inc(:,1:Nch(i),i) = squeeze(nanmean(Pinc,3));
    P_vc(:,1:Nch(i),i)  = squeeze(nanmean(Pvc,3));
    P_nvc(:,1:Nch(i),i) = squeeze(nanmean(Pnvc,3));

    % LFP bipolar
    x = LFP{i}(:,combo(:,1))-LFP{i}(:,combo(:,2));
    a = cell(Nbc(i),1);
    for j=1:Nbc(i);
        a{j} = [ArtLFP{i}{combo(j,1)}; ArtLFP{i}{combo(j,2)}];
    end
    [~,W2,coi2,P_bip(:,1:Nbc(i),i)] = procdata(x, 'freq', f, 'w0', w0, ...
        'filter', [], 'art', a);
    [Cs, Wxys] = wave_cohere ( W2, scale, nsig, ds );
    coi2s = coi2(1:ds:end,:)/nsig;
    cLL_bip(:,1:Nbc(i),1:Nbc(i),i) ...
        = acoh1d_uptri( f, Cs, coi2s, Wxys, phthres, sig );      

%     % EMG
%     if Nemg(i)>0
%         [~,WE,coiE,P_EMG(:,1:Nemg(i),i)] = procdata(EMG{i}, 'freq', f, 'w0', w0, ...
%             'filter', [0 60; 90 110], 'art', ArtEMG{i}, 'rect', 1);
%         Cs = wave_coh3 ( W1, WE, scale, nsig, ds );
%         coiEs = coiE(1:ds:end,:)/nsig;
%         cLE_tot(:,1:Nch(i),1:Nemg(i),i) = coh1d( f, Cs, coi1s, coiEs );
%         Cs = wave_coh3 ( W2, WE, scale, nsig, ds );
%         cLE_bip(:,1:Nbc(i),1:Nemg(i),i) = coh1d( f, Cs, coi2s, coiEs );
%     end

end
    

clear x a W1 W2 WE Wxy C coi1 coi2 coiE combo Cs W1s W2s Wxys coi1s coi2s coiEs



%% Save Variables
save(['pow_coh_' tag '_' patient '.mat'], ...
    'cLL_bip', 'cLL_nvc', 'patient', 'Nex', 'Nbc', 'BiPolCh', ...
    'cLL_tot', 'P_nvc', 'P_bip', 'P_tot', 'P_coh', 'P_vc', ...
    'f', 'P_inc', 'LFPact', 'w0', 'scale', 'STNelec', ...
    'nsig', 'usenan', 'ds', 'rigor', 'zero', 'T');
%     'cLE_bip', 'cLE_tot', 'P_EMG', 'Nemg'