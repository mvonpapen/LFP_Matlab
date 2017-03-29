%% Test mono- and bipolar PSD changes and coherence from OFF to ON

function calc_pow_coh_4_bip(task)

%% Load LFP data
load(['data_' task '_v3']);

%% Parameters
f       = [1:40 45:5:90];
Nf      = length(f);
sig     = 0.43;
w0      = 12;
nsig    = 6;
ds      = 10;
tag     = 'bip-cohinc';
Nex     = length(LFP_OFF);


%% Calculate spectral powers: total, coherent, incoherent and volume
%% conducted

P_EMG_off    = NaN(Nf, 3, Nex);
P_inc_off    = NaN(Nf, 6, Nex);
P_coh_off    = NaN(Nf, 6, Nex);
P_bip_off    = NaN(Nf, 6, Nex);
cLL_bip_off  = NaN(Nf, 6, 6, Nex);
cLE_bip_off  = NaN(Nf, 6, 3, Nex);

P_EMG_on     = NaN(Nf, 3, Nex);
P_coh_on     = NaN(Nf, 6, Nex);
P_inc_on     = NaN(Nf, 6, Nex);
P_bip_on     = NaN(Nf, 6, Nex);
cLE_tot_on   = NaN(Nf, 6, 3, Nex);
cLE_bip_on   = NaN(Nf, 6, 3, Nex);

STNelec      = NaN(6, Nex);
Nbc          = NaN(Nex, 1);
Nemg         = NaN(Nex, 1);
BiPolCh      = cell(Nex, 1);

scale   = (w0+sqrt(2+w0^2))/4/pi ./ f;

for i=1:Nex
    fprintf('Processing data set %u/%u.\n', i, Nex);
    Nemg(i)=size(EMG_OFF{i},2);

    % Determine bipolar combinations % w/o pure nonSTN channels
    combo = nchoosek(1:Nch(i),2);
%     STNelec(1:size(combo,1),i) = sum(LFPact{i}(combo)>0,2);
%     combo = combo( STNelec(:,i)>0, : ); %no bipolar PSD of nonSTN electrodes
    Nbc(i)= size(combo,1);
    BiPolCh{i} = Channel{i}(combo);



    %% OFF

    % LFP bipolar
    x = LFP_OFF{i}(:,combo(:,1))-LFP_OFF{i}(:,combo(:,2));
    a = cell(Nbc(i),1);
    for j=1:Nbc(i);
        a{j} = [ArtLFP_OFF{i}{combo(j,1)}; ArtLFP_OFF{i}{combo(j,2)}];
    end
    [~,W2,coi2,P_bip_off(:,1:Nbc(i),i)] = procdata(x, 'freq', f, 'w0', w0, ...
        'filter', [], 'art', a);
    [Cs, Wxys, W2s] = wave_cohere ( W2, scale, nsig, ds );
    coi2s = coi2(1:ds:end,:)/nsig;
    cLL_bip_off(:,1:Nbc(i),1:Nbc(i),i) = acoh1d_uptri( f, Cs, coi2s );
    [Pcoh, Pinc] = psd_acoh ( f, W2s, Cs, coi2s, sig, Wxys, 0, 0 );
    P_coh_off(:,1:Nbc(i),i) = squeeze(nanmean(Pcoh,3));
    P_inc_off(:,1:Nbc(i),i) = squeeze(nanmean(Pinc,3));

    % EMG
    if Nemg(i)>0
        [~,WE,coiE,P_EMG_off(:,1:Nemg(i),i)] = procdata(EMG_OFF{i}, 'freq', f, 'w0', w0, ...
            'filter', [0 60; 90 110], 'art', ArtEMG_OFF{i}, 'rect', 1);
        coiEs = coiE(1:ds:end,:)/nsig;
        Cs = wave_coh3 ( W2, WE, scale, nsig, ds );
        cLE_bip_off(:,1:Nbc(i),1:Nemg(i),i) = coh1d( f, Cs, coi2s, coiEs );
    end


    %% ON

    % LFP bipolar
    x = LFP_ON{i}(:,combo(:,1))-LFP_ON{i}(:,combo(:,2));
    a = cell(Nbc(i),1);
    for j=1:Nbc(i);
        a{j} = [ArtLFP_ON{i}{combo(j,1)}; ArtLFP_ON{i}{combo(j,2)}];
    end
    [~,W2,coi2,P_bip_on(:,1:Nbc(i),i)] = procdata(x, 'freq', f, ...
        'w0', w0, 'filter', [], 'art', a);
    coi2s = coi2(1:ds:end,:)/nsig;
    [Cs, Wxys, W2s] = wave_cohere ( W2, scale, nsig, ds );
    cLL_bip_on(:,1:Nbc(i),1:Nbc(i),i) = acoh1d_uptri( f, Cs, coi2s );
    [Pcoh, Pinc] = psd_acoh ( f, W2s, Cs, coi2s, sig, Wxys, 0, 0 );
    P_coh_on(:,1:Nbc(i),i) = squeeze(nanmean(Pcoh,3));
    P_inc_on(:,1:Nbc(i),i) = squeeze(nanmean(Pinc,3));

    % EMG
    if Nemg(i)>0
        [~,WE,coiE,P_EMG_on(:,1:Nemg(i),i)] = procdata(EMG_ON{i}, 'freq', f, 'w0', w0, ...
            'filter', [0 60; 90 110], 'art', ArtEMG_ON{i}, 'rect', 1);
        coiEs = coiE(1:ds:end,:)/nsig;
        Cs = wave_coh3 ( W2, WE, scale, nsig, ds );
        cLE_bip_on(:,1:Nbc(i),1:Nemg(i),i) = coh1d( f, Cs, coi2s, coiEs );
    end

end

clear x a W1 W2 WE Wxy C coi1 coi2 coiE combo Cs W1s W2s Wxys coi1s coi2s coiEs



%% Save Variables
save(['pow_coh_' tag '_' task '.mat'], ...
    'P_bip_off', 'P_bip_on', 'f', 'Patient', 'LFPact', 'w0', 'scale', ...
    'P_EMG_off', 'P_EMG_on', 'P_coh_off', 'P_coh_on', 'P_inc_off', 'P_inc_on', ...
    'cLL_bip_off', 'cLL_bip_on', 'cLE_bip_off', 'cLE_bip_on', ...
    'task', 'Nex', 'nsig', 'ds', 'rigor')