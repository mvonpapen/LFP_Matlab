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
zero    = false; % define volume-conduction by 0°--phthres only (not 180°-phthres--180°)


LFPch = {'C', 'L', 'A', 'M', 'P'};
for i=1:Nex
    iLFP = find(ismember(Channel{i}, LFPch));

    LFP{i}=dat{i}(:,iLFP);

    for j=1:length(iLFP)
        ArtLFP{i}{j}=art{i}{iLFP(j)};
    end
    Nch(i) = length(iLFP);
    T(i) = max(t{i}) - min(t{i});

end

clear i dat art t LFPch EMGch


%% Calculate spectral powers: total, coherent, incoherent and volume
%% conducted
nch = max(Nch);

P_tot    = NaN(Nf,nch,Nex);
P_vc     = NaN(Nf,nch,Nex);
P_nvc    = NaN(Nf,nch,Nex);
P_inc    = NaN(Nf,nch,Nex);
P_coh    = NaN(Nf,nch,Nex);
P_bip    = NaN(Nf, 6 ,Nex);

STNelec  = NaN(6,Nex);
Nbc      = NaN(Nex,1);
BiPolCh  = cell(Nex,1);

scale   = (w0+sqrt(2+w0^2))/4/pi ./ f;

for i=1:Nex
    fprintf('Processing data set %u/%u.\n', i, Nex);

    % Determine bipolar combinations w/o pure nonSTN channels
    combo = nchoosek(1:Nch(i),2);
    STNelec(1:size(combo,1),i) = sum(LFPact{i}(combo)>0,2);
    Nbc(i)= size(combo,1);
    BiPolCh{i} = Channel{i}(combo);


    % LFP monopolar
    [~,W,coi]       = procdata(LFP{i}, 'freq', f, 'w0', w0, ...
        'filter', [], 'art', ArtLFP{i}); %'filter', [44 56; 90 110], 
    [Cs, Wxy, W]    = wave_cohere ( W, scale, nsig, ds );
    coi             = coi(1:ds:end,:)/nsig;
    [Pcoh, Pinc, Pvc] = pcc ( f, W, Cs, coi, sig, Wxy, usenan, phthres, zero );
    P_coh(:,1:Nch(i),i) = squeeze(nanmean(Pcoh,3));
    P_inc(:,1:Nch(i),i) = squeeze(nanmean(Pinc,3));
    P_vc(:,1:Nch(i),i)  = squeeze(nanmean(Pvc,3));

    % LFP bipolar
    x = LFP{i}(:,combo(:,1))-LFP{i}(:,combo(:,2));
    a = cell(Nbc(i),1);
    for j=1:Nbc(i);
        a{j} = [ArtLFP{i}{combo(j,1)}; ArtLFP{i}{combo(j,2)}];
    end
    [~,~,~,P_bip(:,1:Nbc(i),i)] = procdata(x, 'freq', f, 'w0', w0, ...
        'filter', [], 'art', a);
end

P_tot = P_inc + P_coh + P_vc;

clear x a W Wxy Cs coi combo



%% Save Variables
save(['LFP_pow_' tag '_' patient '.mat'], ...
    'patient', 'Nex', 'Nbc', 'BiPolCh', 'sig', 'phthres', 'ds', ...
    'P_nvc', 'P_bip', 'P_tot', 'P_coh', 'P_vc', ...
    'f', 'P_inc', 'LFPact', 'w0', 'scale', 'STNelec', ...
    'nsig', 'usenan', 'ds', 'rigor', 'zero', 'T');