%% Test mono- and bipolar PSD changes and coherence from OFF to ON


%% Load LFP data

for type = {'Faust', 'Halte', 'Ruhe'};
type = type{1}

[t, dat, CH, art, rigor, ndat, x, LFPactivity] = ...
    load_dat_from_DATA ( 'experiment', type, 'activity', 0, 'commontime', 1, 'noEMG', 0);

load DATA_last.mat

%% Parameters
f       = [1:30 60:5:90 250:25:350];
Nf      = length(f);
sig     = 0.495;
phthres = 10;
t_smo   = 3;
sc_smo  = 0;
uptri   = 1; %only upper triangle of 1D-coherence
usenan  = 0;




%% Determine sites and analysis
switch type
    case 'Ruhe'
        % Sites with rest ON OFF, Nelec>1 and rigor>=30
        i0   = [1  6 11 26 32 41 56 65]; % last entry: no EMG channels
        iend = [5 10 17 31 68 46 62 67];
    case 'Halte'
        % Sites with hold ON OFF and nLFP>1
        i0   = [1 3 5  9 11 17 21];
        iend = [2 4 6 10 23 18 22];
    case 'Faust'
        % Sites with fist ON OFF and rigor>=30
        i0   = [1 3 8 10 16 18];
        iend = [2 4 9 20 17 19];
end

Nex = length(i0);

LFPch = {'C', 'L', 'A', 'M', 'P'};
EMGch = {'EDCre', 'FDLre', 'FDIre', 'EDCli', 'FDLli', 'FDIli'};
for i=1:Nex
    Channel{i}=CH{iend(i)};
    iLFP = find(ismember(Channel{i}, LFPch));
    iEMG = find(ismember(Channel{i}, EMGch));
    
    LFP_OFF{i}=dat{i0(i)}(:,iLFP);
    EMG_OFF{i}=dat{i0(i)}(:,iEMG);
    LFP_ON{i}=dat{iend(i)}(:,iLFP);
    EMG_ON{i}=dat{iend(i)}(:,iEMG);
    
    for j=1:length(iLFP)
        ArtLFP_OFF{i}{j}=art{i0(i)}{iLFP(j)};
        ArtLFP_ON{i}{j}=art{iend(i)}{iLFP(j)};
    end
    for j=1:length(iEMG)
        ArtEMG_OFF{i}{j}=art{i0(i)}{iEMG(j)};
        ArtEMG_ON{i}{j}=art{iend(i)}{iEMG(j)};
    end
        
    Patient{i}=DATA(ndat(i0(i))).patient;
    Patient2{i}=DATA(ndat(iend(i))).patient;
%     Nch(i)=size(dat{iend(i)},2);
    Nch(i) = length(iLFP);
    LFPact{i}=LFPactivity{iend(i)};
    T_off(i) = max(t{i0(i)}) - min(t{i0(i)});
    T_on(i)  = max(t{iend(i)}) - min(t{iend(i)});
    
    
    if ~strcmp(Patient{i},Patient2{i})
        error('ERROR: i0 and iend refer to different patients for i=%i!', i)
    end
end
rigor=rigor(iend);

clear i i0 iend DATA dat art CH LFPactivity t ndat pat Patient2 x


for w0 = [6, 12, 18]
w0

%% Calculate spectral powers: total, coherent, incoherent and volume
%% conducted
nch = max(Nch);

P_EMG_off    = NaN(Nf, 3 ,Nex);
P_tot_off    = NaN(Nf,nch,Nex);
P_vc_off     = NaN(Nf,nch,Nex);
P_inc_off    = NaN(Nf,nch,Nex);
P_coh_off    = NaN(Nf,nch,Nex);
P_bip_off    = NaN(Nf, 6 ,Nex);
cLL_tot_off  = NaN(Nf,nch,nch,Nex);
cLL_nvc_off  = NaN(Nf,nch,nch,Nex);
cLL_bip_off  = NaN(Nf, 6 , 6 ,Nex);
cLE_tot_off  = NaN(Nf,nch,3,Nex);
% cLE_nvc_off  = NaN(Nf,nch,3,Nex);
cLE_bip_off  = NaN(Nf, 6 ,3,Nex);

P_EMG_on     = NaN(Nf, 3 ,Nex);
P_tot_on     = NaN(Nf,nch,Nex);
P_vc_on      = NaN(Nf,nch,Nex);
P_coh_on     = NaN(Nf,nch,Nex);
P_inc_on     = NaN(Nf,nch,Nex);
P_bip_on     = NaN(Nf, 6 ,Nex);
cLL_tot_on   = NaN(Nf,nch,nch,Nex);
cLL_nvc_on   = NaN(Nf,nch,nch,Nex);
cLL_bip_on   = NaN(Nf, 6 , 6 ,Nex);
cLE_tot_on   = NaN(Nf,nch,3,Nex);
% cLE_nvc_on   = NaN(Nf,nch,3,Nex);
cLE_bip_on   = NaN(Nf, 6 ,3,Nex);

STNelec      = NaN(6,Nex);
Nbc          = NaN(Nex,1);
Nemg         = NaN(Nex,1);
BiPolCh      = cell(Nex,1);

scale   = (w0+sqrt(2+w0^2))/4/pi ./ f;

for i=1:Nex
    [i Nex]
    Nemg(i)=size(EMG_OFF{i},2);
    
    % Determine bipolar combinations w/o pure nonSTN channels
    combo = nchoosek(1:Nch(i),2);
    STNelec(1:size(combo,1),i) = sum(LFPact{i}(combo)>0,2);
    combo = combo( STNelec(:,i)>0, : ); %no bipolar PSD of nonSTN electrodes
    Nbc(i)= size(combo,1);
    BiPolCh{i} = Channel{i}(combo);
    
    
    
    %% OFF
    
    % LFP monopolar
    [x,W1,coi1,P_tot_off(:,1:Nch(i),i)] = procdata(LFP_OFF{i}, 'freq', f, 'w0', w0, ...
        'filter', [], 'art', ArtLFP_OFF{i});
    [C, Wxy] = wave_acoh2 ( W1, scale, t_smo, sc_smo );
    [cLL_tot_off(:,1:Nch(i),1:Nch(i),i), cLL_nvc_off(:,1:Nch(i),1:Nch(i),i)] ...
            = acoh1d( f, C, coi1, Wxy, phthres, uptri );
    [Pcoh, Pinc, Pvc] = psd_acoh ( f, W1, C, coi1, sig, Wxy, usenan );
    P_vc_off(:,1:Nch(i),i)  = squeeze(nanmean(Pvc,3));
    P_coh_off(:,1:Nch(i),i) = squeeze(nanmean(Pcoh,3));
    P_inc_off(:,1:Nch(i),i) = squeeze(nanmean(Pinc,3));
%     [cLL_tot_off(:,1:Nch(i),1:Nch(i),i), cLL_nvc_off(:,1:Nch(i),1:Nch(i),i)] ...
%         = acohere1d ( f, W1, coi1, Wxy, C, phthres );
        
    % LFP bipolar
    x = LFP_OFF{i}(:,combo(:,1))-LFP_OFF{i}(:,combo(:,2));
    for j=1:Nbc(i);
        a{j} = [ArtLFP_OFF{i}{combo(j,1)}; ArtLFP_OFF{i}{combo(j,2)}];
    end
    [x,W2,coi2,P_bip_off(:,1:Nbc(i),i)] = procdata(x, 'freq', f, 'w0', w0, 'filter', [], 'art', a);
    [C, Wxy] = wave_acoh2 ( W2, scale, t_smo, sc_smo );
    cLL_bip_off(:,1:Nbc(i),1:Nbc(i),i) = acoh1d( f, C, coi2, Wxy, phthres, uptri );
%     cLL_bip_off(:,1:Nbc(i),1:Nbc(i),i) = acohere1d ( f, W2, coi2 );
    
    % EMG
    if Nemg(i)>0
        [x,WE,coiE,P_EMG_off(:,1:Nemg(i),i)] = procdata(EMG_OFF{i}, 'freq', f, 'w0', w0, ...
            'filter', [0 60; 90 110], 'art', ArtEMG_OFF{i}, 'rect', 1);
        C = wave_coh3 ( W1, WE, scale, t_smo );
        cLE_tot_off(:,1:Nch(i),1:Nemg(i),i) = coh1d( f, C, coi1, coiE );
        C = wave_coh3 ( W2, WE, scale, t_smo );
        cLE_bip_off(:,1:Nbc(i),1:Nemg(i),i) = coh1d( f, C, coi2, coiE );
%         cLE_tot_off(:,1:Nch(i),1:Nemg(i),i) = cohere1d( f, W1, WE, coi1, coiE );
%         cLE_bip_off(:,1:Nbc(i),1:Nemg(i),i) = cohere1d( f, W2, WE, coi2, coiE );
    end
    
    
    %% ON
    
    % LFP monopolar
    [x,W1,coi1,P_tot_on(:,1:Nch(i),i)] = procdata(LFP_ON{i}, 'freq', f, 'w0', w0, ...
        'filter', [], 'art', ArtLFP_ON{i});
    [C, Wxy] = wave_acoh2 ( W1, scale, t_smo, sc_smo );
    [cLL_tot_on(:,1:Nch(i),1:Nch(i),i), cLL_nvc_on(:,1:Nch(i),1:Nch(i),i)] ...
            = acoh1d( f, C, coi1, Wxy, phthres, uptri );
    [Pcoh, Pinc, Pvc] = psd_acoh ( f, W1, C, coi1, sig, Wxy, usenan );
    P_vc_on(:,1:Nch(i),i)  = squeeze(nanmean(Pvc,3));
    P_coh_on(:,1:Nch(i),i) = squeeze(nanmean(Pcoh,3));
    P_inc_on(:,1:Nch(i),i) = squeeze(nanmean(Pinc,3));
%     [cLL_tot_on(:,1:Nch(i),1:Nch(i),i), cLL_nvc_on(:,1:Nch(i),1:Nch(i),i)] ...
%         = acohere1d ( f, W1, coi1, Wxy, C, phthres );
        
    % LFP bipolar
    x = LFP_ON{i}(:,combo(:,1))-LFP_ON{i}(:,combo(:,2));
    for j=1:Nbc(i);
        a{j} = [ArtLFP_ON{i}{combo(j,1)}; ArtLFP_ON{i}{combo(j,2)}];
    end
    [x,W2,coi2,P_bip_on(:,1:Nbc(i),i)] = procdata(x, 'freq', f, 'w0', w0, 'filter', [], 'art', a);
    [C, Wxy] = wave_acoh2 ( W2, scale, t_smo, sc_smo );
    cLL_bip_on(:,1:Nbc(i),1:Nbc(i),i) = acoh1d( f, C, coi2, Wxy, phthres, uptri );
%     cLL_bip_on(:,1:Nbc(i),1:Nbc(i),i) = acohere1d ( f, W2, coi2 );
    
    % EMG
    if Nemg(i)>0
        [x,WE,coiE,P_EMG_on(:,1:Nemg(i),i)] = procdata(EMG_ON{i}, 'freq', f, 'w0', w0, ...
            'filter', [0 60; 90 110], 'art', ArtEMG_ON{i}, 'rect', 1);
        C = wave_coh3 ( W1, WE, scale, t_smo );
        cLE_tot_on(:,1:Nch(i),1:Nemg(i),i) = coh1d( f, C, coi1, coiE );
        C = wave_coh3 ( W2, WE, scale, t_smo );
        cLE_bip_on(:,1:Nbc(i),1:Nemg(i),i) = coh1d( f, C, coi2, coiE );
%         cLE_tot_on(:,1:Nch(i),1:Nemg(i),i) = cohere1d( f, W1, WE, coi1, coiE );
%         cLE_bip_on(:,1:Nbc(i),1:Nemg(i),i) = cohere1d( f, W2, WE, coi2, coiE );
    end
    
end

P_nvc_off = P_tot_off - P_vc_off;
P_nvc_on  = P_tot_on  - P_vc_on ;
P_nvc_off(P_nvc_off< 0) = 0;
P_nvc_on(P_nvc_on  < 0) = 0;

clear x a W1 W2 WE Wxy C coi1 coi2 coiE combo





%% Total power in frequency bands

fband{1} = find(f>=1  & f<=4 );
fband{2} = find(f>=4  & f<=7 );
fband{3} = find(f>=7  & f<=13);
fband{4} = find(f>=13 & f<=30);
fband{5} = find(f>=13 & f<=20);
fband{6} = find(f>=20 & f<=30);
fband{7} = find(f>=60 & f<=90);
fband{8} = find(f>=250 & f<=350);
Nb = length(fband);

PL_tot_off = NaN(nch*Nex, Nb);
PL_tot_on  = NaN(nch*Nex, Nb);
PL_nvc_off = NaN(nch*Nex, Nb);
PL_nvc_on  = NaN(nch*Nex, Nb);
PL_coh_off = NaN(nch*Nex, Nb);
PL_coh_on  = NaN(nch*Nex, Nb);
PL_inc_off = NaN(nch*Nex, Nb);
PL_inc_on  = NaN(nch*Nex, Nb);
PL_bip_off = NaN(6*Nex, Nb);
PL_bip_on  = NaN(6*Nex, Nb);
PL_vc_off  = NaN(nch*Nex, Nb);
PL_vc_on   = NaN(nch*Nex, Nb);
cLL_nvc_off_band = NaN(nch^2*Nex, Nb);
cLL_nvc_on_band  = NaN(nch^2*Nex, Nb);
cLL_bip_off_band = NaN(6^2*Nex, Nb);
cLL_bip_on_band  = NaN(6^2*Nex, Nb);


% Total monopolar
x=P_tot_off(:,:);
y=P_tot_on(:,:);
for i=1:Nb
    PL_tot_off(:,i) = squeeze(nanmean(x(fband{i},:),1));
    PL_tot_on(:,i)  = squeeze(nanmean(y(fband{i},:),1));
end

% No volume conduction
x=P_nvc_off(:,:);
y=P_nvc_on(:,:);
for i=1:Nb
    PL_nvc_off(:,i) = squeeze(nanmean(x(fband{i},:),1));
    PL_nvc_on(:,i)  = squeeze(nanmean(y(fband{i},:),1));
end

% Coherent
x=P_coh_off(:,:);
y=P_coh_on(:,:);
for i=1:Nb
    PL_coh_off(:,i) = squeeze(nanmean(x(fband{i},:),1));
    PL_coh_on(:,i)  = squeeze(nanmean(y(fband{i},:),1));
end

% Volume Conduction
x=P_vc_off(:,:);
y=P_vc_on(:,:);
for i=1:Nb
    PL_vc_off(:,i) = squeeze(nanmean(x(fband{i},:),1));
    PL_vc_on(:,i)  = squeeze(nanmean(y(fband{i},:),1));
end

% Incoherent
x=P_inc_off(:,:);
y=P_inc_on(:,:);
for i=1:Nb
    PL_inc_off(:,i) = squeeze(nanmean(x(fband{i},:),1));
    PL_inc_on(:,i)  = squeeze(nanmean(y(fband{i},:),1));
end

% Bipolar
x=P_bip_off(:,:);
y=P_bip_on(:,:);
for i=1:Nb
    PL_bip_off(:,i) = squeeze(nanmean(x(fband{i},:),1));
    PL_bip_on(:,i)  = squeeze(nanmean(y(fband{i},:),1));
end
clear x y

% Coh nvc
x=cLL_nvc_off(:,:);
y=cLL_nvc_on(:,:);
for i=1:Nb
    cLL_nvc_off_band(:,i) = squeeze(nanmean(x(fband{i},:),1));
    cLL_nvc_on_band(:,i)  = squeeze(nanmean(y(fband{i},:),1));
end
clear x y

% Coh bip
x=cLL_bip_off(:,:);
y=cLL_bip_on(:,:);
for i=1:Nb
    cLL_bip_off_band(:,i) = squeeze(nanmean(x(fband{i},:),1));
    cLL_bip_on_band(:,i)  = squeeze(nanmean(y(fband{i},:),1));
end
clear x y


% Wilcoxon sign-rank test
for i=1:Nb
    sig_nvc(i) = signrank ( log10(PL_nvc_off(:,i)./PL_nvc_on(:,i) ) );
    sig_tot(i) = signrank ( log10(PL_tot_off(:,i)./PL_tot_on(:,i) ) );
    sig_coh(i) = signrank ( log10(PL_coh_off(:,i)./PL_coh_on(:,i) ) );
    sig_inc(i) = signrank ( log10(PL_inc_off(:,i)./PL_inc_on(:,i) ) );
    sig_vc(i)  = signrank ( log10(PL_vc_off(:,i) ./PL_vc_on(:,i)  ) );
    sig_bip(i) = signrank ( log10(PL_bip_off(:,i)./PL_bip_on(:,i) ) );
    sig_coh_nvc(i) = signrank ( cLL_nvc_off_band(:,i)-cLL_nvc_on_band(:,i) );
    sig_coh_bip(i) = signrank ( cLL_bip_off_band(:,i)-cLL_bip_on_band(:,i) );
end





%% Plot results for each patient
% for i=1:Nex
%     
%     % Power and coherence
%     fig = figure('PaperUnits','centimeters', 'visi', 'off', 'PaperSize', [21 29.7]);
%     subplot(4,2,1)
%     semilogy(f,P_nvc_off(:,:,i),'-',f,P_nvc_on(:,:,i),'--')
%     title('Mono - vol.cond.'), xlim([min(f) max(f)])%, ylim([1e-8 1e-6])
%     xlabel('f [Hz]'), ylabel('PSD [V^2/Hz]')
%     subplot(4,2,2)
%     semilogy(f,P_bip_off(:,:,i),'-',f,P_bip_on(:,:,i),'--')
%     title('Bipolar'), xlim([min(f) max(f)])%, ylim([1e-8 1e-6])
%     xlabel('f [Hz]'), ylabel('PSD [V^2/Hz]')
%     subplot(4,2,3)
%     semilogy(f,P_tot_off(:,:,i),'-',f,P_tot_on(:,:,i),'--')
%     title('mono'), xlim([min(f) max(f)])%, ylim([1e-8 1e-6])
%     xlabel('f [Hz]'), ylabel('PSD [V^2/Hz]')
%     subplot(4,2,4)
%     semilogy(f,P_EMG_off(:,:,i),'-',f,P_EMG_on(:,:,i),'--')
%     title('EMG'), xlim([min(f) max(f)])%, ylim([1e-11 1e-8])
%     xlabel('f [Hz]'), ylabel('PSD [V^2/Hz]')
%     for j=1:Nch(i)-1
%         subplot(4,2,5)
%         plot(f,cLL_nvc_off(:,j+1:end,j,i),'-',f,cLL_nvc_on(:,j+1:end,j,i),'--')
%         title('Coh Mono - vol.cond.'), xlim([min(f) max(f)]), ylim([0 1])
%         xlabel('f [Hz]'), ylabel('Coherence'), hold all
%         subplot(4,2,6)
%         plot(f,cLL_bip_off(:,j+1:end,j,i),'-',f,cLL_bip_on(:,j+1:end,j,i),'--')
%         title('Coh Bipolar'), xlim([min(f) max(f)]), ylim([0 1])
%         xlabel('f [Hz]'), ylabel('Coherence'), hold all
%     end
%     for j=1:Nemg(i)
%         subplot(4,2,7)
%         plot(f,cLE_tot_off(:,:,j,i),'-',f,cLE_tot_on(:,j+1:end,j,i),'--')
%         title('Coh LE Mono'), xlim([min(f) max(f)]), ylim([0 1])
%         xlabel('f [Hz]'), ylabel('Coherence'), hold all
%         subplot(4,2,8)
%         plot(f,cLE_bip_off(:,:,j,i),'-',f,cLE_bip_on(:,j+1:end,j,i),'--')
%         title('Coh LE Bipolar'), xlim([min(f) max(f)]), ylim([0 1])
%         xlabel('f [Hz]'), ylabel('Coherence'), hold all
%     end
%     print(fig, ['pow_coh_P' num2str(i,'%02i') 'w0_' num2str(w0,'%02i') ...
%         '_' type], '-dpsc');
%     saveas(fig, ['pow_coh_P' num2str(i,'%02i') 'w0_' num2str(w0,'%02i') ...
%         '_' type '.fig']);
%     close(fig)
%     
% end



    
% Boxplots
fig = figure('PaperUnits','centimeters', 'visi', 'off', 'PaperSize', [21 29.7]);

% Monopolar total
subplot(3,2,1)
boxplot( log10(PL_tot_on./PL_tot_off), ...
    'labels', {'\delta', '\theta', '\alpha', '\beta', '\beta1', '\beta2', ...
    '\gamma', 'HF'})
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex');
ylim([-1 1]), grid on
ylabel('log_{10}(P_{ON}/P_{OFF})')
title('Total monopolar')
hold all
for i=1:Nb
    if sig_tot(i)<0.01;        
        plot([i-0.05 i+0.05],[0.9 0.9],'*k');
    elseif sig_tot(i)<0.05;
        plot(i,0.9,'*k');
    end
end
% Monopolar inc
subplot(3,2,2)
boxplot( log10(PL_inc_on./PL_inc_off), ...
    'labels', {'\delta', '\theta', '\alpha', '\beta', '\beta1', '\beta2', ...
    '\gamma', 'HF'})
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex');
ylim([-1 1]), grid on
ylabel('log_{10}(P_{ON}/P_{OFF})')
title('Monopolar incoherent')
hold all
for i=1:Nb
    if sig_inc(i)<0.01;        
        plot([i-0.05 i+0.05],[0.9 0.9],'*k');
    elseif sig_inc(i)<0.05;
        plot(i,0.9,'*k');
    end
end
% Monopolar tot-vc
subplot(3,2,3)
boxplot( log10(PL_nvc_on./PL_nvc_off), ...
    'labels', {'\delta', '\theta', '\alpha', '\beta', '\beta1', '\beta2', ...
    '\gamma', 'HF'})
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex');
ylim([-1 1]), grid on
ylabel('log_{10}(P_{ON}/P_{OFF})')
title('Monopolar tot-vc')
hold all
for i=1:Nb
    if sig_nvc(i)<0.01;        
        plot([i-0.05 i+0.05],[0.9 0.9],'*k');
    elseif sig_nvc(i)<0.05;
        plot(i,0.9,'*k');
    end
end
% Monopolar coherent
subplot(3,2,4)
boxplot( log10(PL_coh_on./PL_coh_off), ...
    'labels', {'\delta', '\theta', '\alpha', '\beta', '\beta1', '\beta2', ...
    '\gamma', 'HF'})
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex');
ylim([-1 1]), grid on
ylabel('log_{10}(P_{ON}/P_{OFF})')
title('Monopolar coherent')
hold all
for i=1:Nb
    if sig_coh(i)<0.01;        
        plot([i-0.05 i+0.05],[0.9 0.9],'*k');
    elseif sig_coh(i)<0.05;
        plot(i,0.9,'*k');
    end
end
% Bipolar
subplot(3,2,5)
boxplot( log10(PL_bip_on./PL_bip_off), ...
    'labels', {'\delta', '\theta', '\alpha', '\beta', '\beta1', '\beta2', ...
    '\gamma', 'HF'})
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex');
ylim([-1 1]), grid on
ylabel('log_{10}(P_{ON}/P_{OFF})')
title('Total bipolar')
hold all
for i=1:Nb
    if sig_bip(i)<0.01;        
        plot([i-0.05 i+0.05],[0.9 0.9],'*k');
    elseif sig_bip(i)<0.05;
        plot(i,0.9,'*k');
    end
end
% Monopolar volume conduction
subplot(3,2,6)
boxplot( log10(PL_vc_on./PL_vc_off), ...
    'labels', {'\delta', '\theta', '\alpha', '\beta', '\beta1', '\beta2', ...
    '\gamma', 'HF'})
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex');
ylim([-1 1]), grid on
ylabel('log_{10}(P_{ON}/P_{OFF})')
title('Monopolar, volume cond.')
hold all
for i=1:Nb
    if sig_vc(i)<0.01;        
        plot([i-0.05 i+0.05],[0.9 0.9],'*k');
    elseif sig_vc(i)<0.05;
        plot(i,0.9,'*k');
    end
end

print(fig, ['pow_coh_w0_' num2str(w0,'%02i') '_' type '_box'], '-dpsc');
saveas(fig, ['pow_coh_w0_' num2str(w0,'%02i') '_' type '_box.fig']);
close(fig)



    
%% Boxplots of coherence
fig = figure('visi', 'off', 'PaperSize', [10 10]);

% Monopolar no vc
subplot(2,1,1)
boxplot( cLL_nvc_on_band-cLL_nvc_off_band, ...
    'labels', {'\delta', '\theta', '\alpha', '\beta', '\beta1', '\beta2', ...
    '\gamma', 'HF'})
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex');
ylim([-1 1]), grid on
ylabel('C_{ON}-C_{OFF})')
title('Monopolar - vc')
hold all
for i=1:Nb
    if sig_coh_nvc(i)<0.01;        
        plot([i-0.05 i+0.05],[0.9 0.9],'*k');
    elseif sig_coh_nvc(i)<0.05;
        plot(i,0.9,'*k');
    end
end
% Bipolar
subplot(2,1,2)
boxplot( cLL_bip_on_band-cLL_bip_off_band, ...
    'labels', {'\delta', '\theta', '\alpha', '\beta', '\beta1', '\beta2', ...
    '\gamma', 'HF'})
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex');
ylim([-1 1]), grid on
ylabel('C_{ON}-C_{OFF})')
title('Bipolar')
hold all
for i=1:Nb
    if sig_coh_bip(i)<0.01;        
        plot([i-0.05 i+0.05],[0.9 0.9],'*k');
    elseif sig_coh_bip(i)<0.05;
        plot(i,0.9,'*k');
    end
end

print(fig, ['pow_coh_w0_' num2str(w0,'%02i') '_' type '_cohbox'], '-dpsc');
saveas(fig, ['pow_coh_w0_' num2str(w0,'%02i') '_' type '_cohbox.fig']);
close(fig)
 

end %w0 for-loop end


%% Save Variables
save(['pow_coh_w0_' num2str(w0,'%02i') '_' type '.mat'], ...
    'cLL_bip_off', 'cLL_bip_on', 'cLL_nvc_off', 'cLL_nvc_on', ...
    'cLL_tot_off', 'cLL_tot_on', 'P_nvc_off', 'P_nvc_on', 'P_bip_off', ...
    'P_bip_on', 'P_tot_off', 'P_tot_on', 'P_coh_off', 'P_coh_on', ...
    'P_vc_off', 'P_vc_on', 'P_EMG_off', 'P_EMG_on', 'fband', 'type', 'Nex', ...
    'cLE_bip_off', 'cLE_bip_on', 'cLE_tot_off', 'cLE_tot_on', 'f', ...
    'PL_bip_off', 'PL_bip_on', 'PL_nvc_off', 'PL_nvc_on', 'PL_vc_off', ...
    'PL_tot_off', 'PL_tot_on', 'PL_coh_off', 'PL_coh_on', 'PL_vc_on', ...
    'sig_bip', 'sig_inc', 'sig_nvc', 'sig_tot', 'sig_coh', 'sig_vc', 'sig_coh_bip', ...
    'sig_coh_nvc', 'LFPact', 'Nb', 'cLL_bip_off_band', 'cLL_bip_on_band', ...
    'cLL_nvc_off_band', 'cLL_nvc_on_band', 'w0', 'scale', 't_smo')


clear a sig_bip sig_nvc sig_coh sig_inc sig_tot sig_vc LFP_ON LFP_OFF ...
    EMG_ON EMG_OFF Channel ArtLFP_OFF ArtLFP_ON ArtEMG_OFF ArtEMG_ON ...
    Patient Nch LFPact fband

end %type for loop end
