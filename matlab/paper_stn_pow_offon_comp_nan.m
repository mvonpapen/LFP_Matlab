%% Test monopolar PSD changes from OFF to ON for incoherent, coherent and
%% volume conduction signals using average power (usenan=1)


%% Load LFP data

for type = {'Ruhe'};
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
usenan  = 1;
w0      = 12;




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
for i=1:Nex
    Channel{i}=CH{iend(i)};
    iLFP = find(ismember(Channel{i}, LFPch));
    
    LFP_OFF{i}=dat{i0(i)}(:,iLFP);
    LFP_ON{i}=dat{iend(i)}(:,iLFP);
    
    for j=1:length(iLFP)
        ArtLFP_OFF{i}{j}=art{i0(i)}{iLFP(j)};
        ArtLFP_ON{i}{j}=art{iend(i)}{iLFP(j)};
    end
        
    Patient{i}=DATA(ndat(i0(i))).patient;
    Patient2{i}=DATA(ndat(iend(i))).patient;
    Nch(i) = length(iLFP);
    LFPact{i}=LFPactivity{iend(i)};
    
    
    if ~strcmp(Patient{i},Patient2{i})
        error('ERROR: i0 and iend refer to different patients for i=%i!', i)
    end
end
rigor=rigor(iend);

clear i i0 iend DATA dat art CH LFPactivity t ndat pat Patient2 x




%% Calculate spectral powers: total, coherent, incoherent and volume
%% conducted
nch = max(Nch);

P_vc_off     = NaN(Nf,nch,Nex);
P_inc_off    = NaN(Nf,nch,Nex);
P_coh_off    = NaN(Nf,nch,Nex);

P_vc_on      = NaN(Nf,nch,Nex);
P_coh_on     = NaN(Nf,nch,Nex);
P_inc_on     = NaN(Nf,nch,Nex);

scale   = (w0+sqrt(2+w0^2))/4/pi ./ f;


for i=1:Nex
    [i Nex]
    
    
    %% OFF
    [x,W,coi] = procdata(LFP_OFF{i}, 'freq', f, 'w0', w0, ...
        'filter', [], 'art', ArtLFP_OFF{i});
    [C, Wxy] = wave_acoh2 ( W, scale, t_smo, sc_smo );
    [Pcoh, Pinc, Pvc] = psd_acoh ( f, W, C, coi, sig, Wxy, usenan );
    P_vc_off(:,1:Nch(i),i)  = squeeze(nanmean(Pvc,3));
    P_coh_off(:,1:Nch(i),i) = squeeze(nanmean(Pcoh,3));
    P_inc_off(:,1:Nch(i),i) = squeeze(nanmean(Pinc,3));

    
    %% ON
    [x,W,coi] = procdata(LFP_ON{i}, 'freq', f, 'w0', w0, ...
        'filter', [], 'art', ArtLFP_ON{i});
    [C, Wxy] = wave_acoh2 ( W, scale, t_smo, sc_smo );
    [Pcoh, Pinc, Pvc] = psd_acoh ( f, W, C, coi, sig, Wxy, usenan );
    P_vc_on(:,1:Nch(i),i)  = squeeze(nanmean(Pvc,3));
    P_coh_on(:,1:Nch(i),i) = squeeze(nanmean(Pcoh,3));
    P_inc_on(:,1:Nch(i),i) = squeeze(nanmean(Pinc,3));
    
end

clear i x W Wxy C coi





%% Total power in frequency bands

fband{1} = find(f>=1  & f<=4 );
fband{2} = find(f>=4  & f<=7 );
fband{3} = find(f>=7  & f<=13);
fband{4} = find(f>=13 & f<=30);
fband{5} = find(f>=13 & f<=20);
fband{6} = find(f>=20 & f<=30);
fband{7} = find(f>=60 & f<=90);
fband{8} = find(f>=250 & f<=350);
Nb       = length(fband);

PL_coh_off = NaN(nch*Nex, Nb);
PL_coh_on  = NaN(nch*Nex, Nb);
PL_inc_off = NaN(nch*Nex, Nb);
PL_inc_on  = NaN(nch*Nex, Nb);
PL_vc_off  = NaN(nch*Nex, Nb);
PL_vc_on   = NaN(nch*Nex, Nb);


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


% Wilcoxon sign-rank test
for i=1:Nb
    sig_coh(i) = signrank ( log10(PL_coh_off(:,i)./PL_coh_on(:,i) ) );
    sig_vc(i)  = signrank ( log10(PL_vc_off(:,i) ./PL_vc_on(:,i)  ) );
    if ~any(~isnan(PL_inc_off(:,i)./PL_inc_on(:,i)))
        sig_inc(i)=1;
    else
        sig_inc(i) = signrank ( log10(PL_inc_off(:,i)./PL_inc_on(:,i) ) );
    end
end
%




%% Plot results for each patient

fig = figure('PaperUnits','centimeters', 'visi', 'off', 'PaperSize', [12 15]);
%     'PaperPosition', [1 1 10 13], ...
    %, ...
%     'Position',[0 0 12 15]/2);
    
% Monopolar inc
subplot(3,1,1)
boxplot( log10(PL_inc_on./PL_inc_off), ...
    'labels', {'\delta', '\theta', '\alpha', '\beta', '\beta1', '\beta2', ...
    '\gamma', 'HF'})
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex');
ylim([-1 1]), grid on
ylabel('log_{10}(P_{ON}/P_{OFF})')
title('Incoherent')
hold all
for i=1:Nb
    if sig_inc(i)<0.01;        
        plot([i-0.05 i+0.05],[0.9 0.9],'*k');
    elseif sig_inc(i)<0.05;
        plot(i,0.9,'*k');
    end
end
% Monopolar coherent
subplot(3,1,2)
boxplot( log10(PL_coh_on./PL_coh_off), ...
    'labels', {'\delta', '\theta', '\alpha', '\beta', '\beta1', '\beta2', ...
    '\gamma', 'HF'})
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex');
ylim([-1 1]), grid on
ylabel('log_{10}(P_{ON}/P_{OFF})')
title('Coherent')
hold all
for i=1:Nb
    if sig_coh(i)<0.01;        
        plot([i-0.05 i+0.05],[0.9 0.9],'*k');
    elseif sig_coh(i)<0.05;
        plot(i,0.9,'*k');
    end
end
% Monopolar volume conduction
subplot(3,1,3)
boxplot( log10(PL_vc_on./PL_vc_off), ...
    'labels', {'\delta', '\theta', '\alpha', '\beta', '\beta1', '\beta2', ...
    '\gamma', 'HF'})
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex');
ylim([-1 1]), grid on
ylabel('log_{10}(P_{ON}/P_{OFF})')
title('Volume Conduction')
hold all
for i=1:Nb
    if sig_vc(i)<0.01;        
        plot([i-0.05 i+0.05],[0.9 0.9],'*k');
    elseif sig_vc(i)<0.05;
        plot(i,0.9,'*k');
    end
end

print(fig, ['pow_w0_' num2str(w0,'%02i') '_' type '_nan_box'], '-dpsc');
saveas(fig, ['pow_w0_' num2str(w0,'%02i') '_' type '_nan_box.fig']);
close(fig)



%% Save Variables
save(['pow_w0_' num2str(w0,'%02i') '_' type '_nan.mat'], ...
    'P_inc_off', 'P_inc_on', 'P_coh_off', 'P_coh_on', 'P_vc_off', 'P_vc_on', ...
    'PL_inc_off', 'PL_inc_on', 'PL_inc_off', 'PL_inc_on', 'PL_vc_off', 'PL_vc_on', ...
    'f', 'fband', 'sig_inc', 'sig_coh', 'sig_vc', 'type', ...
    'LFPact', 'Nb', 'Nex', 'w0', 'scale', 't_smo')


clear sig_coh sig_inc sig_vc LFP_ON LFP_OFF Channel ArtLFP_OFF ArtLFP_ON ...
    Patient Nch LFPact fband

end %type for-loop end
