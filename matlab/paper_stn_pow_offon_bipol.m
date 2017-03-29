%% Test bipolar(!)-PSD changes and coherence from OFF to ON


%% Load LFP data
[x, dat, CH, art, rigor, ndat, x, LFPactivity] = ...
    load_dat_from_DATA ( 'experiment', 'Faust', 'activity', 0, 'commontime', 1, 'noEMG', 1);

load DATA_last.mat

%% Parameters
Nf  = 60;
f   = logspace(0, log10(350), Nf);
sig = 0.495;
phthres = 10;

%% Determine sites and analysis
% % Sites with rest ON OFF, Nelec>1 and rigor>=30
% i0   = [1  6 11 26 32 41 65]; %last entry: no EMG channels
% iend = [5 10 17 31 68 46 67];
% ref  = [1  3  1  2  3  2  1]; %hand-picked reference electrodes (act=0) 
% % Sites with hold ON OFF and nLFP>1
% i0   = [1 3 5  9 11 21];
% iend = [2 4 6 10 23 22];
% ref  = [1 3 1  2  3  1]; %hand-picked reference electrodes (act=0) 
% Sites with fist ON OFF and rigor>=30
i0   = [1 3 8 10 18];
iend = [2 4 9 20 19];
ref  = [1 3 2  3  1]; %hand-picked reference electrodes (act=0)

Nex = length(i0);

LFPch = {'C', 'L', 'A', 'M', 'P'};
for i=1:Nex
    Channel{i}=CH{i0(i)};
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

clear i i0 iend DATA dat art CH LFPactivity t ndat pat Patient2



%% Calculate spectral powers: total, coherent, incoherent and volume
%% conducted
P_off        = NaN(Nf,6,Nex);
P_on         = NaN(Nf,6,Nex);
STNelec      = NaN(6,Nex);
Nbc          = NaN(Nex,1);
BiPolCh      = cell(Nex);

for i=1:Nex
    
    % Determine bipolar combinations w/o pure nonSTN channels
%     combo = nchoosek(1:Nch(i),2);
    combo = setdiff(1:Nch(i), ref(i));
    combo = [combo' ones(length(combo),1)*ref(i)];
    STNelec(1:size(combo,1),i) = sum(LFPact{i}(combo)>0,2);
    combo = combo( STNelec(:,i)>0, : ); %no bipolar PSD of nonSTN electrodes
    Nbc(i)= size(combo,1);
    BiPolCh{i} = Channel{i}(combo);
    
    %% OFF
    x = LFP_OFF{i}(:,combo(:,1))-LFP_OFF{i}(:,combo(:,2));
    for j=1:Nbc(i);
        a{j} = [ArtLFP_OFF{i}{combo(j,1)}; ArtLFP_OFF{i}{combo(j,2)}];
    end
    [x,W,coi,P_off(:,1:Nbc(i),i)] = procdata(x, 'freq', f, 'filter', [], 'art', a);
    %% ON
    x = LFP_ON{i}(:,combo(:,1))-LFP_ON{i}(:,combo(:,2));
    for j=1:Nbc(i);
        a{j} = [ArtLFP_ON{i}{combo(j,1)}; ArtLFP_ON{i}{combo(j,2)}];
    end
    [x,W,coi,P_on(:,1:Nbc(i),i)] = procdata(x, 'freq', f, 'filter', [], 'art', a);
end

clear LFP_OFF LFP_ON ArtOFF ArtON x a W Wxy C coi combo


%% Reshape matrices and one-sided coherence
P_bp_off = P_off(:,:);
P_bp_on  = P_on(:,:);

    

%% Significant power of coherence changes?
for i=1:60;
    sig_pow(i)=signrank( P_bp_off(i,:) ./ P_bp_on(i,:), 1 );
end


% %% Plot results
% figure
% subplot(3,3,1)
% loglog( f, pow3_off ./ pow3_on )
% title('Total power change during wrist'), ylim([0.05 20])
% ylabel('P_{OFF}/P_{ON}'), xlabel('f [Hz]')
% hold all
% plot( f, nanmean( mpsd_tot_off ./ mpsd_tot_on, 2 ), '--k', 'linew', 2 )
% i= sig_tot<0.05; if any(i); plot(f(i),10,'*k','markers',5); end
% i= sig_tot<0.01; if any(i); plot(f(i),13,'+k','markers',5); end
% 
% subplot(2,1,2)
% loglog(f,(mpsd_tot_off-mpsd_vc_off) ./ ...
%     (mpsd_tot_on-mpsd_vc_on))
% title('Change of power w/o vol. cond. during wrist'), ylim([0.05 20])
% ylabel('P_{OFF}/P_{ON}'), xlabel('f [Hz]')
% hold all
% plot(f, nanmean((mpsd_tot_off-mpsd_vc_off) ./ ...
%     (mpsd_tot_on-mpsd_vc_on),2), '--k', 'linew', 2)
% i= sig_novc<0.05; if any(i); plot(f(i),10,'*k','markers',5); end
% i= sig_novc<0.01; if any(i); plot(f(i),13,'+k','markers',5); end