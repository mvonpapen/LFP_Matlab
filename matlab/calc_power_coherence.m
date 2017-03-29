%% Test mono- and bipolar PSD changes and coherence from OFF to ON

function calc_power_coherence(task)

%% Load LFP data
% mypool = parpool(4);
% addAttachedFiles(mypool, {'smooth_mat.m', 'wave_cohere.m', 'acoh1d_uptri.m'});


    load(['data_' task '_v3']);


% % task = {'Fist', 'Hold', 'Rest'};
% % for ind = 1:3
% %     type = task{ind};
% %     LFP_OFF = load(['data_' type], 'LFP_OFF');
% %     LFP_ON  = load(['data_' type], 'LFP_ON');
% %     EMG_OFF = load(['data_' type], 'EMG_OFF');
% %     EMG_ON  = load(['data_' type], 'EMG_OFF');
%     % 
%     % [t, dat, CH, art, rigor, ndat, x, LFPactivity] = ...
%     %     load_dat_from_DATA ( 'experiment', type, 'activity', 0, 'commontime', 1, 'noEMG', 0);
%     % 
%     % load DATA_last.mat

    %% Parameters
    f       = [1:45:5:100];
    Nf      = length(f);
    phthres = 15;
    w0      = 12;
    nsig    = 6;
    sig     = sig_coh_thresh(w0, nsig);
    usenan  = 0; % 1 = take mean only over times where signal is present
    ds      = 5;
    tag     = 'v160416';
    Nex     = length(LFP_OFF);
    zero    = true; % define volume-conduction by 0°--phthres only (not 180°-phthres--180°)


    %% Calculate spectral powers: total, coherent, incoherent and volume
    %% conducted
    nch = max(Nch);

    P_EMG_off    = NaN(Nf, 3 ,Nex);
    P_tot_off    = NaN(Nf,nch,Nex);
    P_vc_off     = NaN(Nf,nch,Nex);
    P_nvc_off    = NaN(Nf,nch,Nex);
    P_inc_off    = NaN(Nf,nch,Nex);
    P_coh_off    = NaN(Nf,nch,Nex);
    P_bip_off    = NaN(Nf, 6 ,Nex);
    cLL_tot_off  = NaN(Nf,nch,nch,Nex);
    cLL_nvc_off  = NaN(Nf,nch,nch,Nex);
    cLL_bip_off  = NaN(Nf, 6 , 6 ,Nex);
    cLE_tot_off  = NaN(Nf,nch,3,Nex);
    cLE_bip_off  = NaN(Nf, 6 ,3,Nex);

    P_EMG_on     = NaN(Nf, 3 ,Nex);
    P_tot_on     = NaN(Nf,nch,Nex);
    P_vc_on      = NaN(Nf,nch,Nex);
    P_nvc_on     = NaN(Nf,nch,Nex);
    P_coh_on     = NaN(Nf,nch,Nex);
    P_inc_on     = NaN(Nf,nch,Nex);
    P_bip_on     = NaN(Nf, 6 ,Nex);
    cLL_tot_on   = NaN(Nf,nch,nch,Nex);
    cLL_nvc_on   = NaN(Nf,nch,nch,Nex);
    cLL_bip_on   = NaN(Nf, 6 , 6 ,Nex);
    cLE_tot_on   = NaN(Nf,nch,3,Nex);
    cLE_bip_on   = NaN(Nf, 6 ,3,Nex);

    STNelec      = NaN(6,Nex);
    Nbc          = NaN(Nex,1);
    Nemg         = NaN(Nex,1);
    BiPolCh      = cell(Nex,1);

    scale   = (w0+sqrt(2+w0^2))/4/pi ./ f;

    for i=1:Nex
        fprintf('Processing data set %u/%u.\n', i, Nex);
        Nemg(i)=size(EMG_OFF{i},2);

        % Determine bipolar combinations w/o pure nonSTN channels
        combo = nchoosek(1:Nch(i),2);
        STNelec(1:size(combo,1),i) = sum(LFPact{i}(combo)>0,2);
    %     combo = combo( STNelec(:,i)>0, : ); %no bipolar PSD of nonSTN electrodes
        Nbc(i)= size(combo,1);
        BiPolCh{i} = Channel{i}(combo);



        %% OFF

        % LFP monopolar
        [~,W1,coi1,P_tot_off(:,1:Nch(i),i)] = procdata(LFP_OFF{i}, 'freq', f, 'w0', w0, ...
            'filter', [44 56; 90 110], 'art', ArtLFP_OFF{i}); %'filter', [44 56; 90 110], 
        [Cs, Wxys, W1s] = wave_cohere ( W1, scale, nsig, ds );
        coi1s = coi1(1:ds:end,:)/nsig;
        [cLL_tot_off(:,1:Nch(i),1:Nch(i),i), cLL_nvc_off(:,1:Nch(i),1:Nch(i),i)] ...
                = acoh1d_uptri( f, Cs, coi1s, Wxys, phthres, sig );
        [Pcoh, Pinc, Pvc, Pnvc] = psd_acoh ( f, W1s, Cs, coi1s, sig, Wxys, usenan, phthres, zero );
        P_coh_off(:,1:Nch(i),i) = squeeze(nanmean(Pcoh,3));
        P_inc_off(:,1:Nch(i),i) = squeeze(nanmean(Pinc,3));
        P_vc_off(:,1:Nch(i),i)  = squeeze(nanmean(Pvc,3));
        P_nvc_off(:,1:Nch(i),i) = squeeze(nanmean(Pnvc,3));

        % LFP bipolar
        x = LFP_OFF{i}(:,combo(:,1))-LFP_OFF{i}(:,combo(:,2));
        a = cell(Nbc(i),1);
        for j=1:Nbc(i);
            a{j} = [ArtLFP_OFF{i}{combo(j,1)}; ArtLFP_OFF{i}{combo(j,2)}];
        end
        [~,W2,coi2,P_bip_off(:,1:Nbc(i),i)] = procdata(x, 'freq', f, 'w0', w0, ...
            'filter', [44 56; 90 110], 'art', a);
        [Cs, Wxys] = wave_cohere ( W2, scale, nsig, ds );
        coi2s = coi2(1:ds:end,:)/nsig;
        cLL_bip_off(:,1:Nbc(i),1:Nbc(i),i) ...
            = acoh1d_uptri( f, Cs, coi2s, Wxys, phthres, sig );      
        
%         % EMG
%         if Nemg(i)>0
%             [~,WE,coiE,P_EMG_off(:,1:Nemg(i),i)] = procdata(EMG_OFF{i}, 'freq', f, 'w0', w0, ...
%                 'filter', [0 60; 90 110], 'art', ArtEMG_OFF{i}, 'rect', 1);
%             Cs = wave_coh3 ( W1, WE, scale, nsig, ds );
%             coiEs = coiE(1:ds:end,:)/nsig;
%             cLE_tot_off(:,1:Nch(i),1:Nemg(i),i) = coh1d( f, Cs, coi1s, coiEs );
%             Cs = wave_coh3 ( W2, WE, scale, nsig, ds );
%             cLE_bip_off(:,1:Nbc(i),1:Nemg(i),i) = coh1d( f, Cs, coi2s, coiEs );
%         end


        %% ON
        
        % LFP monopolar
        [~,W1,coi1,P_tot_on(:,1:Nch(i),i)] = procdata(LFP_ON{i}, 'freq', f, 'w0', w0, ...
            'filter', [44 56; 90 110], 'art', ArtLFP_ON{i});
        [Cs, Wxys, W1s] = wave_cohere ( W1, scale, nsig, ds );
        coi1s = coi1(1:ds:end,:)/nsig;
        [cLL_tot_on(:,1:Nch(i),1:Nch(i),i), cLL_nvc_on(:,1:Nch(i),1:Nch(i),i)] ...
                = acoh1d_uptri( f, Cs, coi1s, Wxys, phthres, sig );
        [Pcoh, Pinc, Pvc, Pnvc] = psd_acoh ( f, W1s, Cs, coi1s, sig, Wxys, usenan, phthres, zero );
        P_coh_on(:,1:Nch(i),i) = squeeze(nanmean(Pcoh,3));
        P_inc_on(:,1:Nch(i),i) = squeeze(nanmean(Pinc,3));
        P_vc_on(:,1:Nch(i),i)  = squeeze(nanmean(Pvc,3));
        P_nvc_on(:,1:Nch(i),i) = squeeze(nanmean(Pnvc,3));

        % LFP bipolar
        x = LFP_ON{i}(:,combo(:,1))-LFP_ON{i}(:,combo(:,2));
        a = cell(Nbc(i),1);
        for j=1:Nbc(i);
            a{j} = [ArtLFP_ON{i}{combo(j,1)}; ArtLFP_ON{i}{combo(j,2)}];
        end
        [~,W2,coi2,P_bip_on(:,1:Nbc(i),i)] = procdata(x, 'freq', f, ...
            'w0', w0, 'filter', [44 56; 90 110], 'art', a);
        coi2s = coi2(1:ds:end,:)/nsig;
        [Cs, Wxys] = wave_cohere ( W2, scale, nsig, ds );
        cLL_bip_on(:,1:Nbc(i),1:Nbc(i),i) = acoh1d_uptri( f, Cs, coi2s, Wxys, phthres, sig );

%         % EMG
%         if Nemg(i)>0
%             [~,WE,coiE,P_EMG_on(:,1:Nemg(i),i)] = procdata(EMG_ON{i}, 'freq', f, 'w0', w0, ...
%                 'filter', [0 60; 90 110], 'art', ArtEMG_ON{i}, 'rect', 1);
%             coiEs = coiE(1:ds:end,:)/nsig;
%             Cs = wave_coh3 ( W1, WE, scale, nsig, ds );
%             cLE_tot_on(:,1:Nch(i),1:Nemg(i),i) = coh1d( f, Cs, coi1s, coiEs );
%             Cs = wave_coh3 ( W2, WE, scale, nsig, ds );
%             cLE_bip_on(:,1:Nbc(i),1:Nemg(i),i) = coh1d( f, Cs, coi2s, coiEs );
%         end

    end

    clear x a W1 W2 WE Wxy C coi1 coi2 coiE combo Cs W1s W2s Wxys coi1s coi2s coiEs



    %% Save Variables
    save(['pow_coh_' tag '_w0_' num2str(w0,'%02i') '_' task '.mat'], ...
        'cLL_bip_off', 'cLL_bip_on', 'cLL_nvc_off', 'cLL_nvc_on', ...
        'cLL_tot_off', 'cLL_tot_on', 'P_nvc_off', 'P_nvc_on', 'P_bip_off', ...
        'P_bip_on', 'P_tot_off', 'P_tot_on', 'P_coh_off', 'P_coh_on', ...
        'P_vc_off', 'P_vc_on', 'P_EMG_off', 'P_EMG_on', 'task', 'Nex', ...
        'cLE_bip_off', 'cLE_bip_on', 'cLE_tot_off', 'cLE_tot_on', 'f', ...
        'P_inc_off', 'P_inc_on', 'Patient', 'LFPact', 'w0', 'scale', ...
        'nsig', 'usenan', 'ds', 'rigor', 'Nch', 'Nbc', 'BiPolCh', 'Nemg')


%     clear a sig_bip sig_nvc sig_coh sig_inc sig_tot sig_vc LFP_ON LFP_OFF ...
%         EMG_ON EMG_OFF Channel ArtLFP_OFF ArtLFP_ON ArtEMG_OFF ArtEMG_ON ...
%         Patient Nch LFPact fband

% end %type for loop end

% delete(mypool)
