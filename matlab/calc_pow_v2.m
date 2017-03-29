%% Calculate mono- and bipolar changes in PSD of LFP

function calc_pow_v2(task)

    %% Load LFP data
    load(['data_' task '_v3']);

    
    %% Parameters
    f       = linspace(1,100,100);
    Nf      = length(f);
    w0      = 12;
    nsig    = 6;
    sig     = sig_coh_thresh(w0, nsig);
    phthres = sig_phi_thresh(w0, nsig);
    usenan  = 0; % 1 = take mean only over times where signal is present
    ds      = 5;
    tag     = 'v3_bip';
    STNact  = 1;
    Nex     = length(LFP_OFF);
    zero    = false; % define volume-conduction by 0°--phthres only (not 180°-phthres--180°)


    %% Calculate spectral powers: total, coherent, incoherent and volume
    %% conducted
    if STNact>0
        nch  = 3;
        nch2 = 3;
    else
        nch  = 4;
        nch2 = 6;
    end

    P_inc_off    = NaN(Nf,nch,Nex);
    P_coh_off    = NaN(Nf,nch,Nex);
    P_vc_off     = NaN(Nf,nch,Nex);
    P_bip_off    = NaN(Nf,nch2,Nex);

    P_inc_on     = NaN(Nf,nch,Nex);
    P_coh_on     = NaN(Nf,nch,Nex);
    P_vc_on      = NaN(Nf,nch,Nex);
    P_bip_on     = NaN(Nf,nch2,Nex);

    scale   = (w0+sqrt(2+w0^2))/4/pi ./ f;

    for i=1:Nex
        fprintf('Processing data set %u/%u.\n', i, Nex);

        % Check for STN activity
        ind = find(LFPact{i}>=STNact);
        Nch(i) = length(ind);
        LFPch{i} = Channel{i}(ind);
        lfp = LFP_OFF{i}(:,LFPact{i}>=STNact);
        clear art
        k = 1;
        for j=ind
            art{k} = ArtLFP_OFF{i}{j};
            k = k+1;
        end
        
        % Determine bipolar combinations w/o pure nonSTN channels
        combo = nchoosek(1:Nch(i),2);
        Nbc(i)= size(combo,1);
        BiPolCh{i} = LFPch{i}(combo);


        %% OFF
        % LFP monopolar
        [~,W,coi] = procdata(lfp, 'freq', f, 'w0', w0, ...
            'filter', [], 'art', art); %'filter', [44 56; 90 110], 
        [Cs, Wxy, W] = wave_cohere ( W, scale, nsig, ds );
        coi = coi(1:ds:end,:)/nsig;
        [Pcoh, Pinc, Pvc] = pcc ( f, W, Cs, coi, sig, Wxy, usenan, phthres, zero );
        P_coh_off(:,1:Nch(i),i) = squeeze(nanmean(Pcoh,3));
        P_inc_off(:,1:Nch(i),i) = squeeze(nanmean(Pinc,3));
        P_vc_off(:,1:Nch(i),i)  = squeeze(nanmean(Pvc,3));
        
        % LFP bipolar
        x = lfp(:,combo(:,1))-lfp(:,combo(:,2));
        a = cell(Nbc(i),1);
        for j=1:Nbc(i);
            a{j} = [art{combo(j,1)}; art{combo(j,2)}];
        end
        [~,~,~,P_bip_off(:,1:Nbc(i),i)] = procdata(x, 'freq', f, 'w0', w0, ...
            'filter', [], 'art', a);
                

        %% ON
        
        % LFP monopolar
        lfp = LFP_ON{i}(:,LFPact{i}>=STNact);
        clear art
        k = 1;
        for j=ind
            art{k} = ArtLFP_ON{i}{j};
            k = k+1;
        end
        [~,W,coi] = procdata(lfp, 'freq', f, 'w0', w0, ...
            'filter', [], 'art', art);
        [Cs, Wxy, W] = wave_cohere ( W, scale, nsig, ds );
        coi = coi(1:ds:end,:)/nsig;
        [Pcoh, Pinc, Pvc] = pcc ( f, W, Cs, coi, sig, Wxy, usenan, phthres, zero );
        P_coh_on(:,1:Nch(i),i) = squeeze(nanmean(Pcoh,3));
        P_inc_on(:,1:Nch(i),i) = squeeze(nanmean(Pinc,3));
        P_vc_on(:,1:Nch(i),i)  = squeeze(nanmean(Pvc,3));
        
        % LFP bipolar
        x = lfp(:,combo(:,1))-lfp(:,combo(:,2));
        a = cell(Nbc(i),1);
        for j=1:Nbc(i);
            a{j} = [art{combo(j,1)}; art{combo(j,2)}];
        end
        [~,~,~,P_bip_on(:,1:Nbc(i),i)] = procdata(x, 'freq', f, 'w0', w0, ...
            'filter', [], 'art', a);


    end

    P_tot_off = P_inc_off + P_coh_off + P_vc_off;
    P_tot_on  = P_inc_on  + P_coh_on  + P_vc_on;
    
    clear x a W Wxy C coi1 Cs W1s Wxys coi1s



    %% Save Variables
    save(['LFP_pow_' tag '_' task '.mat'], 'BiPolCh', 'Nbc', ...
        'P_tot_on', 'P_tot_off', 'STNact', 'Nch', 'P_bip_off', 'P_bip_on', ...
        'P_coh_off', 'P_coh_on', 'P_vc_off', 'P_vc_on', 'LFPch', ...
        'task', 'Nex', 'f', 'P_inc_off', 'P_inc_on', 'Patient', ...
        'LFPact', 'w0', 'scale', 'nsig', 'usenan', 'ds', 'rigor')
    
    
%% Plot
% load xyz
% B(:,:,3)=squeeze(nanmean(P_vc_on./P_tot_on-P_vc_off./P_tot_off,2));
% B(:,:,1)=squeeze(nanmean(P_inc_on./P_tot_on-P_inc_off./P_tot_off,2));
% B(:,:,2)=squeeze(nanmean(P_coh_on./P_tot_on-P_coh_off./P_tot_off,2));
% M = [mean(B(:,:,1),2) mean(B(:,:,2),2) mean(B(:,:,3),2)];
% E = [std(B(:,:,1),[],2) std(B(:,:,2),[],2) std(B(:,:,3),[],2)];
% mseb(f,M',E'/sqrt(8))