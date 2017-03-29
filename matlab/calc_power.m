%% Calculate mono- and bipolar changes in PSD of LFP

function calc_power(task)

    %% Load LFP data
    load(['data_' task '_v3']);

    
    %% Parameters
    f       = 1:100;
    Nf      = length(f);
    phthres = 15;
    w0      = 12;
    nsig    = 6;
    sig     = sig_coh_thresh(w0, nsig);
    usenan  = 0; % 1 = take mean only over times where signal is present
    ds      = 5;
    tag     = 'v1';
    Nex     = length(LFP_OFF);
    zero    = false; % define volume-conduction by 0°--phthres only (not 180°-phthres--180°)


    %% Calculate spectral powers: total, coherent, incoherent and volume
    %% conducted
    nch = max(Nch);

    P_vc_off     = NaN(Nf,nch,Nex);
    P_inc_off    = NaN(Nf,nch,Nex);
    P_coh_off    = NaN(Nf,nch,Nex);
    P_bip_off    = NaN(Nf, 6 ,Nex);

    P_vc_on      = NaN(Nf,nch,Nex);
    P_coh_on     = NaN(Nf,nch,Nex);
    P_inc_on     = NaN(Nf,nch,Nex);
    P_bip_on     = NaN(Nf, 6 ,Nex);

    STNelec      = NaN(6,Nex);
    Nbc          = NaN(Nex,1);
    BiPolCh      = cell(Nex,1);

    scale   = (w0+sqrt(2+w0^2))/4/pi ./ f;

    for i=1:Nex
        fprintf('Processing data set %u/%u.\n', i, Nex);

        % Determine bipolar combinations w/o pure nonSTN channels
        combo = nchoosek(1:Nch(i),2);
        STNelec(1:size(combo,1),i) = sum(LFPact{i}(combo)>0,2);
        Nbc(i)= size(combo,1);
        BiPolCh{i} = Channel{i}(combo);



        %% OFF

        % LFP monopolar
        [~,W,coi] = procdata(LFP_OFF{i}, 'freq', f, 'w0', w0, ...
            'filter', [], 'art', ArtLFP_OFF{i}); %'filter', [44 56; 90 110], 
        [Cs, Wxy, W] = wave_cohere ( W, scale, nsig, ds );
        coi = coi(1:ds:end,:)/nsig;
        [Pcoh, Pinc, Pvc] = pcc ( f, W, Cs, coi, sig, Wxy, usenan, phthres, zero );
        P_coh_off(:,1:Nch(i),i) = squeeze(nanmean(Pcoh,3));
        P_inc_off(:,1:Nch(i),i) = squeeze(nanmean(Pinc,3));
        P_vc_off(:,1:Nch(i),i)  = squeeze(nanmean(Pvc,3));

        % LFP bipolar
        x = LFP_OFF{i}(:,combo(:,1))-LFP_OFF{i}(:,combo(:,2));
        a = cell(Nbc(i),1);
        for j=1:Nbc(i);
            a{j} = [ArtLFP_OFF{i}{combo(j,1)}; ArtLFP_OFF{i}{combo(j,2)}];
        end
        [~,~,~,P_bip_off(:,1:Nbc(i),i)] = procdata(x, 'freq', f, 'w0', w0, ...
            'filter', [], 'art', a);
        
        

        %% ON
        
        % LFP monopolar
        [~,W,coi] = procdata(LFP_ON{i}, 'freq', f, 'w0', w0, ...
            'filter', [], 'art', ArtLFP_ON{i});
        [Cs, Wxy, W] = wave_cohere ( W, scale, nsig, ds );
        coi = coi(1:ds:end,:)/nsig;
        [Pcoh, Pinc, Pvc] = pcc ( f, W, Cs, coi, sig, Wxy, usenan, phthres, zero );
        P_coh_on(:,1:Nch(i),i) = squeeze(nanmean(Pcoh,3));
        P_inc_on(:,1:Nch(i),i) = squeeze(nanmean(Pinc,3));
        P_vc_on(:,1:Nch(i),i)  = squeeze(nanmean(Pvc,3));

        % LFP bipolar
        x = LFP_ON{i}(:,combo(:,1))-LFP_ON{i}(:,combo(:,2));
        a = cell(Nbc(i),1);
        for j=1:Nbc(i);
            a{j} = [ArtLFP_ON{i}{combo(j,1)}; ArtLFP_ON{i}{combo(j,2)}];
        end
        [~,~,~,P_bip_on(:,1:Nbc(i),i)] = procdata(x, 'freq', f, ...
            'w0', w0, 'filter', [], 'art', a);


    end

    P_tot_off = P_inc_off + P_coh_off + P_vc_off;
    P_tot_on  = P_inc_on  + P_coh_on  + P_vc_on;
    
    clear x a W1 W2 WE Wxy C coi1 combo Cs W1s Wxys coi1s



    %% Save Variables
    save(['LFP_pow_' tag '_' task '.mat'], ...
        'P_bip_off', 'P_bip_on', 'P_tot_on', 'P_tot_off', ...
        'P_coh_off', 'P_coh_on', 'P_vc_off', 'P_vc_on', ...
        'task', 'Nex', 'f', 'P_inc_off', 'P_inc_on', 'Patient', ...
        'LFPact', 'w0', 'scale', 'nsig', 'usenan', 'ds', 'rigor', ...
        'Nch', 'Nbc', 'BiPolCh')