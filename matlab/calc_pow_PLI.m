%% Calculate mono- and bipolar changes in PSD of LFP

function calc_pow_PLI(task)

    %% Load LFP data
    load(['data_' task '_v3']);

    
    %% Parameters
    f       = linspace(1,100,100);
    Nf      = length(f);
    w0      = 12;
    nsig    = 6;
    tag     = 'v1_filt';
    STNact  = 1;
    Nex     = length(LFP_OFF);
    wPLI    = true;
    sig     = 0.705;
    filt    = [45 55; 90 110];


    %% Calculate spectral powers: total, coherent, incoherent and volume
    %% conducted
    if STNact>0
        nch  = 3;
    else
        nch  = 4;
    end

    P_inc_off    = NaN(Nf,nch,Nex);
    P_coh_off    = NaN(Nf,nch,Nex);

    P_inc_on     = NaN(Nf,nch,Nex);
    P_coh_on     = NaN(Nf,nch,Nex);

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


        %% OFF
        % LFP monopolar
        [~,W,coi] = procdata(lfp, 'freq', f, 'w0', w0, ...
            'filter', filt, 'art', art); %'filter', [44 56; 90 110], 
        coi = coi(1:ds:end,:)/nsig;
        Wxy = W.*conj(W);
        PLI = zeros(size(Wxy));
        for j=1:nch-1
            for k=j+1:nch
                PLI(:,:,j,k) = phase_lag_index(imag(Wxy(:,:,j,k)), ...
                    scale, nsig, dt, 'wPLI', wPLI);
                PLI(:,:,k,j) = PLI(:,:,j,k);
            end
        end
        PLIthreshd = PLI>sig;
        [Pcoh, Pinc] = pcc ( f, W, PLIthreshd, coi );
        P_coh_off(:,1:Nch(i),i) = squeeze(nanmean(Pcoh,3));
        P_inc_off(:,1:Nch(i),i) = squeeze(nanmean(Pinc,3));
                

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
            'filter', filt, 'art', art);
        coi = coi(1:ds:end,:)/nsig;
        Wxy = W.*conj(W);
        PLI = zeros(size(Wxy));
        for j=1:nch-1
            for k=j+1:nch
                PLI(:,:,j,k) = phase_lag_index(imag(Wxy(:,:,j,k)), ...
                    scale, nsig, dt, 'wPLI', wPLI);
                PLI(:,:,k,j) = PLI(:,:,j,k);
            end
        end
        PLIthreshd = PLI>sig;
        [Pcoh, Pinc] = pcc ( f, W, PLIthreshd, coi );
        P_coh_on(:,1:Nch(i),i) = squeeze(nanmean(Pcoh,3));
        P_inc_on(:,1:Nch(i),i) = squeeze(nanmean(Pinc,3));


    end

    P_tot_off = P_inc_off + P_coh_off;
    P_tot_on  = P_inc_on  + P_coh_on;
    
    clear x a W Wxy C coi1 Cs W1s Wxys coi1s



    %% Save Variables
    save(['LFP_pow_wPLI_p01_' tag '_' task '.mat'], ...
        'P_tot_on', 'P_tot_off', 'STNact', 'Nch', 'filt', ...
        'P_coh_off', 'P_coh_on', 'LFPch', 'sig', 'wPLI', ...
        'task', 'Nex', 'f', 'P_inc_off', 'P_inc_on', 'Patient', ...
        'LFPact', 'w0', 'scale', 'nsig', 'ds', 'rigor')
    
    
%% Plot
% load xyz
% B(:,:,3)=squeeze(nanmean(P_vc_on./P_tot_on-P_vc_off./P_tot_off,2));
% B(:,:,1)=squeeze(nanmean(P_inc_on./P_tot_on-P_inc_off./P_tot_off,2));
% B(:,:,2)=squeeze(nanmean(P_coh_on./P_tot_on-P_coh_off./P_tot_off,2));
% M = [mean(B(:,:,1),2) mean(B(:,:,2),2) mean(B(:,:,3),2)];
% E = [std(B(:,:,1),[],2) std(B(:,:,2),[],2) std(B(:,:,3),[],2)];
% mseb(f,M',E'/sqrt(8))