%% Compute coh, inc, vc for cheops output

%% Load LFP data
Npat = 2;
% for Npat = 1:8;
    type = 'rest'; %'Halte', 'Ruhe';
    load(['data_' type '_v3.mat']);
    load(['granger_' type mat2str(Npat) '_dec50.mat'])

    %% Parameters
    Nf      = length(f);
    sig     = 0.43;
    phthres = 15;



    %% OFF
    [~,W1,coi_off] = procdata(LFP_OFF{Npat}, 'freq', f, 'w0', w0, ...
        'filter', [], 'art', ArtLFP_OFF{Npat});
    coi_off = coi_off(dec:dec:end,:); %/t_smo;
    [C, Wxy] = wave_cohere ( W1, scale, t_smo );
    Ph = abs(angle(Wxy)/pi*180);

    inc_off = C<=sig;
    coh_off = C>sig & Ph>phthres & Ph<180-phthres;
    vc_off  = C>sig & (Ph<=phthres | Ph>=180-phthres);
    clear Wxy x C

    %% ON
    [~,W1,coi_on] = procdata(LFP_ON{Npat}, 'freq', f, 'w0', w0, ...
        'filter', [], 'art', ArtLFP_ON{Npat});
    coi_on = coi_on(dec:dec:end,:); %/t_smo;
    [C, Wxy] = wave_cohere ( W1, scale, t_smo );
    Ph = abs(angle(Wxy)/pi*180);

    inc_on = C<=sig;
    coh_on = C>sig & Ph>phthres & Ph<180-phthres;
    vc_on  = C>sig & (Ph<=phthres | Ph>=180-phthres);
    clear Wxy x C


    %% Apply cone of interest
    [Nc,i] = size(c_off);
    [Nc2,i2] = size(c_on);
    if Nc~=Nc2 || i~=i2
        error('Sizes of coh_off and coh_on do not fit!')
    end
    clear Nc2 i2        
    for i=1:Nc
        for j=1:2
            G_off(:,:,i,j) = coi2nan(f,G_off(:,:,i,j),min(coi_off(:,c_off(i,:)),[],2));
            G_on(:,:,i,j) = coi2nan(f,G_on(:,:,i,j),min(coi_on(:,c_on(i,:)),[],2));
        end
    end

    %% Plot Granger causality for all electrode combinations
    back = 0; % use a background of NaN or 0 to fill the tmp-matrices
    fig1=figure('Papertype', 'A4', ...
        'Paperunits', 'normalized', 'PaperPosition', [0 0 1 1], ...
        'PaperPositionmode', 'manual', 'Visible', 'off', ...
        'PaperOrientation', 'portrait');
    annotation(fig1,'textbox', [0.4 0.95 0.2 0.03], 'String', [Patient ' ' type], ...
        'Fontsize', 14, 'HorizontalAlign', 'center');
    for i=1:Nc
        for j=1:2
            nj = setxor([1 2], j);
            subplot(Nc,2,(i-1)*2+j)
            %% OFF
            tmp = G_off(:,:,i,j);
            tmp(~inc_off(:,:,c_off(i,1),c_off(i,2))) = back;
            plot(f, nanmean(G_off(:,:,i,j),2), ...
                '-k', 'linew', 1)
            hold all
            plot(f, nanmean(tmp,2), '-b')
            tmp = G_off(:,:,i,j);
            tmp(~coh_off(:,:,c_off(i,1),c_off(i,2))) = back;
            plot(f, nanmean(tmp,2), '-r')
            tmp = G_off(:,:,i,j);
            tmp(~vc_off(:,:,c_off(i,1),c_off(i,2))) = back;
            plot(f, nanmean(tmp,2), '-g')
            %% ON
            tmp = G_on(:,:,i,j);
            tmp(~inc_on(:,:,c_on(i,1),c_on(i,2))) = back;
            tmp = coi2nan(f,tmp,min(coi_on(:,c_on(i,:)),[],2));
            plot(f, nanmean(coi2nan(f,G_on(:,:,i,j),min(coi_on(:,c_on(i,:)),[],2)),2), ...
                '--k', 'linew', 1)
            plot(f, nanmean(tmp,2), '--b')
            tmp = G_on(:,:,i,j);
            tmp(~coh_on(:,:,c_on(i,1),c_on(i,2))) = back;
            tmp = coi2nan(f,tmp,min(coi_on(:,c_on(i,:)),[],2));
            plot(f, nanmean(tmp,2), '--r')
            tmp = G_on(:,:,i,j);
            tmp(~vc_on(:,:,c_on(i,1),c_on(i,2))) = back;
            tmp = coi2nan(f,tmp,min(coi_on(:,c_on(i,:)),[],2));
            plot(f, nanmean(tmp,2), '--g')
            %% Annotations
            xlabel('f [Hz]')
            ylabel('Granger causality')
            ylim([0 2])
            title(['LFP ' mat2str(c_off(i,j)) '->' mat2str(c_off(i,nj))])
        end
    end
    legend('Total OFF', 'inc OFF', 'coh OFF', 'vc OFF', ...
           'Total ON',  'inc ON',  'coh ON',  'vc ON')
    saveas(fig1, ['caus_' Patient '_' type '.pdf'])
    close(fig1)
% end