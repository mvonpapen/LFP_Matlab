%% Compute coh, inc, vc for cheops output with EMG data

GLEcoh = cell(6,2);
GLEall = cell(6,2);

%% Load LFP data
for Npat = 2;
    ds = 10;
    load('data_fist_v3')
    type = 'Fist';
    load(['gra_' type mat2str(Npat) '_w12n3.mat'])

    %% Parameters
    Nf      = length(f);
    sig     = sig_coh_thresh(w0, nsig);
    phthres = 15;

    
    %% Declare channels
    LFPch = {'C', 'L', 'A', 'M', 'P'};
    EMGch = {'EDCre', 'FDLre', 'FDIre', 'EDCli', 'FDLli', 'FDIli'};
    iLFP = find(ismember(Ch, LFPch));
    iEMG = find(ismember(Ch, EMGch));
    clear LFPch EMGch


    %% OFF
    [~,Woff,coi_off] = procdata(LFP_OFF{Npat}, 'freq', f, 'w0', w0, ...
        'filter', [], 'art', ArtLFP_OFF{Npat});
    if ~isempty(EMG_OFF{Npat})
%         [~,WE,coi_off_E,PE] = procdata(EMG_OFF{Npat}, 'freq', f, 'w0', w0, ...
%             'filter', [0 60; 90 110], 'art', ArtEMG_OFF{Npat}, 'rect', 1);
        [~,WE,coi_off_E,PE] = procdata(EMG_OFF{Npat}, 'freq', f, 'w0', w0, ...
            'filter', [], 'art', ArtEMG_OFF{Npat}, 'rect', 0);
        coi_off = [coi_off coi_off_E]/nsig;
        Woff(:,:,iEMG) = WE;
%     else
%         coi_off = coi_off/nsig;
    end
    [C, Wxy] = wave_cohere ( Woff, scale, nsig );
    
    % Downsampling
    C   = C(:,1:ds:end,:,:);
    Wxy = Wxy(:,1:ds:end,:,:);
    coi_off = coi_off(1:ds:end,:);
    
    % Apply cone of interest
    for i=1:size(C,3)
        for j=1:size(C,4)
            C(:,:,i,j) = coi2nan(f, C(:,:,i,j), min([coi_off(:,i)'; coi_off(:,j)']));
        end
    end
        
    Ph = abs(angle(Wxy)/pi*180);
%     inc_off = C<=sig;
    coh_off = C>sig & Ph>phthres & Ph<180-phthres;
%     vc_off  = C>sig & (Ph<=phthres | Ph>=180-phthres);
    clear Wxy x C coi_off_E
    
     
    %% Determine average causality between LFP and LFP
    GLLcoh{Npat,1} = zeros(Nf, 2, length(iLFP), length(iLFP));
    GLLall{Npat,1} = zeros(Nf, 2, length(iLFP), length(iLFP));
    l = find(ismember(c_off(:,1),iLFP) & ismember(c_off(:,2),iLFP));
    tmp = c_off(l,:);
    for j=1:length(l)
        for i=1:2;
            GLLcoh{Npat,1}(:,i,j) = squeeze( nanmean( smoothts( ...
                (G_off(:,:,l(j),i).*coh_off(:,:,tmp(j,1),tmp(j,2)))','b',2)));
            GLLall{Npat,1}(:,i,j) = squeeze( nanmean( smoothts( ...
                G_off(:,:,l(j),i)','b',2) ) );
        end
    end    
    
    %% Determine average causality between LFP and EMG
    GLEcoh{Npat,1} = zeros(Nf, 2, length(iLFP), length(iEMG));
    GLEall{Npat,1} = zeros(Nf, 2, length(iLFP), length(iEMG));
    l = find(ismember(c_off(:,1),iLFP) & ismember(c_off(:,2),iEMG));
    n = 1;
    for j=iLFP
        for k=1:length(iEMG)
            for i=1:2;
                GLEcoh{Npat,1}(:,i,j,k) = squeeze( nanmean( smoothts( ...
                    (G_off(:,:,l(n),i).*coh_off(:,:,j,k+iLFP(end)))','b',2)));
                GLEall{Npat,1}(:,i,j,k) = squeeze( nanmean( smoothts( ...
                    G_off(:,:,l(n),i)','b',2) ) );
            end
            n = n+1;
        end
    end
    
    
    %% ON
    [~,Won,coi_on,PL] = procdata(LFP_ON{Npat}, 'freq', f, 'w0', w0, ...
        'filter', [], 'art', ArtLFP_ON{Npat});
    if ~isempty(EMG_ON{Npat})
%         [~,WE,coi_on_E] = procdata(EMG_ON{Npat}, 'freq', f, 'w0', w0, ...
%             'filter', [0 60; 90 110], 'art', ArtEMG_ON{Npat}, 'rect', 1);
        [~,WE,coi_on_E] = procdata(EMG_ON{Npat}, 'freq', f, 'w0', w0, ...
            'filter', [], 'art', ArtEMG_ON{Npat}, 'rect', 0);
        coi_on = [coi_on coi_on_E]/nsig;
        Won(:,:,iEMG) = WE;
%     else
%         coi_on = coi_on/nsig;
    end
    [C, Wxy] = wave_cohere ( Won, scale, nsig );
    
    % Downsampling
    C      = C(:,1:ds:end,:,:);
    Wxy    = Wxy(:,1:ds:end,:,:);
    coi_on = coi_on(1:ds:end,:);
    
    % Apply cone of interest
    for i=1:size(C,3)
        for j=1:size(C,4)
            C(:,:,i,j) = coi2nan(f, C(:,:,i,j), min([coi_on(:,i)'; coi_on(:,j)']));
        end
    end
    
    Ph = abs(angle(Wxy)/pi*180);
%     inc_on = C<=sig;
    coh_on = C>sig & Ph>phthres & Ph<180-phthres;
%     vc_on  = C>sig & (Ph<=phthres | Ph>=180-phthres);
    clear Wxy x C W WE Ph coi_on_E
    
     
    %% Determine average causality between LFP and LFP
    GLLcoh{Npat,2} = zeros(Nf, 2, length(iLFP), length(iLFP));
    GLLall{Npat,2} = zeros(Nf, 2, length(iLFP), length(iLFP));
    l = find(ismember(c_on(:,1),iLFP) & ismember(c_on(:,2),iLFP));
    tmp = c_on(l,:);
    for j=1:length(l)
        for i=1:2;
            GLLcoh{Npat,2}(:,i,j) = squeeze( nanmean( smoothts( ...
                (G_on(:,:,l(j),i).*coh_on(:,:,tmp(j,1),tmp(j,2)))','b',2)));
            GLLall{Npat,2}(:,i,j) = squeeze( nanmean( smoothts( ...
                G_on(:,:,l(j),i)','b',2) ) );
        end
    end
    
     
    %% Determine average causality between LFP and EMG
    GLEcoh{Npat,2} = zeros(Nf, 2, length(iLFP), length(iEMG));
    GLEall{Npat,2} = zeros(Nf, 2, length(iLFP), length(iEMG));
    l = find(ismember(c_on(:,1),iLFP) & ismember(c_on(:,2),iEMG));
    n = 1;
    for j=iLFP
        for k=1:length(iEMG)
            for i=1:2;
                GLEcoh{Npat,2}(:,i,j,k) = squeeze( nanmean( smoothts( ...
                    (G_on(:,:,l(n),i).*coh_on(:,:,j,k+iLFP(end)))','b',2)));
                GLEall{Npat,2}(:,i,j,k) = squeeze( nanmean( smoothts( ...
                    G_on(:,:,l(n),i)','b',2) ) );
            end
            n = n+1;
        end
    end
% 
% 
%     %% Apply cone of interest
%     [Nc,i] = size(c_off);
%     [Nc2,i2] = size(c_on);
%     if Nc~=Nc2 || i~=i2
%         error('Sizes of coh_off and coh_on do not fit!')
%     end
%     clear Nc2 i2
%     for i=1:Nc
%         for j=1:2
%             G_off(:,:,i,j) = coi2nan(f,G_off(:,:,i,j),min(coi_off(:,c_off(i,:)),[],2));
%             G_on(:,:,i,j) = coi2nan(f,G_on(:,:,i,j),min(coi_on(:,c_on(i,:)),[],2));
%         end
%     end

%     %% Plot Granger causality for all electrode combinations
%     back = 0; % use a background of NaN or 0 to fill the tmp-matrices
%     fig1=figure('Papertype', 'A4', ...
%         'Paperunits', 'normalized', 'PaperPosition', [0 0 1 1], ...
%         'PaperPositionmode', 'manual', 'Visible', 'off', ...
%         'PaperOrientation', 'portrait');
%     annotation(fig1,'textbox', [0.4 0.95 0.2 0.03], 'String', [Patient ' ' type], ...
%         'Fontsize', 14, 'HorizontalAlign', 'center');
%     for i=1:Nc
%         for j=1:2
%             nj = setxor([1 2], j);
%             subplot(Nc,2,(i-1)*2+j)
%             %% OFF
%             plot(f, nanmean(G_off(:,:,i,j),2), ...
%                 '-k', 'linew', 1)
%             hold all
%             if ismember(i, iLFP)
%                 tmp = G_off(:,:,i,j);
%                 tmp(~inc_off(:,:,c_off(i,1),c_off(i,2))) = back;
%                 plot(f, nanmean(tmp,2), '-b')
%                 tmp = G_off(:,:,i,j);
%                 tmp(~coh_off(:,:,c_off(i,1),c_off(i,2))) = back;
%                 plot(f, nanmean(tmp,2), '-r')
%                 tmp = G_off(:,:,i,j);
%                 tmp(~vc_off(:,:,c_off(i,1),c_off(i,2))) = back;
%                 plot(f, nanmean(tmp,2), '-g')
%             else
%                 tmp = G_off(:,:,i,j);
%                 tmp(~inc_off(:,:,c_off(i,1),c_off(i,2))) = back;
%                 plot(f, nanmean(tmp,2), '-b')
%                 tmp = G_off(:,:,i,j);
%                 tmp(~coh_off(:,:,c_off(i,1),c_off(i,2)) & ...
%                      ~vc_off(:,:,c_off(i,1),c_off(i,2))) = back;
%                 plot(f, nanmean(tmp,2), '-r')
%             end
%             if i==1 && j==2
%                 legend('Total OFF', 'inc OFF', 'coh OFF', 'vc OFF', ...
%                        'Total ON',  'inc ON',  'coh ON',  'vc ON')
%             end
%                 
%             %% ON
%             plot(f, nanmean(coi2nan(f,G_on(:,:,i,j),min(coi_on(:,c_on(i,:)),[],2)),2), ...
%                 '--k', 'linew', 1)
%             if ismember(i, iLFP)
%                 tmp = G_on(:,:,i,j);
%                 tmp(~inc_on(:,:,c_on(i,1),c_on(i,2))) = back;
%                 tmp = coi2nan(f,tmp,min(coi_on(:,c_on(i,:)),[],2));
%                 plot(f, nanmean(tmp,2), '--b')
%                 tmp = G_on(:,:,i,j);
%                 tmp(~coh_on(:,:,c_on(i,1),c_on(i,2))) = back;
%                 tmp = coi2nan(f,tmp,min(coi_on(:,c_on(i,:)),[],2));
%                 plot(f, nanmean(tmp,2), '--r')
%                 tmp = G_on(:,:,i,j);
%                 tmp(~vc_on(:,:,c_on(i,1),c_on(i,2))) = back;
%                 tmp = coi2nan(f,tmp,min(coi_on(:,c_on(i,:)),[],2));
%                 plot(f, nanmean(tmp,2), '--g')
%             else
%                 tmp = G_on(:,:,i,j);
%                 tmp(~inc_on(:,:,c_on(i,1),c_on(i,2))) = back;
%                 plot(f, nanmean(tmp,2), '-b')
%                 tmp = G_on(:,:,i,j);
%                 tmp(~coh_on(:,:,c_on(i,1),c_on(i,2)) & ...
%                      ~vc_on(:,:,c_on(i,1),c_on(i,2))) = back;
%                 plot(f, nanmean(tmp,2), '-r')
%             end
%             %% Annotations
%             xlabel('f [Hz]')
%             ylabel('Granger causality')
%             ylim([0 2])
%             title([Channel{Npat}{c_off(i,j)} '->' Channel{Npat}{c_off(i,nj)}])
%         end
%     end
%     clear tmp i j nj
%     saveas(fig1, ['caus_' Patient '_' type '.fig'])
%     close(fig1)
end
for i=2; for j=1:2; Gmeancoh(:,:,i,j) = nanmean(nanmean(GLEcoh{i,j},4),3); end; end
for i=2; for j=1:2; Gmeanall(:,:,i,j) = nanmean(nanmean(GLEall{i,j},4),3); end; end

OFFm = mean(GLEall{2,1}(:,:,:),3);
OFFs = std(GLEall{2,1}(:,:,:),[],3);
ONm = mean(GLEall{2,2}(:,:,:),3);
ONs = std(GLEall{2,2}(:,:,:),[],3);
mseb(f,OFFm',OFFs'), xlim([0 17])
figure, mseb(f,ONm',ONs'), xlim([0 17])