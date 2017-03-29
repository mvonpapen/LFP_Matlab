%%  Plot relative and absolute (percentual) ON-OFF difference of PSDs

clear sig_off sig_on Won Woff

task = 'hold';


% Take out patients w/o task
if strcmp(task, 'fist')
    pat1 = [1 2 4 5 7 8];
    pat2 = 1:6;
elseif strcmp(task, 'hold')
    pat1 = [1 2 3 4 5 7 8];
    pat2 = 1:7;
end

textlab = {'Inc', 'Coh', 'VC'};
tag     = 'v3_bip'; %denotes which file to load
usemean = false; %false = median
Bonnf   = 6; % 1: no Bonnferroni correction
difftype = 'relRest'; % 'rel'=perc.diff, 'relRest'=normed to restpower, 'abs'=abs.diff
patmean = false;
% axis limits
limy    = [-40 40];
limy_bip= [-40 40];

P = NaN(6,6,2);

for ind = 1:2
    
    switch ind
        case 1
            fig1    = figure('Papersize', [7 6], 'PaperPosition', [0.75 0.5 5.5 5], ...
               'PaperPositionmode', 'manual', 'Visible', 'off');
            pnum = 1;
            fb     =  [13 20];
            fr_str = '\beta_1';
%             fnam   = ['psd_diff_mean_' tag '_' task '2rest_beta1.eps'];
        case 2
            pnum = 2;
            fb     =  [20 30];
            fr_str = '\beta_2';
            fnam   = ['psd_diff_box_' tag '_' task '2rest_beta12.eps'];
        case 3
            fig1    = figure('Papersize', [6 5], 'PaperPosition', [0.75 0.5 4.5 4], ...
               'PaperPositionmode', 'manual', 'Visible', 'off');
            pnum = 1;
            fb     =  [30 40];
            fr_str = '\gamma_1';
%             fnam   = ['psd_diff_mean_' tag '_' task '2rest_gamma1.eps'];
        case 4
            pnum = 2;
            fb     =  [60 90];
            fr_str = '\gamma_2';
            fnam   = ['psd_diff_box_' tag '_' task '2rest_gamma12.eps'];
    end

    lg_pos = [0.14, 0.48, 0, 0];

           
           
    %%%%%%%% vvv Rest vvv %%%%%%%%%%%%%%%

    load(['LFP_pow_' tag '_rest.mat'], 'P_tot_off', 'P_tot_on', ...
        'P_coh_off', 'P_coh_on', 'P_vc_off', 'P_vc_on', 'f', ...
        'P_inc_off', 'P_inc_on', 'P_bip_off', 'P_bip_on')
    j  = find(f>=fb(1) & f<=fb(2) & ~(f>45 & f<55));


    % Specify patients and average over LFP channels
    if patmean
        totx_rest = squeeze(nanmean(P_tot_off(:,:,pat1),2));
        toty_rest = squeeze(nanmean(P_tot_on(:,:,pat1),2));
        xinc_rest = squeeze(nanmean(P_inc_off(:,:,pat1),2));
        yinc_rest = squeeze(nanmean(P_inc_on(:,:,pat1),2));
        xcoh_rest = squeeze(nanmean(P_coh_off(:,:,pat1),2));
        ycoh_rest = squeeze(nanmean(P_coh_on(:,:,pat1),2));
        xvc_rest  = squeeze(nanmean(P_vc_off(:,:,pat1),2));
        yvc_rest  = squeeze(nanmean(P_vc_on(:,:,pat1),2));
        xbip_rest  = squeeze(nanmean(P_bip_off(:,:,pat1),2));
        ybip_rest  = squeeze(nanmean(P_bip_on(:,:,pat1),2));
        xaco_rest = squeeze(nanmean(P_coh_off(:,:,pat1)+P_vc_off(:,:,pat1),2));
        yaco_rest = squeeze(nanmean(P_coh_on(:,:,pat1)+P_vc_on(:,:,pat1),2));
    else
        totx_rest = P_tot_off(:,:,pat1); totx_rest=totx_rest(:,:);
        toty_rest = P_tot_on(:,:,pat1); toty_rest=toty_rest(:,:);
        xinc_rest = P_inc_off(:,:,pat1); xinc_rest=xinc_rest(:,:);
        yinc_rest = P_inc_on(:,:,pat1); yinc_rest=yinc_rest(:,:);
        xcoh_rest = P_coh_off(:,:,pat1); xcoh_rest=xcoh_rest(:,:);
        ycoh_rest = P_coh_on(:,:,pat1); ycoh_rest=ycoh_rest(:,:);
        xvc_rest = P_vc_off(:,:,pat1); xvc_rest=xvc_rest(:,:);
        yvc_rest = P_vc_on(:,:,pat1); yvc_rest=yvc_rest(:,:);
        xbip_rest = P_bip_off(:,:,pat1); xbip_rest=xbip_rest(:,:);
        ybip_rest = P_bip_on(:,:,pat1); ybip_rest=ybip_rest(:,:);
        xaco_rest = P_coh_off(:,:,pat1)+P_vc_off(:,:,pat1); xaco_rest=xaco_rest(:,:);
        yaco_rest = P_coh_on(:,:,pat1)+P_vc_on(:,:,pat1); yaco_rest=yaco_rest(:,:);
    end

    

    %%%%%%%% vvv Task vvv %%%%%%%%%%%%%%%

    load(['LFP_pow_' tag '_' task '.mat'], 'P_tot_off', 'P_tot_on', ...
        'P_coh_off', 'P_coh_on', 'P_vc_off', 'P_vc_on', 'f', ...
        'P_inc_off', 'P_inc_on', 'P_bip_off', 'P_bip_on')


    % Specify patients and average over LFP channels
    if patmean
        totx_task = squeeze(nanmean(P_tot_off(:,:,pat2),2));
        toty_task = squeeze(nanmean(P_tot_on(:,:,pat2),2));
        xinc_task = squeeze(nanmean(P_inc_off(:,:,pat2),2));
        yinc_task = squeeze(nanmean(P_inc_on(:,:,pat2),2));
        xcoh_task = squeeze(nanmean(P_coh_off(:,:,pat2),2));
        ycoh_task = squeeze(nanmean(P_coh_on(:,:,pat2),2));
        xvc_task  = squeeze(nanmean(P_vc_off(:,:,pat2),2));
        yvc_task  = squeeze(nanmean(P_vc_on(:,:,pat2),2));
        xbip_task  = squeeze(nanmean(P_bip_off(:,:,pat2),2));
        ybip_task  = squeeze(nanmean(P_bip_on(:,:,pat2),2));
        xaco_task = squeeze(nanmean(P_coh_off(:,:,pat2)+P_vc_off(:,:,pat2),2));
        yaco_task = squeeze(nanmean(P_coh_on(:,:,pat2)+P_vc_on(:,:,pat2),2));
    else
        totx_task = P_tot_off(:,:,pat2); totx_task=totx_task(:,:);
        toty_task = P_tot_on(:,:,pat2); toty_task=toty_task(:,:);
        xinc_task = P_inc_off(:,:,pat2); xinc_task=xinc_task(:,:);
        yinc_task = P_inc_on(:,:,pat2); yinc_task=yinc_task(:,:);
        xcoh_task = P_coh_off(:,:,pat2); xcoh_task=xcoh_task(:,:);
        ycoh_task = P_coh_on(:,:,pat2); ycoh_task=ycoh_task(:,:);
        xvc_task = P_vc_off(:,:,pat2); xvc_task=xvc_task(:,:);
        yvc_task = P_vc_on(:,:,pat2); yvc_task=yvc_task(:,:);
        xbip_task = P_bip_off(:,:,pat2); xbip_task=xbip_task(:,:);
        ybip_task = P_bip_on(:,:,pat2); ybip_task=ybip_task(:,:);
        xaco_task = P_coh_off(:,:,pat2)+P_vc_off(:,:,pat2); xaco_task=xaco_task(:,:);
        yaco_task = P_coh_on(:,:,pat2)+P_vc_on(:,:,pat2); yaco_task=yaco_task(:,:);
    end



    %% Differences of Relative Power
    if strcmp(difftype, 'rel')
        % OFF: Fist - Rest
        dPinc_off = squeeze( nanmean(...
            xinc_task(j,:)./totx_task(j,:) - xinc_rest(j,:)./totx_rest(j,:) ) )*100;
        dPcoh_off = squeeze( nanmean(...
            xcoh_task(j,:)./totx_task(j,:) - xcoh_rest(j,:)./totx_rest(j,:) ) )*100;
        dPvc_off  = squeeze( nanmean(...
            xvc_task(j,:)./totx_task(j,:) - xvc_rest(j,:)./totx_rest(j,:) ) )*100;
        dPaco_off = squeeze( nanmean(...
            xaco_task(j,:)./totx_task(j,:) - xaco_rest(j,:)./totx_rest(j,:) ) )*100;
        % ON: Fist - Rest
        dPinc_on = squeeze( nanmean(...
            yinc_task(j,:)./toty_task(j,:) - yinc_rest(j,:)./toty_rest(j,:) ) )*100;
        dPcoh_on = squeeze( nanmean(...
            ycoh_task(j,:)./toty_task(j,:) - ycoh_rest(j,:)./toty_rest(j,:) ) )*100;
        dPvc_on  = squeeze( nanmean(...
            yvc_task(j,:)./toty_task(j,:) - yvc_rest(j,:)./toty_rest(j,:) ) )*100;
        dPaco_on = squeeze( nanmean(...
            yaco_task(j,:)./toty_task(j,:) - yaco_rest(j,:)./toty_rest(j,:) ) )*100;
    end


    %% Differences of Power Relative to Total Rest
    if strcmp(difftype, 'relRest')
        % OFF: Fist - Rest
        dPinc_off = squeeze( nanmean( (xinc_task(j,:) - xinc_rest(j,:)) ./ ...
            totx_rest(j,:)) )*100;
        dPcoh_off = squeeze( nanmean( (xcoh_task(j,:) - xcoh_rest(j,:)) ./ ...
            totx_rest(j,:)) )*100;
        dPvc_off  = squeeze( nanmean( (xvc_task(j,:)  - xvc_rest(j,:) ) ./ ...
            totx_rest(j,:)) )*100;
        dPbip_off = squeeze( nanmean( (xbip_task(j,:) - xbip_rest(j,:)) ./ ...
            xbip_rest(j,:)) )*100;
        dPaco_off = squeeze( nanmean( (xaco_task(j,:) - xaco_rest(j,:)) ./ ...
            totx_rest(j,:)) )*100;
        % ON: Fist - Rest
        dPinc_on = squeeze( nanmean( (yinc_task(j,:) - yinc_rest(j,:)) ./ ...
            toty_rest(j,:)) )*100;
        dPcoh_on = squeeze( nanmean( (ycoh_task(j,:) - ycoh_rest(j,:)) ./ ...
            toty_rest(j,:)) )*100;
        dPvc_on  = squeeze( nanmean( (yvc_task(j,:)  - yvc_rest(j,:) ) ./ ...
            toty_rest(j,:)) )*100;
        dPbip_on = squeeze( nanmean( (ybip_task(j,:) - ybip_rest(j,:)) ./ ...
            ybip_rest(j,:)) )*100;
        dPaco_on = squeeze( nanmean( (yaco_task(j,:) - yaco_rest(j,:)) ./ ...
            toty_rest(j,:)) )*100;
    end
                        
                        

    %% Differences of Total Power
    if strcmp(difftype, 'abs')
        % axis limits
        limy    = [-20 20]*3e-9;
        % OFF: Fist - Rest
        dPinc_off = squeeze( nanmean( xinc_task(j,:) - xinc_rest(j,:) ));
        dPcoh_off = squeeze( nanmean( xcoh_task(j,:) - xcoh_rest(j,:) ));
        dPvc_off  = squeeze( nanmean( xvc_task(j,:)  - xvc_rest(j,:)  ));
        dPaco_off = squeeze( nanmean( xaco_task(j,:) - xaco_rest(j,:) ));
        dPbip_off = squeeze( nanmean( xbip_task(j,:) - xbip_rest(j,:) ));
        % ON: Fist - Rest
        dPinc_on = squeeze( nanmean( yinc_task(j,:) - yinc_rest(j,:) ));
        dPcoh_on = squeeze( nanmean( ycoh_task(j,:) - ycoh_rest(j,:) ));
        dPvc_on  = squeeze( nanmean( yvc_task(j,:)  - yvc_rest(j,:)  ));
        dPaco_on = squeeze( nanmean( yaco_task(j,:) - yaco_rest(j,:) ));
        dPbip_on = squeeze( nanmean( ybip_task(j,:) - ybip_rest(j,:) ));
    end
                        
                        
    clear xinc yinc xcoh ycoh xvc yvc totx toty



    % Wilcoxon sign-rank test
    [sig_off(1,ind),~,tmp] = signrank ( dPinc_off  );
    Woff(1,ind) = tmp.signedrank;
    [sig_off(2,ind),~,tmp]  = signrank ( dPcoh_off );
    Woff(2,ind) = tmp.signedrank;
    [sig_off(3,ind),~,tmp]  = signrank ( dPvc_off );
    Woff(3,ind) = tmp.signedrank;

    [sig_on(1,ind),~,tmp] = signrank ( dPinc_on  );
    Won(1,ind) = tmp.signedrank;
    [sig_on(2,ind),~,tmp]  = signrank ( dPcoh_on );
    Won(2,ind) = tmp.signedrank;
    [sig_on(3,ind),~,tmp]  = signrank ( dPvc_on );
    Won(3,ind) = tmp.signedrank;
    
    [sig_off(4,ind),~,tmp]  = signrank ( dPbip_off );
    Won(4,ind) = tmp.signedrank;
    [sig_on(4,ind),~,tmp]  = signrank ( dPbip_on );
    Won(4,ind) = tmp.signedrank;
        


    sig_off
    sig_on

    
%     %% Variables for superbar
%     N = sum(~isnan([dPinc_off; dPcoh_off; dPvc_off]),2);
%     M(1,:,pnum) = nanmedian([dPinc_off; dPcoh_off; dPvc_off],2);
%     M(2,:,pnum) = nanmedian([dPinc_on; dPcoh_on; dPvc_on],2);
%     E(1,:,pnum) = nanstd([dPinc_off; dPcoh_off; dPvc_off],[],2)./sqrt(N);
%     E(2,:,pnum) = nanstd([dPinc_on; dPcoh_on; dPvc_on],[],2)./sqrt(N);
%     P(1,2,pnum) = ranksum(dPinc_off, dPinc_on); P(2,1,pnum)=P(1,2,pnum);
%     P(3,4,pnum) = ranksum(dPcoh_off, dPcoh_on); P(4,3,pnum)=P(3,4,pnum);
%     P(5,6,pnum) = ranksum(dPvc_off, dPvc_on); P(6,5,pnum)=P(5,6,pnum);
% %     C = [.8 .2 .2; .2 .2 .8; .2 .8 .2];
% %     superbar(M(:,:,i), 'E', E(:,:,i), 'P', P(:,:,i), 'BarFaceColor', C);
    
    %% Group results

    % PCC
    subplot(2,1,pnum)
    clear XX
    xx = 1:2;
    XX(:,1,1) = dPinc_off; XX(:,1,2) = dPcoh_off; XX(:,1,3) = dPvc_off;
    XX(:,2,1) = dPinc_on; XX(:,2,2) = dPcoh_on; XX(:,2,3) = dPvc_on;
    boxPlot(xx, XX, 'showScatter', true, 'groupwidth', 0.7, 'xSpacing', 'x', ...
        'symbolMarker', '+', 'useMean', usemean, ...
        'showlegend', true, 'grouplabels', {'Incoherent', 'Coherent', 'Vol.-cond.'}, ...
        'boxcolor', {[1 0.5 0.5], [0.5 0.5 1], [0.5 1 0.5]}, ...
        'LineStyle', '--', 'Linewidth', 0.7, 'MedianColor', 'k', ...
        'ScatterMarker', '.', 'ScatterColor', 'k', ...
        'symbolcolor', {[0.7 0 0], [0 0 0.7], [0 0.7 0]}, 'outlierSize', 13);
    box on; grid on
    hold all
    set(gca, 'XTick', [1 2], 'XTickLabel', {'OFF', 'ON'})
    ylim(limy), grid on
    ylabel('$\Delta P$', 'Interpreter', 'latex')
    title(['PCC, ' fr_str])
    if pnum==1
        annotation(fig1,'textbox', [0.01 0.93 0.12 0.07], 'String', task, ...
            'Fontsize', 14, 'Linestyle', 'none', 'Units', 'Normalized', ...
            'Fontweight', 'bold', 'Backgroundcolor', 'none');
    end

    for ii=1:length(sig_off)
        if sig_off(ii,ind)*Bonnf<0.05;        
            plot(0.6+ii*0.2+[-0.03 0.03],0.93*limy([2 2]),'*k');
        elseif sig_off(ii,ind)<0.05;
            plot(0.6+ii*0.2,0.93*limy(2),'*k');
        end
    end
    for ii=1:length(sig_on)-1
        if sig_on(ii,ind)*Bonnf<0.05;        
            plot(1.6+ii*0.2+[-0.03 0.03],0.93*limy([2 2]),'*k');
        elseif sig_on(ii,ind)<0.05;
            plot(1.6+ii*0.2,0.93*limy(2),'*k');
        end
    end
    
    % Bipolar
    clear YY
    figure
    xx = 1:2;
    YY(:,1) = dPbip_off; YY(:,2) = dPbip_on;
    boxPlot(xx, YY, 'showScatter', true, 'groupwidth', 0.7, 'xSpacing', 'x', ...
        'symbolMarker', '+', 'useMean', usemean, ...
        'LineStyle', '--', 'Linewidth', 0.7, 'MedianColor', 'k', ...
        'ScatterMarker', '.', 'ScatterColor', 'k', ...
        'outlierSize', 13);
    box on; grid on
    hold all
    set(gca, 'XTick', [1 2], 'XTickLabel', {'OFF', 'ON'})
    ylim(limy_bip), grid on
    ylabel('$\Delta P$', 'Interpreter', 'latex')
    title(['Bipolar, ' fr_str])
    if pnum==1
        annotation(fig1,'textbox', [0.01 0.93 0.12 0.07], 'String', task, ...
            'Fontsize', 14, 'Linestyle', 'none', 'Units', 'Normalized', ...
            'Fontweight', 'bold', 'Backgroundcolor', 'none');
    end

    if sig_off(4,ind)*Bonnf<0.05;        
        plot(0.8+[-0.03 0.03],0.93*limy([2 2]),'*k');
    elseif sig_off(4,ind)<0.05;
        plot(0.8,0.93*limy(2),'*k');
    end
    if sig_on(4,ind)*Bonnf<0.05;        
        plot(1.8+[-0.03 0.03],0.93*limy([2 2]),'*k');
    elseif sig_on(4,ind)<0.05;
        plot(1.8,0.93*limy(2),'*k');
    end


    set(fig1, 'Visi', 'on')
%     if pnum==2
%         print(fig1, fnam, '-depsc')
%         close(fig1)
%     end
    
end
