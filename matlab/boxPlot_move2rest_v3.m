%%  Plot relative and absolute (percentual) ON-OFF difference of PSDs


TASKS = {'hold', 'fist'};
fb = [20 30];
fr_str = '\beta_2';
cutoff  = 0;
textlab = {'Inc', 'Coh', 'VC', 'Bip'};
STNact  = 1;
tag     = 'PCC_v1_filt'; %denotes which file to load
usemean = false; %false = median
patmean = false;
Bonnf   = 6; % 1: no Bonnferroni correction

fig1    = figure('Papersize', [7 6], 'PaperPosition', [0.75 0.5 5.5 5], ...
   'PaperPositionmode', 'manual', 'Visible', 'off');
fnam   = ['psd_diff_box_' tag '_move2rest_beta2.eps'];




for ind=1:2
    
    task = TASKS{ind};

    % Take out patients w/o task
    if strcmp(task, 'fist')
        pat1 = [1 2 4 5 7 8];
        pat2 = [1:6];
    elseif strcmp(task, 'hold')
        pat1 = [1 2 3 4 5 7 8];
        pat2 = [1:7];
    end

    % % Log-averages?
    logavg = 0;
    if logavg==1
        limy   = [-1 1];
        limy2  = [-1 1];
        ytic2  = [-1 -0.5 0 0.5 1];
    else
        limy    = [-40 40];
        limy2   = [-40 40];
        ytic2   = [-40 -20 0 20 40];
    end


    lg_pos = [0.14, 0.48, 0, 0];

           
           
    %%%%%%%% vvv Rest vvv %%%%%%%%%%%%%%%

    load(['LFP_pow_' tag '_rest.mat'], ...
        'P_coh_off', 'P_coh_on', 'P_vc_off', 'P_vc_on', 'f', 'LFPact', ...
        'P_inc_off', 'P_inc_on', 'P_bip_off', 'P_bip_on')
    j  = find(f>=fb(1) & f<=fb(2) & ~(f>40 & f<60));


    val = cutoff;
    P_tot_off = P_coh_off + P_inc_off + P_vc_off;
    P_tot_on  = P_coh_on + P_inc_on + P_vc_on;
    P_inc_off(P_inc_off<=cutoff) = val;
    P_inc_on(P_inc_on<=cutoff)   = val;
    P_coh_off(P_coh_off<=cutoff) = val;
    P_coh_on(P_coh_on<=cutoff)   = val;
    P_vc_off(P_vc_off<=cutoff)   = val;
    P_vc_on(P_vc_on<=cutoff)     = val;


    % Only use channels with STN activity
    for i=1:length(LFPact)
        P_tot_off(:,LFPact{i}<STNact,i) = NaN;
        P_inc_off(:,LFPact{i}<STNact,i) = NaN;
        P_coh_off(:,LFPact{i}<STNact,i) = NaN;
        P_vc_off(:,LFPact{i}<STNact,i) = NaN;
    end
   
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

    % Logarithmic averages?
    if logavg==1
        totx_rest = log10(P_tot_off(:,:));
        toty_rest = log10(P_tot_on(:,:));
        xinc_rest = log10(P_inc_off(:,:));
        yinc_rest = log10(P_inc_on(:,:));
        xcoh_rest = log10(P_coh_off(:,:));
        ycoh_rest = log10(P_coh_on(:,:));
        xvc_rest  = log10(P_vc_off(:,:));
        yvc_rest  = log10(P_vc_on(:,:));
        xbip_rest  = log10(P_bip_off(:,:));
        ybip_rest  = log10(P_bip_on(:,:));
    end

    %%%%%%%% vvv Fist vvv %%%%%%%%%%%%%%%

    load(['LFP_pow_' tag '_' task '.mat'], 'P_tot_off', 'P_tot_on', ...
        'P_coh_off', 'P_coh_on', 'P_vc_off', 'P_vc_on', 'f', 'LFPact', ...
        'P_inc_off', 'P_inc_on', 'P_bip_off', 'P_bip_on')


    val = cutoff;
    P_inc_off(P_inc_off<=cutoff) = val;
    P_inc_on(P_inc_on<=cutoff)   = val;
    P_coh_off(P_coh_off<=cutoff) = val;
    P_coh_on(P_coh_on<=cutoff)   = val;
    P_vc_off(P_vc_off<=cutoff)   = val;
    P_vc_on(P_vc_on<=cutoff)     = val;


    
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



    %% Differences (normalized by total power)
    % OFF
    dPtot_off = squeeze( nanmean( (totx_task(j,:)-totx_rest(j,:) ) ...
                            ./ totx_rest(j,:) ) )*100;
    dPinc_off = squeeze( nanmean( (xinc_task(j,:)-xinc_rest(j,:) ) ...
                            ./ totx_rest(j,:) ) )*100;
    dPcoh_off = squeeze( nanmean( (xcoh_task(j,:)-xcoh_rest(j,:) ) ...
                            ./ totx_rest(j,:) ) )*100;
    dPvc_off  = squeeze( nanmean( (xvc_task(j,:) -xvc_rest(j,:) ) ...
                            ./ totx_rest(j,:) ) )*100;
    dPbip_off = squeeze( nanmean( (xbip_task(j,:) -xbip_rest(j,:) ) ...
                            ./ xbip_rest(j,:) ) )*100;
    % ON
    dPtot_on = squeeze( nanmean( (toty_task(j,:)-toty_rest(j,:) ) ...
                            ./ toty_rest(j,:) ) )*100;
    dPinc_on = squeeze( nanmean( (yinc_task(j,:)-yinc_rest(j,:) ) ...
                            ./ toty_rest(j,:) ) )*100;
    dPcoh_on = squeeze( nanmean( (ycoh_task(j,:)-ycoh_rest(j,:) ) ...
                            ./ toty_rest(j,:) ) )*100;
    dPvc_on  = squeeze( nanmean( (yvc_task(j,:) -yvc_rest(j,:) ) ...
                            ./ toty_rest(j,:) ) )*100;
    dPbip_on = squeeze( nanmean( (ybip_task(j,:) -ybip_rest(j,:) ) ...
                            ./ ybip_rest(j,:) ) )*100;
                        

    clear xinc yinc xcoh ycoh xvc yvc totx toty
    dPinc_off(dPinc_off==Inf) = 100;
    dPinc_on(dPinc_on==Inf) = 100;
    dPcoh_off(dPcoh_off==Inf) = 100;
    dPcoh_on(dPcoh_on==Inf) = 100;



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



    %% Group results

    % PCC
    subplot(2,3,[1 2]+(ind-1)*3)
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
    ylabel(['$\hbox{\fontsize{16}{20}\selectfont\(\mathbf{\frac{P^\mathrm{' task '}_i-P^\mathrm{rest}_i}{P^\mathrm{rest}_\mathrm{tot}}}\)}\,[\%]$'], ...
        'Interpreter', 'latex', 'Fontsize', 16)
    title(['Changes in $' fr_str '$ induced by ' task], 'Interpreter', 'Latex')

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
    subplot(2,3,ind*3)
    xx = 1:2;
    clear XX
    XX(:,1,1) = dPbip_off;
    XX(:,2,1) = dPbip_on;
    boxPlot(xx, XX, 'showScatter', true, 'groupwidth', 0.5, 'xSpacing', 'x', ...
        'symbolMarker', '+', 'useMean', usemean, ...
        'LineStyle', '--', 'Linewidth', 0.7, 'MedianColor', 'k', ...
        'ScatterMarker', '.', 'outlierSize', 13, 'symbolColor', 'k');
    box on; grid on
    hold all
    set(gca, 'XTick', [1 2], 'XTickLabel', {'OFF', 'ON'}, 'YTick', ytic2)
    ylim(limy2), grid on
    title(['Bipolar, ' task], 'Interpreter', 'Latex')
    ylabel(['$\hbox{\fontsize{16}{20}\selectfont\(\mathbf{\frac{P^\mathrm{' task '}-P^\mathrm{rest}}{P^\mathrm{rest}}}\)}\,[\%]$'], ...
        'Interpreter', 'latex', 'Fontsize', 14)

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
%     if ind==2
%         print(fig1, fnam, '-depsc')
%         close(fig1)
%     end
    
end