%%  Plot relative and absolute (percentual) ON-OFF difference of PSDs


task = 'fist';


% Take out patients w/o task
if strcmp(task, 'fist')
    pat1 = [1 2 4 5 7 8];
    pat2 = 1:6;
elseif strcmp(task, 'hold')
    pat1 = [1 2 3 4 5 7 8];
    pat2 = 1:7;
end

cutoff  = 1e-10;
textlab = {'Inc', 'Coh', 'VC', 'Bip'};
STNact  = 1;
tag     = 'v3_bip'; %denotes which file to load
usemean = false; %false = median
Bonnf   = 8; % 1: no Bonnferroni correction

% % Log-averages?
logavg = 0;
if logavg==1
    limy   = [-1 1];
    limy2  = [-1 1];
    ytic2  = [-1 -0.5 0 0.5 1];
else
    limy    = [-20 20];
    limy2   = [-60 60];
    ytic2   = [-60 -30 0 30 60];
end


for ind = 1:2%:4
    
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

    % Specify patients
    P_tot_off = P_tot_off(:,:,pat1);
    P_tot_on = P_tot_on(:,:,pat1);
    P_inc_off = P_inc_off(:,:,pat1);
    P_inc_on = P_inc_on(:,:,pat1);
    P_coh_off = P_coh_off(:,:,pat1);
    P_coh_on = P_coh_on(:,:,pat1);
    P_vc_off = P_vc_off(:,:,pat1);
    P_vc_on = P_vc_on(:,:,pat1);
    P_bip_off = P_bip_off(:,:,pat1);
    P_bip_on = P_bip_on(:,:,pat1);

    % Normalized absolute differences of PSDs
    totx_rest = P_tot_off(:,:);
    toty_rest = P_tot_on(:,:);
    xinc_rest = P_inc_off(:,:);
    yinc_rest = P_inc_on(:,:);
    xcoh_rest = P_coh_off(:,:);
    ycoh_rest = P_coh_on(:,:);
    xvc_rest  = P_vc_off(:,:);
    yvc_rest  = P_vc_on(:,:);
    xbip_rest  = P_bip_off(:,:);
    ybip_rest  = P_bip_on(:,:);

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


    % Specify patients
    P_tot_off = P_tot_off(:,:,pat2);
    P_tot_on = P_tot_on(:,:,pat2);
    P_inc_off = P_inc_off(:,:,pat2);
    P_inc_on = P_inc_on(:,:,pat2);
    P_coh_off = P_coh_off(:,:,pat2);
    P_coh_on = P_coh_on(:,:,pat2);
    P_vc_off = P_vc_off(:,:,pat2);
    P_vc_on = P_vc_on(:,:,pat2);
    P_bip_off = P_bip_off(:,:,pat2);
    P_bip_on = P_bip_on(:,:,pat2);


    % Normalized absolute differences of PSDs
    totx_fist = P_tot_off(:,:);
    toty_fist = P_tot_on(:,:);
    xinc_fist = P_inc_off(:,:);
    yinc_fist = P_inc_on(:,:);
    xcoh_fist = P_coh_off(:,:);
    ycoh_fist = P_coh_on(:,:);
    xvc_fist  = P_vc_off(:,:);
    yvc_fist  = P_vc_on(:,:);
    xbip_fist  = P_bip_off(:,:);
    ybip_fist  = P_bip_on(:,:);

    % Logarithmic averages?
    if logavg==1
        totx_fist = log10(P_tot_off(:,:));
        toty_fist = log10(P_tot_on(:,:));
        xinc_fist = log10(P_inc_off(:,:));
        yinc_fist = log10(P_inc_on(:,:));
        xcoh_fist = log10(P_coh_off(:,:));
        ycoh_fist = log10(P_coh_on(:,:));
        xvc_fist  = log10(P_vc_off(:,:));
        yvc_fist  = log10(P_vc_on(:,:));
        xbip_fist  = log10(P_bip_off(:,:));
        ybip_fist  = log10(P_bip_on(:,:));
    end



    %% Differences (normalized by total power)
    % OFF
    dPabs_tot_off = squeeze( nanmean( totx_fist(j,:)-totx_rest(j,:) ) ...
                            ./ nanmean( totx_rest(j,:) ) )*100;
    dPabs_inc_off = squeeze( nanmean( xinc_fist(j,:)-xinc_rest(j,:) ) ...
                            ./ nanmean(totx_rest(j,:) ) )*100;
    dPabs_coh_off = squeeze( nanmean( xcoh_fist(j,:)-xcoh_rest(j,:) ) ...
                            ./ nanmean(totx_rest(j,:) ) )*100;
    dPabs_vc_off  = squeeze( nanmean( xvc_fist(j,:) -xvc_rest(j,:) ) ...
                            ./ nanmean(totx_rest(j,:) ) )*100;
    dPabs_bip_off = squeeze( nanmean( xbip_fist(j,:) -xbip_rest(j,:) ) ...
                            ./ nanmean(xbip_rest(j,:) ) )*100;
    % ON
    dPabs_tot_on = squeeze( nanmean( toty_fist(j,:)-toty_rest(j,:) ) ...
                            ./ nanmean(toty_rest(j,:) ) )*100;
    dPabs_inc_on = squeeze( nanmean( yinc_fist(j,:)-yinc_rest(j,:) ) ...
                            ./ nanmean(toty_rest(j,:) ) )*100;
    dPabs_coh_on = squeeze( nanmean( ycoh_fist(j,:)-ycoh_rest(j,:) ) ...
                            ./ nanmean(toty_rest(j,:) ) )*100;
    dPabs_vc_on  = squeeze( nanmean( yvc_fist(j,:) -yvc_rest(j,:) ) ...
                            ./ nanmean(toty_rest(j,:) ) )*100;
    dPabs_bip_on = squeeze( nanmean( ybip_fist(j,:) -ybip_rest(j,:) ) ...
                            ./ nanmean(ybip_rest(j,:) ) )*100;
    %% Log-Differences (not normalized)                    
    if logavg==1
        % OFF
        dPabs_tot_off = squeeze( nanmean( totx_fist(j,:)-totx_rest(j,:) ) );
        dPabs_inc_off = squeeze( nanmean( xinc_fist(j,:)-xinc_rest(j,:) ) );
        dPabs_coh_off = squeeze( nanmean( xcoh_fist(j,:)-xcoh_rest(j,:) ) );
        dPabs_vc_off  = squeeze( nanmean( xvc_fist(j,:) -xvc_rest(j,:) ) );
        dPabs_bip_off = squeeze( nanmean( xbip_fist(j,:)-xbip_rest(j,:) ) );
        % ON
        dPabs_tot_on = squeeze( nanmean( toty_fist(j,:)-toty_rest(j,:) ) );
        dPabs_inc_on = squeeze( nanmean( yinc_fist(j,:)-yinc_rest(j,:) ) );
        dPabs_coh_on = squeeze( nanmean( ycoh_fist(j,:)-ycoh_rest(j,:) ) );
        dPabs_vc_on  = squeeze( nanmean( yvc_fist(j,:) -yvc_rest(j,:) ) );
        dPabs_bip_on = squeeze( nanmean( ybip_fist(j,:)-ybip_rest(j,:) ) );
    end          

    clear xinc yinc xcoh ycoh xvc yvc totx toty
    dPabs_inc_off(dPabs_inc_off==Inf) = 100;
    dPabs_inc_on(dPabs_inc_on==Inf) = 100;
    dPabs_coh_off(dPabs_coh_off==Inf) = 100;
    dPabs_coh_on(dPabs_coh_on==Inf) = 100;



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
    subplot(2,3,[1 2]+(pnum-1)*3)
    clear XX
    xx = 1:2;
    XX(:,1,1) = dPabs_inc_off; XX(:,1,2) = dPabs_coh_off; XX(:,1,3) = dPabs_vc_off;
    XX(:,2,1) = dPabs_inc_on; XX(:,2,2) = dPabs_coh_on; XX(:,2,3) = dPabs_vc_on;
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
    ylabel(['$\hbox{\fontsize{16}{20}\selectfont\(\mathbf{\frac{P^\mathrm{' task '}_i-P^\mathrm{rest}_i}{P^\mathrm{rest}_\mathrm{tot}}\,[\%]}\)}$'], ...
        'Interpreter', 'latex', 'Fontsize', 16)
    title(['PCC, ' fr_str])
    if pnum==1
        annotation(fig1,'textbox', [0.01 0.93 0.12 0.07], 'String', task, ...
            'Fontsize', 14, 'Linestyle', 'none', 'Units', 'Normalized', ...
            'Fontweight', 'bold', 'Backgroundcolor', 'none');
    end

    for ii=1:length(sig_off)-1
        if sig_off(ii)<0.01;        
            plot(0.6+ii*0.2+[-0.03 0.03],0.93*limy([2 2]),'*k');
        elseif sig_off(ii)<0.05;
            plot(0.6+ii*0.2,0.93*limy(2),'*k');
        end
    end
    for ii=1:length(sig_on)-1
        if sig_on(ii)<0.01;        
            plot(1.6+ii*0.2+[-0.03 0.03],0.93*limy([2 2]),'*k');
        elseif sig_on(ii)<0.05;
            plot(1.6+ii*0.2,0.93*limy(2),'*k');
        end
    end
    % Bipolar
    subplot(2,3,pnum*3)
    xx = 1:2;
    clear XX
    XX(:,1,1) = dPabs_bip_off;
    XX(:,2,1) = dPabs_bip_on;
    boxPlot(xx, XX, 'showScatter', true, 'groupwidth', 0.5, 'xSpacing', 'x', ...
        'symbolMarker', '+', 'useMean', usemean, ...
        'LineStyle', '--', 'Linewidth', 0.7, 'MedianColor', 'k', ...
        'ScatterMarker', '.', 'outlierSize', 13, 'symbolColor', 'k');
    box on; grid on
    hold all
    set(gca, 'XTick', [1 2], 'XTickLabel', {'OFF', 'ON'}, 'YTick', ytic2)
    ylim(limy2), grid on
    title(['Bipolar, ' fr_str])

%     if sig_off(4)<0.01;        
%         plot(1+[-0.05 0.05],0.93*limy2([2 2]),'*k');
%     elseif sig_off(4)<0.05;
%         plot(1,0.93*limy2(2),'*k');
%     end
%     if sig_on(4)<0.01;        
%         plot(2+[-0.05 0.05],0.93*limy2([2 2]),'*k');
%     elseif sig_on(4)<0.05;
%         plot(2,0.93*limy2(2),'*k');
%     end


    set(fig1, 'Visi', 'on')
%     if pnum==2
%         print(fig1, fnam, '-depsc')
%         close(fig1)
%     end
    
end