%%  Plot relative and absolute (percentual) ON-OFF difference of PSDs

nch     = 4;
cutoff  = 1e-10;
exp     = {'rest', 'hold', 'fist'};
textlab = {'\beta_1', '\beta_2', '\gamma_1', '\gamma_2'};
limy    = [-35 35];
limy2   = [-80 80];
xticks  = 1:length(textlab);
yticks  = [-30 -20 -10 0 10 20 30];
yticks2 = [-80 -40 0 40 80];
xx      = 1:length(textlab);
fig1    = figure('PaperSize', [8 7], ...
           'PaperPositionmode', 'manual', 'PaperPosition', [1 0.5 7 6], ...
           'Visible', 'off');
tag     = 'v1';
fnam    = ['dPSD_' tag '.eps'];
STNact  = 0;
usemean = false; %false = median



%%%%%%%% vvv ABSOLUTE DIFFERENCES vvv %%%%%%%%%%%%%%%

clear CC PP fband sig_inc sig_coh sig_vc sig_bip dP_inc dP_coh dP_vc dP_bip

load(['LFP_pow_' tag '_rest.mat'], 'f')
% Define frequency bandsclear fband
fband{1} = find( f>=13 & f<=20 );
fband{2} = find( f>=20 & f<=30 );
fband{3} = find( f>=30 & f<=40 ); % & ~(f>40 & f<60) );
fband{4} = find( f>=60 & f<=90 );
Nb = length(fband);

for ii = 1:3
    
    task = exp{ii};
    load(['LFP_pow_' tag '_' task '.mat'], 'P_tot_off', 'P_tot_on', ...
        'P_coh_off', 'P_coh_on', 'P_vc_off', 'P_vc_on', 'P_bip_off', ...
        'P_bip_on', 'P_inc_off', 'P_inc_on', 'f', 'Patient', 'LFPact')
    
    val = cutoff;
%     P_tot_off = P_coh_off + P_inc_off + P_vc_off;
%     P_tot_on  = P_coh_on + P_inc_on + P_vc_on;
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
        P_vc_off(:,LFPact{i}<STNact,i)  = NaN;
        P_bip_off(:,LFPact{i}<STNact,i) = NaN;
    end
    

    
    % Normalized absolute differences of PSDs
    totx = P_tot_off(:,:);
    xinc = P_inc_off(:,:);
    yinc = P_inc_on(:,:);
    xcoh = P_coh_off(:,:);
    ycoh = P_coh_on(:,:);
    xvc  = P_vc_off(:,:);
    yvc  = P_vc_on(:,:);
    xbip = P_bip_off(:,:);
    ybip = P_bip_on(:,:);
    
    for i = 1:Nb
        j = fband{i};
        dP_inc{ii}(:,i) = squeeze( nanmean( yinc(j,:)-xinc(j,:) ) ...
                                ./ nanmean( totx(j,:)) )*100;
        dP_coh{ii}(:,i) = squeeze( nanmean( ycoh(j,:)-xcoh(j,:) ) ...
                                ./ nanmean( totx(j,:)) )*100;
        dP_vc{ii}(:,i)  = squeeze( nanmean( yvc(j,:)-xvc(j,:) ) ...
                                ./ nanmean( totx(j,:)) )*100;
        dP_bip{ii}(:,i) = squeeze( nanmean( ybip(j,:)-xbip(j,:) ) ...
                                ./ nanmean( xbip(j,:)) )*100;
    end
    
    
    % Wilcoxon sign-rank test
    for i=1:Nb
        sig_inc(i,ii) = signrank ( dP_inc{ii}(:,i) ) * Nb;
        sig_coh(i,ii) = signrank ( dP_coh{ii}(:,i) ) * Nb;
        sig_vc(i,ii)  = signrank ( dP_vc{ii}(:,i)  ) * Nb;
        sig_bip(i,ii) = signrank ( dP_bip{ii}(:,i) ) * Nb;
    end

    % PCC
    clear XX
    XX(:,:,1)=dP_inc{ii}; XX(:,:,2)=dP_coh{ii}; XX(:,:,3)=dP_vc{ii};
    subplot(3,3,(ii-1)*3+[1 2])
    boxPlot(xx, XX, 'showScatter', true, 'groupwidth', 0.7, 'xSpacing', 'x', ...
        'symbolMarker', '+', 'useMean', usemean, ...
        'showlegend', true, 'grouplabels', {'Incoherent', 'Coherent', 'Vol.-cond.'}, ...
        'boxcolor', {[1 0.5 0.5], [0.5 0.5 1], [0.5 1 0.5]}, ...
        'LineStyle', '--', 'Linewidth', 0.7, 'MedianColor', 'k', ...
        'ScatterMarker', '.', 'ScatterColor', 'k', ...
        'symbolcolor', {[0.7 0 0], [0 0 0.7], [0 0.7 0]}, 'outlierSize', 13);
    box on; grid on
    hold all
    set(gca, 'YTick', yticks);
    set(gca, 'XTickLabel', textlab);
    ylim(limy), xlim([0.5 Nb+0.5]), grid on
    ylabel({['{\bf\fontsize{14}' task '}']; '\DeltaP/P_{OFF} [%]'})
    if ii==1
        title('PCC', 'fontw', 'bold', 'fontsi', 14)
    end
%     if ii==3
%         legend('Incoherent', 'Coherent', 'Vol.-cond.')
%     end
    for i=1:Nb
        if sig_inc(i,ii)<0.01;        
            plot([i-0.03 i+0.03]-0.25,0.9*limy([2 2]),'*k');
        elseif sig_inc(i,ii)<0.05;
            plot(i-0.25,0.9*limy(2),'*k');
        end
        if sig_coh(i,ii)<0.01;        
            plot([i-0.03 i+0.03],0.9*limy([2 2]),'*k');
        elseif sig_coh(i,ii)<0.05;
            plot(i,0.9*limy(2),'*k');
        end
        if sig_vc(i,ii)<0.01;        
            plot([i-0.03 i+0.03]+0.25,0.9*limy([2 2]),'*k');
        elseif sig_vc(i,ii)<0.05;
            plot(i+0.25,0.9*limy(2),'*k');
        end
    end
    
    
    % Bipolar
    subplot(3,3,(ii-1)*3+3)
    XX = dP_bip{ii};
    boxPlot(xx, XX, 'showScatter', true, 'groupwidth', 0.5, 'xSpacing', 'x', ...
        'symbolMarker', '+', 'useMean', usemean, ...
        'LineStyle', '--', 'Linewidth', 0.7, 'MedianColor', 'k', ...
        'ScatterMarker', '.', 'outlierSize', 13, 'symbolColor', 'k');
    box on; grid on
    hold all
    set(gca,'TickLabelInterpreter','tex');
    set(gca, 'XTick', xticks);
    set(gca, 'YTick', yticks2);
    set(gca, 'XTickLabel', textlab);
    ylim(limy2), xlim([0.5 Nb+0.5]), grid on
    if ii==1
        title('Bipolar', 'fontw', 'bold', 'fontsi', 14)
    end
    hold all
    for i=1:Nb
        if sig_bip(i,ii)<0.01;        
            plot([i-0.05 i+0.05],0.9*limy2([2 2]),'*k');
        elseif sig_bip(i,ii)<0.05;
            plot(i,0.9*limy2(2),'*k');
        end
    end
    
    
    %% Correlation with rigor
    rigor = [100 100 80 30 100 50 20 50];
    updrs = [38  40  55 41 45  47 30 43];
    score = updrs;
    switch ii
        case 1
            x = score;
        case 2
            x = score([1 2 3 4 5 7 8]);
        case 3
            x = score([1 2 4 5 7 8]);
    end
    X = []; X2 =[];
    for i=1:length(x);
        X  = [X x([i i i i])]; 
        X2 = [X2 x([i i i i i i])]; 
    end
    
    for i=1:Nb
        j=~isnan(dP_inc{ii}(:,i));
        [CC(1,i,ii), PP(1,i,ii)] = corr(X(j)',dP_inc{ii}(j,i),'type', 'Spearman');
    end
    for i=1:Nb
        j=~isnan(dP_coh{ii}(:,i));
        [CC(2,i,ii), PP(2,i,ii)] = corr(X(j)',dP_coh{ii}(j,i),'type', 'Spearman');
    end
    for i=1:Nb
        j=~isnan(dP_vc{ii}(:,i));
        [CC(3,i,ii), PP(3,i,ii)] = corr(X(j)',dP_vc{ii}(j,i), 'type', 'Spearman');
    end
    for i=1:Nb
        j=~isnan(dP_bip{ii}(:,i));
        [CC(4,i,ii), PP(4,i,ii)] = corr(X2(j)',dP_bip{ii}(j,i),'type', 'Spearman');
    end

end

i = fdr(PP);
[i PP(i) CC(i)],

set(fig1, 'Visi', 'on')
% print(fig1, fnam, '-depsc')
% close(fig1)