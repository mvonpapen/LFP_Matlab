%%  Plot relative and absolute (percentual) ON-OFF difference of PSDs

exp     = 'rest';
textlab = {'\beta_1', '\beta_2'};
limy    = [-30 30];
yticks  = [-30 -20 -10 0 10 20 30];
xx      = 1:length(textlab);
fig1    = figure('PaperSize', [8 7], ...
           'PaperPositionmode', 'manual', 'PaperPosition', [1 0.5 7 6], ...
           'Visible', 'off');
tag     = 'v3';
fnam    = ['dPSD_' tag '.eps'];
STNact  = 1;
patmean = false;
difftype = 'rel';
Bonnf   = 6;
usemean = false; %false = median



%%%%%%%% vvv ABSOLUTE DIFFERENCES vvv %%%%%%%%%%%%%%%

clear CC PP fband sig_inc sig_coh sig_vc sig_bip ...
    dP_inc dP_coh dP_vc dP_bip Winc Wcoh Wvc

load(['LFP_pow_' tag '_rest.mat'], 'f')
% Define frequency bandsclear fband
fband{1} = find( f>=13 & f<=20);
fband{2} = find( f>=20 & f<=30);
Nb = length(fband);
xticks  = 1:Nb;


task = exp;
load(['LFP_pow_' tag '_' task '.mat'], 'P_tot_off', 'P_tot_on', ...
    'P_coh_off', 'P_coh_on', 'P_vc_off', 'P_vc_on', ...
    'P_inc_off', 'P_inc_on', 'f', 'Patient', 'LFPact')

if patmean
    % Averaged over patients
    xtot = squeeze(nanmean(P_tot_off,2));
    ytot = squeeze(nanmean(P_tot_on,2));
    xinc = squeeze(nanmean(P_inc_off,2));
    yinc = squeeze(nanmean(P_inc_on,2));
    xcoh = squeeze(nanmean(P_coh_off,2));
    ycoh = squeeze(nanmean(P_coh_on,2));
    xvc  = squeeze(nanmean(P_vc_off,2));
    yvc  = squeeze(nanmean(P_vc_on,2));
else
    % Averaged over patients
    xtot = P_tot_off(:,:);
    ytot = P_tot_on(:,:);
    xinc = P_inc_off(:,:);
    yinc = P_inc_on(:,:);
    xcoh = P_coh_off(:,:);
    ycoh = P_coh_on(:,:);
    xvc  = P_vc_off(:,:);
    yvc  = P_vc_on(:,:);
end

if strcmp(difftype, 'rel')
    % Calculate difference of percentages
    for i = 1:Nb
        j = fband{i};
        dP_inc(:,i) = squeeze( nanmean( ...
            yinc(j,:)./ytot(j,:) - xinc(j,:)./xtot(j,:) ))*100;
        dP_coh(:,i) = squeeze( nanmean( ...
            ycoh(j,:)./ytot(j,:) - xcoh(j,:)./xtot(j,:) ))*100;
        dP_vc(:,i)  = squeeze( nanmean( ...
            yvc(j,:)./ytot(j,:) - xvc(j,:)./xtot(j,:) ))*100;
    end
end

if strcmp(difftype, 'relOFF')
    % Calculate difference as percentage of OFF power
    for i = 1:Nb
        j = fband{i};
        dP_inc(:,i) = squeeze( nanmean( ...
            (yinc(j,:) - xinc(j,:)) ./ xtot(j,:) ))*100;
        dP_coh(:,i) = squeeze( nanmean( ...
            (ycoh(j,:) - xcoh(j,:)) ./ xtot(j,:) ))*100;
        dP_vc(:,i)  = squeeze( nanmean( ...
            (yvc(j,:) - xvc(j,:)) ./ xtot(j,:) ))*100;
    end
end


% Wilcoxon sign-rank test
for i=1:Nb
    [sig_inc(i),~,tmp] = signrank ( dP_inc(:,i)  );
    Winc(i) = tmp.signedrank;
    [sig_coh(i),~,tmp]  = signrank ( dP_coh(:,i) );
    Wcoh(i) = tmp.signedrank;
    [sig_vc(i),~,tmp]  = signrank ( dP_vc(:,i) );
    Wvc(i) = tmp.signedrank;
end

% PCC
clear XX
XX(:,:,1)=dP_inc; XX(:,:,2)=dP_coh; XX(:,:,3)=dP_vc;
%     subplot(1,Nb,ii)
if ii==1
    leg = true;
else
    leg = false;
end
boxPlot(xx, XX, 'showScatter', true, 'groupwidth', 0.7, 'xSpacing', 'x', ...
    'symbolMarker', '+', 'useMean', usemean, ...
    'showlegend', leg, 'grouplabels', {'Incoherent', 'Coherent', 'Vol.-cond.'}, ...
    'boxcolor', {[1 0.5 0.5], [0.5 0.5 1], [0.5 1 0.5]}, ...
    'LineStyle', '--', 'Linewidth', 0.7, 'MedianColor', 'k', ...
    'ScatterMarker', '.', 'ScatterColor', 'k', ...
    'symbolcolor', {[0.7 0 0], [0 0 0.7], [0 0.7 0]}, 'outlierSize', 13);
box on; grid on
hold all
set(gca, 'YTick', yticks);
set(gca, 'XTickLabel', textlab);
ylim(limy), xlim([0.5 Nb+0.5]), grid on
ylabel(['$\Delta(P_i/P_\mathrm{tot})$ [\%] (ON-OFF)'], 'Interp', 'Latex')
title(task, 'Interp', 'Latex');
for i=1:Nb
    if sig_inc(i)*Bonnf<0.05;        
        plot([i-0.03 i+0.03]-0.25,0.9*limy([2 2]),'*k');
    elseif sig_inc(i)<0.05;
        plot(i-0.25,0.9*limy(2),'*k');
    end
    if sig_coh(i)*Bonnf<0.05;        
        plot([i-0.03 i+0.03],0.9*limy([2 2]),'*k');
    elseif sig_coh(i)<0.05;
        plot(i,0.9*limy(2),'*k');
    end
    if sig_vc(i)*Bonnf<0.05;        
        plot([i-0.03 i+0.03]+0.25,0.9*limy([2 2]),'*k');
    elseif sig_vc(i)<0.05;
        plot(i+0.25,0.9*limy(2),'*k');
    end
end
    
    
%     %% Correlation with rigor
%     rigor = [100 100 80 30 100 50 20 50];
%     updrs = [38  40  55 41 45  47 30 43];
%     score = updrs;
%     switch ii
%         case 1
%             x = score;
%         case 2
%             x = score([1 2 3 4 5 7 8]);
%         case 3
%             x = score([1 2 4 5 7 8]);
%     end
%     X = []; X2 =[];
%     for i=1:length(x);
%         X  = [X x([i i i])]; 
%         X2 = [X2 x([i i i i i])]; 
%     end
%     
%     for i=1:Nb
%         j=~isnan(dP_inc{ii}(:,i));
%         [CC(1,i,ii), PP(1,i,ii)] = corr(X(j)',dP_inc{ii}(j,i),'type', 'Spearman');
%     end
%     for i=1:Nb
%         j=~isnan(dP_coh{ii}(:,i));
%         [CC(2,i,ii), PP(2,i,ii)] = corr(X(j)',dP_coh{ii}(j,i),'type', 'Spearman');
%     end
%     for i=1:Nb
%         j=~isnan(dP_vc{ii}(:,i));
%         [CC(3,i,ii), PP(3,i,ii)] = corr(X(j)',dP_vc{ii}(j,i), 'type', 'Spearman');
%     end


% i = fdr(PP);
% [i PP(i) CC(i)],

set(fig1, 'Visi', 'on')
% print(fig1, fnam, '-depsc')
% close(fig1)