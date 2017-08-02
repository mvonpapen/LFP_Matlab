%%  Plot relative and absolute (percentual) ON-OFF difference of PSDs

exp     = 'rest';
textlab = {'\beta_1', '\beta_2'};
xx      = 1:length(textlab);
fig1    = figure('PaperSize', [8 7], ...
           'PaperPositionmode', 'manual', 'PaperPosition', [1 0.5 7 6], ...
           'Visible', 'off');
tag     = 'PCC_v1_filt';
fnam    = ['dPSD_' tag '.eps'];
STNact  = 1;
patmean = false;
difftype = 'rel';
Bonnf   = 4;
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
    'P_coh_off', 'P_coh_on', ...
    'P_inc_off', 'P_inc_on', 'f', 'Patient', 'LFPact')

if patmean
    % Averaged over patients
    xtot = squeeze(nanmean(P_tot_off,2));
    ytot = squeeze(nanmean(P_tot_on,2));
    xinc = squeeze(nanmean(P_inc_off,2));
    yinc = squeeze(nanmean(P_inc_on,2));
    xcoh = squeeze(nanmean(P_coh_off,2));
    ycoh = squeeze(nanmean(P_coh_on,2));
else
    % Averaged over patients
    xtot = P_tot_off(:,:);
    ytot = P_tot_on(:,:);
    xinc = P_inc_off(:,:);
    yinc = P_inc_on(:,:);
    xcoh = P_coh_off(:,:);
    ycoh = P_coh_on(:,:);
end

if strcmp(difftype, 'rel')
    % Calculate difference of percentages
    for i = 1:Nb
        j = fband{i};
        Pi(:,i,1) = squeeze( nanmean( xinc(j,:)./xtot(j,:) ))*100;
        Pi(:,i,2) = squeeze( nanmean( yinc(j,:)./ytot(j,:) ))*100;
        dPinc(:,i) = squeeze( nanmean( ...
            yinc(j,:)./ytot(j,:) - xinc(j,:)./xtot(j,:) ))*100;
        dPcoh(:,i) = squeeze( nanmean( ...
            ycoh(j,:)./ytot(j,:) - xcoh(j,:)./xtot(j,:) ))*100;
        limy    = [-40 40];
        yticks  = -40:10:40;
        ylab = '$\Delta(P_i/P_\mathrm{tot})$ [\%]';
    end
end

if strcmp(difftype, 'relRest')
    % Calculate difference of percentages
    for i = 1:Nb
        j = fband{i};
        dPinc(:,i) = squeeze( nanmean( ...
            (yinc(j,:) - xinc(j,:))./xtot(j,:) ))*100;
        dPcoh(:,i) = squeeze( nanmean( ...
            (ycoh(j,:) - xcoh(j,:))./xtot(j,:) ))*100;
        limy    = [-100 60];
        yticks  = -100:20:60;
        ylab = '$(\Delta P_i) / P_\mathrm{tot,off}$ [\%]';
    end
end

if strcmp(difftype, 'abs')
    % Calculate absolute differences
    for i = 1:Nb
        j = fband{i};
        dPinc(:,i) = squeeze( nanmean(yinc(j,:) - xinc(j,:)) );
        dPcoh(:,i) = squeeze( nanmean(ycoh(j,:) - xcoh(j,:)) );
        limy    = [-10 6]*1e-8;
        yticks  = (-10:2:6)*1e-8;
        ylab = '$\Delta P_i$ [V$^2$/Hz]';
    end
end


% Wilcoxon sign-rank test
for i=1:Nb
    [sig_inc(i),~,tmp] = signrank ( dPinc(:,i)  );
    Winc(i) = tmp.signedrank;
    [sig_coh(i),~,tmp]  = signrank ( dPcoh(:,i) );
    Wcoh(i) = tmp.signedrank;
end

% PCC
clear XX
XX(:,:,1)=dPinc; XX(:,:,2)=dPcoh;
boxPlot(xx, XX, 'showScatter', true, 'groupwidth', 0.7, 'xSpacing', 'x', ...
    'symbolMarker', '+', 'useMean', usemean, ...
    'showlegend', true, 'grouplabels', {'Incoherent', 'Coherent'}, ...
    'boxcolor', {[1 0.5 0.5], [0.5 0.5 1]}, ...
    'LineStyle', '--', 'Linewidth', 0.7, 'MedianColor', 'k', ...
    'ScatterMarker', '.', 'ScatterColor', 'k', ...
    'symbolcolor', {[0.7 0 0], [0 0 0.7]}, 'outlierSize', 13);
box on; grid on
hold all
set(gca, 'YTick', yticks);
set(gca, 'XTickLabel', textlab);
ylim(limy), xlim([0.5 Nb+0.5]), grid on
ylabel(ylab, 'Interp', 'Latex')
title('PCC', 'Interp', 'Latex');
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
end
    

set(fig1, 'Visi', 'on')
% print(fig1, fnam, '-depsc')
% close(fig1)