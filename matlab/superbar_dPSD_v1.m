%%  Plot relative and absolute (percentual) ON-OFF difference of PSDs

exp     = 'rest';
textlab = {'\beta_1', '\beta_2'};
xx      = 1:length(textlab);
fig1    = figure('PaperSize', [8 7], ...
           'PaperPositionmode', 'manual', 'PaperPosition', [1 0.5 7 6], ...
           'Visible', 'off');
tag     = 'v3_bip';
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
    'P_coh_off', 'P_coh_on', 'P_vc_off', 'P_vc_on', 'P_bip_off', 'P_bip_on', ...
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
    xbip  = squeeze(nanmean(P_bip_off,2));
    ybip  = squeeze(nanmean(P_bip_on,2));
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
    xbip  = P_bip_off(:,:);
    ybip  = P_bip_on(:,:);
end

% Calculate difference (percentages for PCC and absolute for bip)
for i = 1:Nb
    j = fband{i};
    Pi(:,1,i) = squeeze( nanmean( xinc(j,:)./xtot(j,:) ))*100;
    Pi(:,2,i) = squeeze( nanmean( yinc(j,:)./ytot(j,:) ))*100;
    Pc(:,1,i) = squeeze( nanmean( xcoh(j,:)./xtot(j,:) ))*100;
    Pc(:,2,i) = squeeze( nanmean( ycoh(j,:)./ytot(j,:) ))*100;
    Pv(:,1,i) = squeeze( nanmean( xvc(j,:)./xtot(j,:) ))*100;
    Pv(:,2,i) = squeeze( nanmean( yvc(j,:)./ytot(j,:) ))*100;
    Pb(:,1,i) = squeeze( nanmean( xbip(j,:) ));
    Pb(:,2,i) = squeeze( nanmean( ybip(j,:) ));
    limy    = [-40 40];
    yticks  = -40:10:40;
    ylab = '$P_i/P_\mathrm{tot}$ [\%]';
    limy2   = [-100 100]*1e-8;
    yticks2 = (-100:25:100)*1e-8;
    ylab2 = '$P_\mathrm{bip}$';
end



%% Mean, SEM, and p-values
M = [nanmean(Pi); nanmean(Pc); nanmean(Pv); nanmean(Pb)];
E = [nanstd(Pi)/sqrt(sum(~isnan(Pi(:,1,1))));...
     nanstd(Pc)/sqrt(sum(~isnan(Pc(:,1,1))));...
     nanstd(Pv)/sqrt(sum(~isnan(Pv(:,1,1))));...
     nanstd(Pb)/sqrt(sum(~isnan(Pb(:,1,1)))) ];
% Wilcoxon sign-rank test
prob_pcc = zeros(6,6,2);
prob_bip = NaN(2,2,2);
for i=1:Nb
    [prob_pcc(1,4,i),~,tmp] = signrank ( Pi(:,1,i)-Pi(:,2,i) );
    Winc(i) = tmp.signedrank;
    [prob_pcc(2,5,i),~,tmp] = signrank ( Pc(:,1,i)-Pc(:,2,i) );
    Wcoh(i) = tmp.signedrank;
    [prob_pcc(3,6,i),~,tmp]  = signrank ( Pv(:,1,i)-Pv(:,2,i) );
    Wvc(i) = tmp.signedrank;
    [prob_bip(1,2,i),~,tmp] = signrank ( Pb(:,1,i)-Pb(:,2,i) );
    Wbip(i) = tmp.signedrank;
    prob_bip(2,1,i) = prob_bip(1,2,i);
    prob_pcc(:,:,i) = prob_pcc(:,:,i) + triu(prob_pcc(:,:,i))';
end
prob_pcc(prob_pcc==0) = NaN;

%% Plot superbar
tit={'$\beta_1 (13-20\,$Hz)'; '$\beta_2 (20-30\,$Hz)'};
for i=1:2
    subplot(2,4,[1:3]+(i-1)*4)
    % Color
    C = [.8 .2 .2; .3 .3 .9; .2 .8 .2];
    superbar(M(1:3,:,i), 'E', E(1:3,:,i), 'P', prob_pcc(:,:,i), ...
        'BarFaceColor', C, 'PStarThreshold', [0.05, 0.05/Bonnf 0 0], ...
        'PStarShowNS', false, 'PStarFontSize', 16, 'PStarColor', [0 0 0]);
    box on
    ylim([0 100]), xlim([0.5 3.5])
    xtl = {'inc', 'coh', 'vc'};
    set(gca, 'XTickLabel', xtl, 'TickLabelInterpreter', 'latex', 'XTick', 1:3)
    ylabel('$P_i/P_\mathrm{tot} [\%]$', 'Interpreter', 'Latex')
    title(tit{i}, 'Interpreter', 'Latex')
end
for i=1:2
    subplot(2,4,4+(i-1)*4)
    % Color
    C = [0 0 0];
    superbar(M(4,:,i), 'E', E(4,:,i), 'P', prob_bip(:,:,i), ...
        'BarFaceColor', C, 'PStarThreshold', [0.05, 0.05/Bonnf 0 0], ...
        'PStarShowNS', false, 'PStarFontSize', 16, 'PStarColor', [0 0 0]);
    box on
    ylim([0 200]*1e-9), xlim([0.5 2.5])
    xtl = {'bip'};
    set(gca, 'XTickLabel', xtl, 'TickLabelInterpreter', 'latex', 'XTick', 1.5)
    ylabel('$P_\mathrm{bip}$ [V$^2$/Hz]', 'Interpreter', 'Latex')
    title(['$\beta_' mat2str(i) '$'], 'Interpreter', 'Latex')
end
% set(gca, 'YTick', yticks);
% set(gca, 'XTickLabel', textlab);
% ylim(limy), xlim([0.5 Nb+0.5]), grid on
% ylabel(ylab, 'Interp', 'Latex')
% title('PCC', 'Interp', 'Latex');
    

set(fig1, 'Visi', 'on')