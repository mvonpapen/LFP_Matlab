%%  Plot relative and absolute (percentual) ON-OFF difference of PSDs

cutoff = 1e-12;
textlab = {'Inc', 'Coh', 'Bip'};
limy = [-50 50];
fig1 = figure('Papersize', [20 8], 'PaperPosition', [0.75 0.5 18.5 7], ...
        'PaperPositionmode', 'manual', 'Visible', 'off');
% fig1 = figure('Papersize', [10 3], 'PaperPosition', [0.75 0.5 8.5 2], ...
%         'PaperPositionmode', 'manual', 'Visible', 'off');
tag    = 'bip-cohinc';

% % Take out patients w/o fist
pat1 = [1 2 4 5 7 8];
pat2 = 1:6;
% pat1 = [2 4 5 8];
% pat2 = [2 3 4 6];
task = 'fist';
% % Take out patients w/o hold
% % pat1 = [1 2 3 4 5 7 8];
% % pat2 = 1:7;
% pat1 = [2 4 5 8];
% pat2 = [2 4 5 7];
% task = 'hold';

fb     =  [13 20];
fb_str = '$\beta_1=13{-}20\,$Hz';
fnam   = ['psd_diff_' tag '_' task '2rest_peak_beta.eps'];
% fb     =  [60 85];
% fb_str = '$\gamma=60{-}90\,$Hz';
% fnam   = ['psd_diff_' tag '_' task '2rest_gamma.eps'];


%%%%%%%% vvv Rest vvv %%%%%%%%%%%%%%%

load(['pow_coh_' tag '_rest.mat'], 'P_bip_off', 'P_bip_on', ...
    'P_coh_off', 'P_coh_on', 'P_inc_off', 'P_inc_on', 'f', 'LFPact')
fband  = find(f>=fb(1) & f<=fb(2));

P_inc_off(P_inc_off<=cutoff) = 0;
P_inc_on(P_inc_on<=cutoff)   = 0;
P_coh_off(P_coh_off<=cutoff) = 0;
P_coh_on(P_coh_on<=cutoff)   = 0;

    
% Specify patients
P_inc_off = P_inc_off(:,:,pat1);
P_inc_on = P_inc_on(:,:,pat1);
P_coh_off = P_coh_off(:,:,pat1);
P_coh_on = P_coh_on(:,:,pat1);
P_bip_off = P_bip_off(:,:,pat1);
P_bip_on = P_bip_on(:,:,pat1);

% Normalized absolute differences of PSDs
totx_rest = P_bip_off(:,:);
toty_rest = P_bip_on(:,:);
xinc_rest = P_inc_off(:,:);
yinc_rest = P_inc_on(:,:);
xcoh_rest = P_coh_off(:,:);
ycoh_rest = P_coh_on(:,:);

%%%%%%%% vvv Fist vvv %%%%%%%%%%%%%%%

load(['pow_coh_' tag '_' task '.mat'], 'P_bip_off', 'P_bip_on', ...
    'P_coh_off', 'P_coh_on', 'P_inc_off', 'P_inc_on', 'f', 'LFPact')

P_inc_off(P_inc_off<=cutoff) = 0;
P_inc_on(P_inc_on<=cutoff)   = 0;
P_coh_off(P_coh_off<=cutoff) = 0;
P_coh_on(P_coh_on<=cutoff)   = 0;

% Specify patients
P_inc_off = P_inc_off(:,:,pat2);
P_inc_on = P_inc_on(:,:,pat2);
P_coh_off = P_coh_off(:,:,pat2);
P_coh_on = P_coh_on(:,:,pat2);
P_bip_off = P_bip_off(:,:,pat2);
P_bip_on = P_bip_on(:,:,pat2);

% Normalized absolute differences of PSDs
totx_fist = P_bip_off(:,:);
toty_fist = P_bip_on(:,:);
xinc_fist = P_inc_off(:,:);
yinc_fist = P_inc_on(:,:);
xcoh_fist = P_coh_off(:,:);
ycoh_fist = P_coh_on(:,:);



%% Differences (normalized by total power)
% OFF
dP_tot_off  = squeeze( nanmean( totx_fist(fband,:)-totx_rest(fband,:), 1 ) ...
                        ./ nanmean(totx_rest(fband,:), 1 ) )*100;
dP_inc_off  = squeeze( nanmean( xinc_fist(fband,:)-xinc_rest(fband,:), 1 ) ...
                        ./ nanmean(totx_rest(fband,:), 1 ) )*100;
dP_coh_off  = squeeze( nanmean( xcoh_fist(fband,:)-xcoh_rest(fband,:), 1 ) ...
                        ./ nanmean(totx_rest(fband,:), 1 ) )*100;
% ON
dP_tot_on  = squeeze( nanmean( toty_fist(fband,:)-toty_rest(fband,:), 1 ) ...
                        ./ nanmean(toty_rest(fband,:), 1 ) )*100;
dP_inc_on  = squeeze( nanmean( yinc_fist(fband,:)-yinc_rest(fband,:), 1 ) ...
                        ./ nanmean(toty_rest(fband,:), 1 ) )*100;
dP_coh_on  = squeeze( nanmean( ycoh_fist(fband,:)-ycoh_rest(fband,:), 1 ) ...
                        ./ nanmean(toty_rest(fband,:), 1 ) )*100;
clear xinc yinc xcoh ycoh totx toty
    
    
% Wilcoxon sign-rank test
sig_off(1)  = signrank ( dP_inc_off );
sig_off(2)  = signrank ( dP_coh_off );
sig_off(3)  = signrank ( dP_tot_off );

sig_on(1)  = signrank ( dP_inc_on );
sig_on(2)  = signrank ( dP_coh_on );
sig_on(3)  = signrank ( dP_tot_on );
    



%% Group results
X = [dP_inc_off dP_coh_off dP_tot_off];
Y = [dP_inc_on  dP_coh_on  dP_tot_on];
G = zeros(1,length(X));
n  = length(dP_inc_off);
G(1:n) = 1;
G(n+1:2*n) = 2;
G(2*n+1:3*n) = 3;

% Plot OFF
subplot(1,2,1);
boxplot( X, G, 'labels', textlab, 'whisker', 100)
ylim(limy), grid on
ylabel(['($P^\mathrm{' task '}_i-P^\mathrm{rest}_i)/P^\mathrm{rest}_\mathrm{tot}\,[\%]$'], 'Interpreter', 'latex')
title(['OFF, ' fb_str], 'fontw', 'bold', 'fontsi', 14, 'Interpreter', 'latex')
hold all
for ii=1:length(sig_off)
    if sig_off(ii)<0.01;        
        plot([ii-0.07 ii+0.07],0.9*limy([2 2]),'*k');
    elseif sig_off(ii)<0.05;
        plot(ii,0.9*limy(2),'*k');
    end
end

% Plot ON
subplot(1,2,2);
boxplot( Y, G, 'labels', textlab, 'whisker', 100)
ylim(limy), grid on
ylabel(['($P^\mathrm{' task '}_i-P^\mathrm{rest}_i)/P^\mathrm{rest}_\mathrm{tot}\,[\%]$'], 'Interpreter', 'latex')
title(['ON, ' fb_str], 'fontw', 'bold', 'fontsi', 14, 'Interpreter', 'latex')
hold all
for ii=1:length(sig_on)
    if sig_on(ii)<0.01;        
        plot([ii-0.07 ii+0.07],0.9*limy([2 2]),'*k');
    elseif sig_on(ii)<0.05;
        plot(ii,0.9*limy(2),'*k');
    end
end
print(fig1, fnam, '-depsc')
close(fig1)