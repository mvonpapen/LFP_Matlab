%%  Plot relative and absolute (percentual) ON-OFF difference of PSDs

cutoff = 1e-10;
textlab = {'Inc', 'Coh', 'VC', 'Bip'};
limy = [-50 50];
fig1 = figure('Papersize', [7 3], 'PaperPosition', [0.75 0.5 5.5 2], ...
        'PaperPositionmode', 'manual', 'Visible', 'off');
STNact = 0;
tag    = 'ds10_pc15';

% Take out patients w/o FIST
pat1 = [1 2 4 5 7 8];
pat2 = 1:6;
% % pat1 = [2 4 5 8];
% % pat2 = [2 3 4 6];
task = 'fist';

% % Take out patients w/o HOLD
% pat1 = [1 2 3 4 5 7 8];
% pat2 = 1:7;
% % pat1 = [2 4 5 8];
% % pat2 = [2 4 5 7];
% task = 'hold';

fb     =  [13 20];
fb_str = '$\beta_1=13{-}20\,$Hz';
fnam   = ['psd_diff_' tag '_' task '2rest_beta.eps'];
% fb     =  [60 85];
% fb_str = '$\gamma=60{-}85\,$Hz';
% fnam   = ['psd_diff_' tag '_' task '2rest_gamma.eps'];


%%%%%%%% vvv Rest vvv %%%%%%%%%%%%%%%

load(['pow_coh_' tag '_w0_12_rest.mat'], 'P_tot_off', 'P_tot_on', ...
    'P_coh_off', 'P_coh_on', 'P_vc_off', 'P_vc_on', 'P_bip_off', ...
    'P_bip_on', 'P_inc_off', 'P_inc_on', 'f', 'LFPact')
fband  = find(f>=fb(1) & f<=fb(2));


% P_inc_off=P_tot_off-P_coh_off-P_vc_off;
% P_inc_on=P_tot_on-P_coh_on-P_vc_on;
val = cutoff;
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
xbip_rest = P_bip_off(:,:);
ybip_rest = P_bip_on(:,:);

%%%%%%%% vvv Fist vvv %%%%%%%%%%%%%%%

load(['pow_coh_' tag '_w0_12_' task '.mat'], 'P_tot_off', 'P_tot_on', ...
    'P_coh_off', 'P_coh_on', 'P_vc_off', 'P_vc_on', 'P_bip_off', ...
    'P_bip_on', 'P_inc_off', 'P_inc_on', 'f')

% P_inc_off=P_tot_off-P_coh_off-P_vc_off;
% P_inc_on=P_tot_on-P_coh_on-P_vc_on;
P_inc_off(P_inc_off<=cutoff) = 0;
P_inc_on(P_inc_on<=cutoff)   = 0;
P_coh_off(P_coh_off<=cutoff) = 0;
P_coh_on(P_coh_on<=cutoff)   = 0;
P_vc_off(P_vc_off<=cutoff)   = 0;
P_vc_on(P_vc_on<=cutoff)     = 0;

    
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
xbip_fist = P_bip_off(:,:);
ybip_fist = P_bip_on(:,:);



% %% Ratios
% % OFF
% dPabs_tot_off = log10(squeeze( nanmean( totx_fist(fband,:)./totx_rest(fband,:), 1 ) ));
% dPabs_inc_off = log10(squeeze( nanmean( xinc_fist(fband,:)./xinc_rest(fband,:), 1 ) ));
% dPabs_coh_off = log10(squeeze( nanmean( xcoh_fist(fband,:)./xcoh_rest(fband,:), 1 ) ));
% dPabs_vc_off  = log10(squeeze( nanmean( xvc_fist(fband,:) ./xvc_rest(fband,:),  1 ) ));
% dPabs_bip_off = log10(squeeze( nanmean( xbip_fist(fband,:)./xbip_rest(fband,:), 1 ) ));
% % ON
% dPabs_tot_on  = log10(squeeze( nanmean( toty_fist(fband,:)./toty_rest(fband,:), 1 ) ));
% dPabs_inc_on  = log10(squeeze( nanmean( yinc_fist(fband,:)./yinc_rest(fband,:), 1 ) ));
% dPabs_coh_on  = log10(squeeze( nanmean( ycoh_fist(fband,:)./ycoh_rest(fband,:), 1 ) ));
% dPabs_vc_on   = log10(squeeze( nanmean( yvc_fist(fband,:) ./yvc_rest(fband,:),  1 ) ));
% dPabs_bip_on  = log10(squeeze( nanmean( ybip_fist(fband,:)./ybip_rest(fband,:), 1 ) ));


%% Differences (normalized by total power)
% OFF
dPabs_tot_off  = squeeze( nanmean( totx_fist(fband,:)-totx_rest(fband,:), 1 ) ...
                        ./ nanmean(totx_rest(fband,:), 1 ) )*100;
dPabs_inc_off  = squeeze( nanmean( xinc_fist(fband,:)-xinc_rest(fband,:), 1 ) ...
                        ./ nanmean(totx_rest(fband,:), 1 ) )*100;
dPabs_coh_off = squeeze( nanmean( xcoh_fist(fband,:)-xcoh_rest(fband,:), 1 ) ...
                        ./ nanmean(totx_rest(fband,:), 1 ) )*100;
dPabs_vc_off = squeeze( nanmean( xvc_fist(fband,:) -xvc_rest(fband,:),  1 ) ...
                        ./ nanmean(totx_rest(fband,:), 1 ) )*100;
dPabs_acoh_off = squeeze( nanmean(xcoh_fist(fband,:)+yvc_fist(fband,:) ...
    -xcoh_rest(fband,:)-xvc_rest(fband,:),1) ...
    ./ nanmean(xcoh_rest(fband,:)+xvc_rest(fband,:), 1 ) )*100;
dPabs_nvc_off  = squeeze( nanmean(xcoh_fist(fband,:)+xinc_fist(fband,:) ...
    -xcoh_rest(fband,:)-xinc_rest(fband,:),1) ...
    ./ nanmean(xcoh_rest(fband,:)+xinc_rest(fband,:), 1 ) )*100;
dPabs_bip_off  = squeeze( nanmean( xbip_fist(fband,:)-xbip_rest(fband,:), 1 ) ...
                        ./ nanmean(xbip_rest(fband,:), 1 ) )*100;
% ON
dPabs_tot_on  = squeeze( nanmean( toty_fist(fband,:)-toty_rest(fband,:), 1 ) ...
                        ./ nanmean(toty_rest(fband,:), 1 ) )*100;
dPabs_inc_on  = squeeze( nanmean( yinc_fist(fband,:)-yinc_rest(fband,:), 1 ) ...
                        ./ nanmean(toty_rest(fband,:), 1 ) )*100;
dPabs_coh_on = squeeze( nanmean( ycoh_fist(fband,:)-ycoh_rest(fband,:), 1 ) ...
                        ./ nanmean(toty_rest(fband,:), 1 ) )*100;
dPabs_vc_on = squeeze( nanmean( yvc_fist(fband,:) -yvc_rest(fband,:),  1 ) ...
                        ./ nanmean(toty_rest(fband,:), 1 ) )*100;
dPabs_acoh_on = squeeze( nanmean(ycoh_fist(fband,:)+yvc_fist(fband,:) ...
    -ycoh_rest(fband,:)-yvc_rest(fband,:),1) ...
    ./ nanmean(ycoh_rest(fband,:)+yvc_rest(fband,:), 1 ) )*100;
dPabs_nvc_on  = squeeze( nanmean(ycoh_fist(fband,:)+yinc_fist(fband,:) ...
    -ycoh_rest(fband,:)-yinc_rest(fband,:),1) ...
    ./ nanmean(ycoh_rest(fband,:)+yinc_rest(fband,:), 1 ) )*100;
dPabs_bip_on  = squeeze( nanmean( ybip_fist(fband,:)-ybip_rest(fband,:), 1 ) ...
                        ./ nanmean(ybip_rest(fband,:), 1 ) )*100;

                    

% %% Differences (normalized by inc,coh,vc)
% % OFF
% dPabs_tot_off  = squeeze( nanmean( totx_fist(fband,:)-totx_rest(fband,:), 1 ) ...
%                         ./ nanmean(totx_rest(fband,:), 1 ) )*100;
% dPabs_inc_off  = squeeze( nanmean( xinc_fist(fband,:)-xinc_rest(fband,:), 1 ) ...
%                         ./ nanmean(xinc_rest(fband,:), 1 ) )*100;
% dPabs_coh_off = squeeze( nanmean( xcoh_fist(fband,:)-xcoh_rest(fband,:), 1 ) ...
%                         ./ nanmean(xcoh_rest(fband,:), 1 ) )*100;
% dPabs_vc_off = squeeze( nanmean( xvc_fist(fband,:) -xvc_rest(fband,:),  1 ) ...
%                         ./ nanmean(xvc_rest(fband,:), 1 ) )*100;
% dPabs_acoh_off = squeeze( nanmean(xcoh_fist(fband,:)+yvc_fist(fband,:) ...
%     -xcoh_rest(fband,:)-xvc_rest(fband,:),1) ./ nanmean(xcoh_rest(fband,:)+xvc_rest(fband,:),1) )*100;
% dPabs_nvc_off  = squeeze( nanmean(xcoh_fist(fband,:)+xinc_fist(fband,:) ...
%     -xcoh_rest(fband,:)-xinc_rest(fband,:),1) ./ nanmean(xcoh_rest(fband,:)-xinc_rest(fband,:),1) )*100;
% dPabs_bip_off  = squeeze( nanmean( xbip_fist(fband,:)-xbip_rest(fband,:), 1 ) ...
%                         ./ nanmean(xbip_rest(fband,:), 1 ) )*100;
% % ON
% dPabs_tot_on  = squeeze( nanmean( toty_fist(fband,:)-toty_rest(fband,:), 1 ) ...
%                         ./ nanmean(toty_rest(fband,:), 1 ) )*100;
% dPabs_inc_on  = squeeze( nanmean( yinc_fist(fband,:)-yinc_rest(fband,:), 1 ) ...
%                         ./ nanmean(yinc_rest(fband,:), 1 ) )*100;
% dPabs_coh_on = squeeze( nanmean( ycoh_fist(fband,:)-ycoh_rest(fband,:), 1 ) ...
%                         ./ nanmean(ycoh_rest(fband,:), 1 ) )*100;
% dPabs_vc_on = squeeze( nanmean( yvc_fist(fband,:) -yvc_rest(fband,:),  1 ) ...
%                         ./ nanmean(yvc_rest(fband,:), 1 ) )*100;
% dPabs_acoh_on = squeeze( nanmean(ycoh_fist(fband,:)+yvc_fist(fband,:) ...
%     -ycoh_rest(fband,:)-yvc_rest(fband,:),1) ./ nanmean(toty_rest(fband,:), 1 ) )*100;
% dPabs_nvc_on  = squeeze( nanmean(ycoh_fist(fband,:)+yinc_fist(fband,:) ...
%     -ycoh_rest(fband,:)-yinc_rest(fband,:),1) ./ nanmean(toty_rest(fband,:), 1 ) )*100;
% dPabs_bip_on  = squeeze( nanmean( ybip_fist(fband,:)-ybip_rest(fband,:), 1 ) ...
%                         ./ nanmean(ybip_rest(fband,:), 1 ) )*100;

clear xinc yinc xcoh ycoh xvc yvc totx toty
dPabs_inc_off(dPabs_inc_off==Inf) = NaN;
dPabs_inc_on(dPabs_inc_on==Inf) = NaN;
dPabs_coh_off(dPabs_coh_off==Inf) = NaN;
dPabs_coh_on(dPabs_coh_on==Inf) = NaN;
    


% Wilcoxon sign-rank test
sig_off(1)  = signrank ( dPabs_inc_off  );
sig_off(2)  = signrank ( dPabs_coh_off );
sig_off(3)  = signrank ( dPabs_vc_off );
sig_off(4)  = signrank ( dPabs_bip_off  );

sig_on(1)  = signrank ( dPabs_inc_on  );
sig_on(2)  = signrank ( dPabs_coh_on );
sig_on(3)  = signrank ( dPabs_vc_on );
sig_on(4)  = signrank ( dPabs_bip_on  );
    



%% Group results
X = [dPabs_inc_off dPabs_coh_off dPabs_vc_off dPabs_bip_off];
Y = [dPabs_inc_on dPabs_coh_on dPabs_vc_on dPabs_bip_on];
G = zeros(1,length(X));
n  = length(dPabs_inc_off);
G(1:n) = 1;
G(n+1:2*n) = 2;
G(2*n+1:3*n) = 3;
G(3*n+1:end) = 4;

% Plot OFF
subplot(1,2,1);
boxplot( X, G, 'labels', textlab, 'whisker', 100)
ylim(limy), grid on
ylabel(['($P^\mathrm{' task '}_i-P^\mathrm{rest}_i)/P^\mathrm{rest}_i\,[\%]$'], 'Interpreter', 'latex')
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
ylabel(['($P^\mathrm{' task '}_i-P^\mathrm{rest}_i)/P^\mathrm{rest}_i\,[\%]$'], 'Interpreter', 'latex')
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