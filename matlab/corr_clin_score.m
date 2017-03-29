%%  Calculate correlation of dP with rigor or updrs

clear CC PP fband sig_inc sig_coh sig_vc sig_bip dP_inc dP_coh dP_vc dP_bip

% Parameters
tag     = 'v1';
STNact  = 0;
cutoff  = 1e-10;
score   = 'rigor'; % rigor/updrs


% Define frequency bandsclear fband
load(['LFP_pow_' tag '_rest.mat'], 'f')
fband{1} = find( f>=13 & f<=20 );
fband{2} = find( f>=20 & f<=30 );
fband{3} = find( f>=30 & f<=40 ); % & ~(f>40 & f<60) );
fband{4} = find( f>=60 & f<=90 );
Nb = length(fband);

% Clinical scores
SCORE(:,1) = [100 100 80 30 100 50 20 50]; % rigor
SCORE(:,2) = [38  40  55 41 45  47 30 43]; % updrs


%%%%%%%% vvv ABSOLUTE DIFFERENCES vvv %%%%%%%%%%%%%%%

exp     = {'rest', 'hold', 'fist'};
% fig1    = figure('PaperSize', [8 7], ...
%            'PaperPositionmode', 'manual', 'PaperPosition', [1 0.5 7 6], ...
%            'Visible', 'on');
       
       
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
    toty = P_tot_on(:,:);
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
%         test_coh_tot{ii}(:,i) = squeeze( nanmean( ycoh(j,:) ) ...
%                                 ./ nanmean( toty(j,:)) )*100;
    end
    
    
    % Wilcoxon sign-rank test
    for i=1:Nb
        sig_inc(i,ii) = signrank ( dP_inc{ii}(:,i) ) * Nb;
        sig_coh(i,ii) = signrank ( dP_coh{ii}(:,i) ) * Nb;
        sig_vc(i,ii)  = signrank ( dP_vc{ii}(:,i)  ) * Nb;
        sig_bip(i,ii) = signrank ( dP_bip{ii}(:,i) ) * Nb;
%         sig_test_coh_tot(i,ii) = signrank ( test_coh_tot{ii}(:,i) ) * Nb;
    end
    
    
    %% Correlation with rigor
    if strcmp(score,'rigor')
        x = SCORE(:,1)';
    elseif strcmp(score,'updrs')
        x = SCORE(:,2)';
    end
    switch ii
        case 2
            x = x([1 2 3 4 5 7 8]);
        case 3
            x = x([1 2 4 5 7 8]);
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
%     for i=1:Nb
%         j=~isnan(test_coh_tot{ii}(:,i));
%         [CC(5,i,ii), PP(5,i,ii)] = corr(X(j)',test_coh_tot{ii}(j,i),'type', 'Spearman');
%     end

end

i = fdr(PP);
[i PP(i) CC(i)],


% %% Plot
% ii=2; %task
% fb=4; %freq.band
% xx=[25:60]; %x-values
% x=SCORE([1 2 3 4 5 7 8],2)'; %score (1=rigor, 2=updrs)
% 
% y=dP_vc{ii}(:,fb); 
% X = [];
% for i=1:length(x)
%     X  = [X x([i i i i i i])];
% end
% figure, plot(X', y, 'k+');
% j=~isnan(y);
% [p,S] = polyfit(X(j),y(j)',1);
% [Y,dY] = polyconf(p,xx,S,'predopt','curve');
% hold all
% plot(xx,Y,'-k',xx,Y-dY,'--k',xx,Y+dY,'--k')