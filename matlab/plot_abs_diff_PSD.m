%%  Plot relative and absolute (percentual) ON-OFF difference of PSDs

nch = 4;
cutoff = 1e-10;
exp = {'rest', 'hold', 'fist'};
textlab = {'\beta_1', '\beta_2', '\gamma'};
% textlab = {'\theta', '\alpha', '\beta_1', '\beta_2', '\gamma'};
% textlab = {'\delta', '\theta', '\alpha', '\beta', '\beta1', '\beta2', ...
%         '\gamma', 'HF'};
limy = [-50 50];
limy2= [-150 150];
yticks = [-50 -25 0 25 50];
fig1 = figure('PaperSize', [14 8], ...
        'PaperPositionmode', 'manual', 'PaperPosition', [1 1 12 6], ...
        'Visible', 'off');
tag  = 'test_ds10_pc15';
fnam = ['psd_diff_' tag '.eps'];
ncomp= 3; %length(textlab)-1;



%%%%%%%% vvv ABSOLUTE DIFFERENCES vvv %%%%%%%%%%%%%%%

for ii = 1:3
    
    type = exp{ii};
    load(['pow_coh_' tag '_w0_12_' type '.mat'], ...
        'P_coh_off', 'P_coh_on', 'P_vc_off', 'P_vc_on', 'P_bip_off', ...
        'P_bip_on', 'P_inc_off', 'P_inc_on', 'f', 'Patient')
%     if ii==1
%         load(['pow_coh_2nd_w0_12_' type '.mat'], 'P_tot_off', 'P_tot_on', ...
%             'P_coh_off', 'P_coh_on', 'P_vc_off', 'P_vc_on', 'P_bip_off', ...
%             'P_bip_on', 'P_inc_off', 'P_inc_on', 'f', 'Patient')
%     end
    
%     P_inc_off = P_tot_off - P_coh_off - P_vc_off;
%     P_inc_on  = P_tot_on - P_coh_on - P_vc_on;
    val = cutoff;
    P_tot_off = P_coh_off + P_inc_off + P_vc_off;
    P_tot_on  = P_coh_on + P_inc_on + P_vc_on;
    P_inc_off(P_inc_off<=cutoff) = val;
    P_inc_on(P_inc_on<=cutoff)   = val;
    P_coh_off(P_coh_off<=cutoff) = val;
    P_coh_on(P_coh_on<=cutoff)   = val;
    P_vc_off(P_vc_off<=cutoff)   = val;
    P_vc_on(P_vc_on<=cutoff)     = val;
    

    % Define frequency bands
    clear fband
%     fband{1} = find(f>=1  & f<=4 );
%     fband{1} = find(f>=4  & f<=7 );
%     fband{2} = find(f>=7  & f<=13);
%     fband{1} = find(f>=13 & f<=30);
    fband{1} = find(f>=30 & f<=45);
    fband{2} = find(f>=20 & f<=35);
    fband{3} = find(f>=60 & f<=85);
%     fband{8} = find(f>=250 & f<=350);
    Nb = length(fband);
    

%     % Patient-specific actions
%     for i=1:length(Patient)
% %         tmp(i)=any(strcmp(Patient{i}, {'P01_L', 'P02_R', 'P04_L'}));
%         tmp(i)=any(strcmp(Patient{i}, {'P03_L', 'P05_L', 'P01_L', 'P02_R', 'P04_L'}));
%     end
%     P_tot_off(:,:,tmp) = NaN;
%     P_bip_off(:,:,tmp) = NaN;
%     
%     load('pow_coh_P02_R', 'P_tot', 'P_inc', 'P_coh', 'P_vc', 'P_bip')
%     [~, nch, ~] = size(P_tot);
%     P_tot_off(:,1:nch,1) = P_tot(:,:,2);
%     P_inc_off(:,1:nch,1) = P_inc(:,:,2);
%     P_coh_off(:,1:nch,1) = P_coh(:,:,2);
%     P_vc_off(:,1:nch,1)  = P_vc(:,:,2);
%     [~, nch, ~] = size(P_bip);
%     P_bip_off(:,1:nch,1) = P_bip(:,:,2);
%     if ii==1 || ii==2
%         load('pow_coh_P04_L', 'P_tot', 'P_inc', 'P_coh', 'P_vc', 'P_bip')
%         [~, nch, ~] = size(P_tot);
%         P_tot_off(:,1:nch,3) = P_tot(:,:,2);
%         P_inc_off(:,1:nch,3) = P_inc(:,:,2);
%         P_coh_off(:,1:nch,3) = P_coh(:,:,2);
%         P_vc_off(:,1:nch,3)  = P_vc(:,:,2);
%         [~, nch, ~] = size(P_bip);
%         P_bip_off(:,1:nch,3) = P_bip(:,:,2);
%     end
%     switch ii
%         case 1
%             n = 7;
%         case 2
%             n = 6;
%         case 3
%             n = 5;
%     end
%     load('pow_coh_P01_L', 'P_tot', 'P_inc', 'P_coh', 'P_vc', 'P_bip')
%     [~, nch, ~] = size(P_tot);
%     P_tot_off(:,1:nch,n) = P_tot(:,:,2);
%     P_inc_off(:,1:nch,n) = P_inc(:,:,2);
%     P_coh_off(:,1:nch,n) = P_coh(:,:,2);
%     P_vc_off(:,1:nch,n)  = P_vc(:,:,2);
%     [~, nch, ~] = size(P_bip);
%     P_bip_off(:,1:nch,n) = P_bip(:,:,2);
%     clear P_tot P_inc P_coh P_vc P_bip n
    
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
%     if ii==1
%         restoff    = P_tot_off(:,:,[1 2 4 5 7 8]);
%         reston     = P_tot_on(:,:,[1 2 4 5 7 8]);
%         restincoff = P_tot_off(:,:);
%         restincon  = P_tot_on(:,:);
%         restcohoff = P_tot_off(:,:,[1 2 4 5 7 8]);
%         restcohon  = P_tot_on(:,:,[1 2 4 5 7 8]);
%         restvcoff  = P_tot_off(:,:);
%         restvcon   = P_tot_on(:,:);
%     end
    clear dPabs_inc dPabs_coh dPabs_vc dPabs_bip dP_incoh_OFF dP_incoh_ON tmp
    for i = 1:Nb
        dPabs_inc(:,i) = squeeze( nanmean( yinc(fband{i},:)-xinc(fband{i},:), 1 ) ...
                                ./ nanmean(totx(fband{i},:), 1 ) )*100;
        dPabs_coh(:,i) = squeeze( nanmean( ycoh(fband{i},:)-xcoh(fband{i},:), 1 ) ...
                                ./ nanmean(totx(fband{i},:), 1 ) )*100;
        dPabs_vc(:,i)  = squeeze( nanmean( yvc(fband{i},:) -xvc(fband{i},:),  1 ) ...
                                ./ nanmean(totx(fband{i},:),  1 ) )*100;
        dPabs_bip(:,i) = squeeze( nanmean( ybip(fband{i},:)-xbip(fband{i},:), 1 ) ...
                                ./ nanmean(xbip(fband{i},:), 1 ) )*100;
        dP_incoh_OFF(:,i)= squeeze( nanmean( xinc(fband{i},:), 1 ) ...
                                ./ nanmean(xcoh(fband{i},:), 1 ) );
        dP_incoh_ON(:,i) = squeeze( nanmean( yinc(fband{i},:), 1 ) ...
                                ./ nanmean(ycoh(fband{i},:), 1 ) );
    end
%     dPabs_inc(abs(dPabs_inc)==Inf) = NaN;
%     dPabs_coh(abs(dPabs_coh)==Inf) = NaN;
%     dPabs_vc(abs(dPabs_vc)==Inf) = NaN;
%     dPabs_bip(abs(dPabs_bip)==Inf) = NaN;
    clear xinc yinc xcoh ycoh xvc yvc totx toty
    
    
    % Wilcoxon sign-rank test
    for i=1:Nb
        sig_inc(i,ii) = signrank ( dPabs_inc(:,i) ) * ncomp;
        sig_coh(i,ii) = signrank ( dPabs_coh(:,i) ) * ncomp;
        sig_vc(i,ii)  = signrank ( dPabs_vc(:,i)  ) * ncomp;
        sig_bip(i,ii) = signrank ( dPabs_bip(:,i) ) * ncomp;
    end

    if ii==1
        Pi = dPabs_inc;
        Pc = dPabs_coh;
        Pv = dPabs_vc;
    end

    % Monopolar inc
    subplot(3,4,(ii-1)*4+1)
    boxplot( dPabs_inc, 'labels', textlab, 'whisker', 50)
    set(gca,'TickLabelInterpreter','tex');
    set(gca, 'YTick', yticks);
    ylim(limy), grid on
    ylabel({['{\bf\fontsize{14}' type '}']; '\DeltaP/P_{OFF} [%]'})
    if ii==1
        title('Incoherent', 'fontw', 'bold', 'fontsi', 14)
    end
    hold all
    for i=1:Nb
        if sig_inc(i,ii)<0.01;        
            plot([i-0.07 i+0.07],0.9*limy([2 2]),'*k');
        elseif sig_inc(i,ii)<0.05;
            plot(i,0.9*limy(2),'*k');
        end
    end
    
    % Monopolar coherent
    subplot(3,4,(ii-1)*4+2)
    boxplot( dPabs_coh, 'labels', textlab, 'whisker', 50)
    set(gca,'TickLabelInterpreter','tex');
    set(gca, 'YTick', yticks);
    ylim(limy), grid on
    if ii==1
        title('Coherent', 'fontw', 'bold', 'fontsi', 14)
    end
    hold all
    for i=1:Nb
        if sig_coh(i,ii)<0.01;        
            plot([i-0.07 i+0.07],0.9*limy([2 2]),'*k');
        elseif sig_coh(i,ii)<0.05;
            plot(i,0.9*limy(2),'*k');
        end
    end
    
    % Monopolar volume conduction
    subplot(3,4,(ii-1)*4+3)
    boxplot( dPabs_vc, 'labels', textlab, 'whisker', 50)
    set(gca,'TickLabelInterpreter','tex');
    set(gca, 'YTick', yticks);
    ylim(limy), grid on
    if ii==1
        title('Volume Conduction', 'fontw', 'bold', 'fontsi', 14)
    end
    hold all
    for i=1:Nb
        if sig_vc(i,ii)<0.01;        
            plot([i-0.07 i+0.07],0.9*limy([2 2]),'*k');
        elseif sig_vc(i,ii)<0.05;
            plot(i,0.9*limy(2),'*k');
        end
    end
    
    % Bipolar
    subplot(3,4,(ii-1)*4+4)
    boxplot( dPabs_bip, 'labels', textlab, 'whisker', 50)
    set(gca,'TickLabelInterpreter','tex');
    set(gca, 'YTick', [-200 -150 -100 -50 0 50 100 150 200]);
    ylim(limy2), grid on
    if ii==1
        title('Bipolar', 'fontw', 'bold', 'fontsi', 14)
    end
    hold all
    for i=1:Nb
        if sig_bip(i,ii)<0.01;        
            plot([i-0.07 i+0.07],0.9*limy2([2 2]),'*k');
        elseif sig_bip(i,ii)<0.05;
            plot(i,0.9*limy2(2),'*k');
        end
    end

end

% text( -3, 0.5, 'Fist','HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'middle', 'Rotation', 90, ...
%     'units', 'normalized', 'Fontsize', 14, 'Fontw', 'bold')
% text( -3.8, 1.9, 'Hold','HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'middle', 'Rotation', 90, ...
%     'units', 'normalized', 'Fontsize', 14, 'Fontw', 'bold')
% text( -3.8, 3.3, 'Rest','HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'middle', 'Rotation', 90, ...
%     'units', 'normalized', 'Fontsize', 14, 'Fontw', 'bold')
% set(fig1, 'visi', 'on')
% saveas(fig1, 'psd_diff_abs_all.eps')
print(fig1, fnam, '-depsc')
close(fig1)