%%  Plot relative and absolute (percentual) ON-OFF difference of PSDs

nch = 4;
cutoff = 1e-20;
exp = {'rest', 'hold', 'fist'};
textlab = {'\theta', '\alpha', '\beta_1', '\beta_2', '\gamma'};
% textlab = {'\delta', '\theta', '\alpha', '\beta', '\beta1', '\beta2', ...
%         '\gamma', 'HF'};
limy = [-200 200];
limy2= [-200 200];
fig1 = figure('PaperSize', [16 10], ...
        'PaperPositionmode', 'manual', 'PaperPosition', [1 1 14 8], ...
        'Visible', 'off');
fnam = 'psd_diff_rel_nopeak_t-g_ds10_pc15.eps';



%%%%%%%% vvv ABSOLUTE DIFFERENCES vvv %%%%%%%%%%%%%%%

for ii = 1:3
    
    type = exp{ii};
    load(['pow_coh_ds10_pc15_w0_12_' type '.mat'], 'P_tot_off', 'P_tot_on', ...
        'P_coh_off', 'P_coh_on', 'P_vc_off', 'P_vc_on', 'P_bip_off', ...
        'P_bip_on', 'P_inc_off', 'P_inc_on', 'f', 'Patient')
    
%     P_inc_off = P_tot_off - P_coh_off - P_vc_off;
%     P_inc_on  = P_tot_on - P_coh_on - P_vc_on;
    P_inc_off(P_inc_off<=cutoff) = 0;
    P_inc_on(P_inc_on<=cutoff)   = 0;
    P_coh_off(P_coh_off<=cutoff) = 0;
    P_coh_on(P_coh_on<=cutoff)   = 0;
    P_vc_off(P_vc_off<=cutoff)   = 0;
    P_vc_on(P_vc_on<=cutoff)     = 0;
    

    % Define frequency bands
    clear fband
%     fband{1} = find(f>=1  & f<=4 );
    fband{1} = find(f>=4  & f<=7 );
    fband{2} = find(f>=7  & f<=13);
%     fband{3} = find(f>=13 & f<=30);
    fband{3} = find(f>=13 & f<=20);
    fband{4} = find(f>=20 & f<=30);
    fband{5} = find(f>=60 & f<=85);
%     fband{8} = find(f>=250 & f<=350);
    Nb = length(fband);
    

    % Normalized absolute differences of PSDs
    for i=1:length(Patient)
        tmp(i)=any(strcmp(Patient{i}, {'P01_L', 'P02_R', 'P04_L'}));
    end
    tmp = ~tmp;
    P_tot_off(:,:,tmp) = NaN;
    P_inc_off(:,:,tmp) = NaN;
    P_coh_off(:,:,tmp) = NaN;
    P_vc_off(:,:,tmp)  = NaN;
    P_bip_off(:,:,tmp) = NaN;
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
    clear dP_inc dP_coh dP_vc dP_bip dP_incoh_OFF dP_incoh_ON tmp
    for i = 1:Nb
        dP_inc(:,i) = squeeze( nanmean( yinc(fband{i},:)-xinc(fband{i},:), 1 ) ...
                                ./ nanmean(xinc(fband{i},:), 1 ) )*100;
        dP_coh(:,i) = squeeze( nanmean( ycoh(fband{i},:)-xcoh(fband{i},:), 1 ) ...
                                ./ nanmean(xcoh(fband{i},:), 1 ) )*100;
        dP_vc(:,i)  = squeeze( nanmean( yvc(fband{i},:) -xvc(fband{i},:),  1 ) ...
                                ./ nanmean(xvc(fband{i},:),  1 ) )*100;
        dP_bip(:,i) = squeeze( nanmean( ybip(fband{i},:)-xbip(fband{i},:), 1 ) ...
                                ./ nanmean(xbip(fband{i},:), 1 ) )*100;
        dP_incoh_OFF(:,i)= squeeze( nanmean( xinc(fband{i},:), 1 ) ...
                                ./ nanmean(xcoh(fband{i},:), 1 ) );
        dP_incoh_ON(:,i) = squeeze( nanmean( yinc(fband{i},:), 1 ) ...
                                ./ nanmean(ycoh(fband{i},:), 1 ) );
    end
    clear xinc yinc xcoh ycoh xvc yvc totx toty
    
    
    % Wilcoxon sign-rank test
    for i=1:Nb
        sig_inc(i) = signrank ( dP_inc(:,i) );
        sig_coh(i) = signrank ( dP_coh(:,i) );
        sig_vc(i)  = signrank ( dP_vc(:,i)  );
        sig_bip(i) = signrank ( dP_bip(:,i) );
    end


    % Monopolar inc
    subplot(3,4,(ii-1)*4+1)
    boxplot( dP_inc, 'labels', textlab, 'whisker', 50)
    set(gca,'TickLabelInterpreter','tex');
    ylim(limy), grid on
    ylabel({['{\bf\fontsize{14}' type '}']; '\DeltaP/P_{OFF} [%]'})
    if ii==1
        title('Incoherent', 'fontw', 'bold', 'fontsi', 14)
    end
    hold all
    for i=1:Nb
        if sig_inc(i)<0.01;        
            plot([i-0.07 i+0.07],0.9*limy([2 2]),'*k');
        elseif sig_inc(i)<0.05;
            plot(i,0.9*limy(2),'*k');
        end
    end
    % Monopolar coherent
    subplot(3,4,(ii-1)*4+2)
    boxplot( dP_coh, 'labels', textlab, 'whisker', 50)
    set(gca,'TickLabelInterpreter','tex');
    ylim(limy), grid on
    if ii==1
        title('Local Coherent', 'fontw', 'bold', 'fontsi', 14)
    end
    hold all
    for i=1:Nb
        if sig_coh(i)<0.01;        
            plot([i-0.07 i+0.07],0.9*limy([2 2]),'*k');
        elseif sig_coh(i)<0.05;
            plot(i,0.9*limy(2),'*k');
        end
    end
    % Monopolar volume conduction
    subplot(3,4,(ii-1)*4+3)
    boxplot( dP_vc, 'labels', textlab, 'whisker', 50)
    set(gca,'TickLabelInterpreter','tex');
    ylim(limy), grid on
    if ii==1
        title('Volume Conduction', 'fontw', 'bold', 'fontsi', 14)
    end
    hold all
    for i=1:Nb
        if sig_vc(i)<0.01;        
            plot([i-0.07 i+0.07],0.9*limy([2 2]),'*k');
        elseif sig_vc(i)<0.05;
            plot(i,0.9*limy(2),'*k');
        end
    end
    % Bipolar
    subplot(3,4,(ii-1)*4+4)
    boxplot( dP_bip, 'labels', textlab, 'whisker', 50)
    set(gca,'TickLabelInterpreter','tex');
    set(gca, 'YTick', [-200 -150 -100 -50 0 50 100 150 200]);
    ylim(limy2), grid on
    if ii==1
        title('Bipolar', 'fontw', 'bold', 'fontsi', 14)
    end
    hold all
    for i=1:Nb
        if sig_bip(i)<0.01;        
            plot([i-0.07 i+0.07],0.9*limy2([2 2]),'*k');
        elseif sig_bip(i)<0.05;
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