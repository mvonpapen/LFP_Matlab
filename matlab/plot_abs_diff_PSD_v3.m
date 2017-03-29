%%  Plot relative and absolute (percentual) ON-OFF difference of PSDs

nch = 4;
cutoff = 1e-12;
exp = {'rest', 'hold', 'fist'};
textlab = {'\beta', '\beta_1', '\beta_2', '\gamma'};
limy = [-50 50];
limy2= [-100 100];
yticks = [-50 -25 0 25 50];
fig1 = figure('PaperSize', [26 20], ...
        'PaperPositionmode', 'manual', 'PaperPosition', [1 1 24 18], ...
        'Visible', 'off');
fnam = 'psd_diff_test_ds10_pc15.eps';



%%%%%%%% vvv ABSOLUTE DIFFERENCES vvv %%%%%%%%%%%%%%%

for ii = 1:3
    
    type = exp{ii};
    load(['pow_coh_ds10_pc15_w0_12_' type '.mat'], ...
        'P_coh_off', 'P_coh_on', 'P_vc_off', 'P_vc_on', 'P_bip_off', ...
        'P_bip_on', 'P_inc_off', 'P_inc_on', 'f', 'Patient')
    
    val = 0;
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
    fband{1} = find(f>=13 & f<=30);
    fband{2} = find(f>=13 & f<=20);
    fband{3} = find(f>=20 & f<=30);
    fband{4} = find(f>=60 & f<=85);
    Nb = length(fband);
    

    % Normalized absolute differences of PSDs
%     for i=1:length(Patient)
%         tmp(i)=any(strcmp(Patient{i}, {'AUIN_L', 'BEMI_R', 'GRFR_L'}));
%     end
%     tmp = ~tmp;
%     P_tot_off(:,:,tmp) = NaN;
%     P_bip_off(:,:,tmp) = NaN;
    totx = squeeze(nanmean(P_tot_off,2));
    xinc = squeeze(nanmean(P_inc_off,2));
    yinc = squeeze(nanmean(P_inc_on,2));
    xcoh = squeeze(nanmean(P_coh_off,2));
    ycoh = squeeze(nanmean(P_coh_on,2));
    xvc  = squeeze(nanmean(P_vc_off,2));
    yvc  = squeeze(nanmean(P_vc_on,2));
    xbip = squeeze(nanmean(P_bip_off,2));
    ybip = squeeze(nanmean(P_bip_on,2));
    
    
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
    dPabs_inc(abs(dPabs_inc)==Inf) = 100;
    dPabs_coh(abs(dPabs_coh)==Inf) = 100;
    dPabs_vc(abs(dPabs_vc)==Inf) = 100;
    dPabs_bip(abs(dPabs_bip)==Inf) = 100;
    clear xinc yinc xcoh ycoh xvc yvc totx toty
    
    
    % Wilcoxon sign-rank test
    for i=1:Nb
        sig_inc(i) = signrank ( dPabs_inc(:,i) );
        sig_coh(i) = signrank ( dPabs_coh(:,i) );
        sig_vc(i)  = signrank ( dPabs_vc(:,i)  );
        sig_bip(i) = signrank ( dPabs_bip(:,i) );
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
        if sig_inc(i)<0.01;        
            plot([i-0.07 i+0.07],0.9*limy([2 2]),'*k');
        elseif sig_inc(i)<0.05;
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
    boxplot( dPabs_vc, 'labels', textlab, 'whisker', 50)
    set(gca,'TickLabelInterpreter','tex');
    set(gca, 'YTick', yticks);
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
    boxplot( dPabs_bip, 'labels', textlab, 'whisker', 50)
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

print(fig1, fnam, '-depsc')
close(fig1)