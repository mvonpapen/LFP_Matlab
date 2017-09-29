%%  Plot relative and absolute (percentual) ON-OFF difference of PSDs

cutoff = 1e-20;
exp = {'rest', 'hold', 'fist'};
textlab = {'\theta', '\alpha', '\beta_1', '\beta_2', '\gamma'};
limy = [-150 150];
limy2= [-150 150];
yticks = [-150 -100 -50 0 50 100 150];
fig1 = figure('PaperSize', [13 8], ...
        'PaperPositionmode', 'manual', 'PaperPosition', [1 1 12 6], ...
        'Visible', 'off');
fnam = 'psd_diff_bip_all_t-g_ds10.eps';
STNact = 0;


%%%%%%%% vvv ABSOLUTE DIFFERENCES vvv %%%%%%%%%%%%%%%

for ii = 1:3
    
    task = exp{ii};
    load(['pow_coh_bip-cohinc_' task '.mat'], 'P_bip_on', 'P_bip_off', ...
        'P_coh_off', 'P_coh_on', 'P_bip_off', 'P_inc_on', 'P_inc_off', ...
        'f', 'Patient', 'LFPact')
    
    P_inc_off(P_inc_off<=cutoff) = 0;
    P_inc_on(P_inc_on<=cutoff)   = 0;
    P_coh_off(P_coh_off<=cutoff) = 0;
    P_coh_on(P_coh_on<=cutoff)   = 0;
    
    % no bipolar PSD of nonSTN electrodes
    for i=1:length(LFPact)
        combo = nchoosek(1:length(LFPact{i}),2);
        STNsum = sum(LFPact{i}(combo),2);
        P_bip_off(:,STNsum<STNact,i) = NaN;
    end

    % Define frequency bands
    clear fband
    fband{1} = find(f>=4  & f<=7 );
    fband{2} = find(f>=7  & f<=13);
    fband{3} = find(f>=13 & f<=20);
    fband{4} = find(f>=20 & f<=30);
    fband{5} = find(f>=60 & f<=85);
    Nb = length(fband);
    

    % Normalized absolute differences of PSDs
%     for i=1:length(Patient)
%         tmp(i)=strcmp(Patient{i},'P01_L');
%     end
%     P_tot_off(:,:,tmp) = NaN;
%     P_bip_off(:,:,tmp) = NaN;
    totx = P_bip_off(:,:);
    toty = P_bip_on(:,:);
    xinc = P_inc_off(:,:);
    yinc = P_inc_on(:,:);
    xcoh = P_coh_off(:,:);
    ycoh = P_coh_on(:,:);
    clear dP_inc dP_coh dP_bip tmp
    for i = 1:Nb
        dP_inc(:,i) = squeeze( nanmean( yinc(fband{i},:)-xinc(fband{i},:), 1 ) ...
                                ./ nanmean(totx(fband{i},:), 1 ) )*100;
        dP_coh(:,i) = squeeze( nanmean( ycoh(fband{i},:)-xcoh(fband{i},:), 1 ) ...
                                ./ nanmean(totx(fband{i},:), 1 ) )*100;
        dP_bip(:,i) = squeeze( nanmean( toty(fband{i},:)-totx(fband{i},:), 1 ) ...
                                ./ nanmean(totx(fband{i},:), 1 ) )*100;
    end
    dP_inc(abs(dP_inc)==Inf) = NaN;
    dP_coh(abs(dP_coh)==Inf) = NaN;
    dP_bip(abs(dP_bip)==Inf) = NaN;
    clear xinc yinc xcoh ycoh totx toty
    
    
    % Wilcoxon sign-rank test
    for i=1:Nb
        sig_inc(i) = signrank ( dP_inc(:,i) );
        sig_coh(i) = signrank ( dP_coh(:,i) );
        sig_bip(i) = signrank ( dP_bip(:,i) );
    end


    % bipolar inc
    subplot(3,3,(ii-1)*3+1)
    boxplot( dP_inc, 'labels', textlab, 'whisker', 50)
    set(gca,'TickLabelInterpreter','tex');
    set(gca, 'YTick', yticks);
    ylim(limy), grid on
    ylabel({['{\bf\fontsize{14}' task '}']; '\DeltaP/P_{OFF} [%]'})
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
    % bipolar coherent
    subplot(3,3,(ii-1)*3+2)
    boxplot( dP_coh, 'labels', textlab, 'whisker', 50)
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
    % Bipolar
    subplot(3,3,(ii-1)*3+3)
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
print(fig1, fnam, '-depsc')
close(fig1)