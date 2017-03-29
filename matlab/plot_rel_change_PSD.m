%%  Plot relative and absolute (percentual) ON-OFF difference of PSDs

nch = 4;
cutoff = 1e-10;
exp = {'rest', 'hold', 'fist'};
textlab = {'\beta_1', '\beta_2', '\gamma'};
limy = [-2 2];
limy2= [-2 2];
yticks = [-2 -1 0 1 2];
fig1 = figure('PaperSize', [14 8], ...
        'PaperPositionmode', 'manual', 'PaperPosition', [1 1 12 6], ...
        'Visible', 'off');
tag  = 'ds10_pc15';
fnam = ['psd_relchg_' tag '.eps'];
ncomp= 3;



%%%%%%%% vvv ABSOLUTE DIFFERENCES vvv %%%%%%%%%%%%%%%

for ii = 1:3
    
    type = exp{ii};
    load(['pow_coh_' tag '_w0_12_' type '.mat'], 'P_tot_off', 'P_tot_on', ...
        'P_coh_off', 'P_coh_on', 'P_vc_off', 'P_vc_on', 'P_bip_off', ...
        'P_bip_on', 'P_inc_off', 'P_inc_on', 'f', 'Patient')
    
    val = cutoff;
    P_inc_off(P_inc_off<=cutoff) = val;
    P_inc_on(P_inc_on<=cutoff)   = val;
    P_coh_off(P_coh_off<=cutoff) = val;
    P_coh_on(P_coh_on<=cutoff)   = val;
    P_vc_off(P_vc_off<=cutoff)   = val;
    P_vc_on(P_vc_on<=cutoff)     = val;
    

    % Define frequency bands
    clear fband
    fband{1} = find(f>=13 & f<=20);
    fband{2} = find(f>=20 & f<=30);
    fband{3} = find(f>=60 & f<=85);
    Nb = length(fband);
    

    
    
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
    
    clear inc coh vc bip tmp
    for i = 1:Nb
        inc(:,i) = log10( squeeze( nanmean( yinc(fband{i},:)./xinc(fband{i},:), 1 ) ) );
        coh(:,i) = log10( squeeze( nanmean( ycoh(fband{i},:)./xcoh(fband{i},:), 1 ) ) );
        vc(:,i)  = log10( squeeze( nanmean( yvc(fband{i},:) ./xvc(fband{i},:),  1 ) ) );
        bip(:,i) = log10( squeeze( nanmean( ybip(fband{i},:)./xbip(fband{i},:), 1 ) ) );
    end
    inc(abs(inc)==Inf) = 10;
    coh(abs(coh)==Inf) = 10;
    vc(abs(vc)==Inf) = 10;
    bip(abs(bip)==Inf) = 10;
    clear xinc yinc xcoh ycoh xvc yvc totx toty
    
    
    % Wilcoxon sign-rank test
    for i=1:Nb
        sig_inc(i,ii) = signrank ( inc(:,i) ) * ncomp;
        sig_coh(i,ii) = signrank ( coh(:,i) ) * ncomp;
        sig_vc(i,ii)  = signrank ( vc(:,i)  ) * ncomp;
        sig_bip(i,ii) = signrank ( bip(:,i) ) * ncomp;
    end

    % Monopolar inc
    subplot(3,4,(ii-1)*4+1)
    boxplot( inc, 'labels', textlab, 'whisker', 50)
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
    boxplot( coh, 'labels', textlab, 'whisker', 50)
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
    boxplot( vc, 'labels', textlab, 'whisker', 50)
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
    boxplot( bip, 'labels', textlab, 'whisker', 50)
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

print(fig1, fnam, '-depsc')
close(fig1)