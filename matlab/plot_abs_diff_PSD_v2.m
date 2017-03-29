%%  Plot relative and absolute (percentual) ON-OFF difference of PSDs

nch = 4;
cutoff = 1e-20;
exp = {'Rest', 'Hold', 'Fist'};
textlab = {'Inc', 'Coh', 'VC', 'Bip'};
limy = [-80 80];
fig1 = figure('Papersize', [10 3], 'PaperPosition', [0.5 0.5 9 2], ...
        'PaperPositionmode', 'manual', 'Visible', 'off');
onlySTN = 0;

%%%%%%%% vvv ABSOLUTE DIFFERENCES vvv %%%%%%%%%%%%%%%

for i = 1:3
    
    task = exp{i};
    load(['pow_coh_ds10_pc15_w0_12_' task '.mat'], 'P_tot_off', 'P_tot_on', ...
        'P_coh_off', 'P_coh_on', 'P_vc_off', 'P_vc_on', 'P_bip_off', ...
        'P_bip_on', 'P_inc_off', 'P_inc_on', 'f', 'LFPact')
    
%     P_inc_off=P_tot_off-P_coh_off-P_vc_off;
%     P_inc_on=P_tot_on-P_coh_on-P_vc_on;
    P_inc_off(P_inc_off<=cutoff) = 0;
    P_inc_on(P_inc_on<=cutoff)   = 0;
    P_coh_off(P_coh_off<=cutoff) = 0;
    P_coh_on(P_coh_on<=cutoff)   = 0;
    P_vc_off(P_vc_off<=cutoff)   = 0;
    P_vc_on(P_vc_on<=cutoff)     = 0;
    
    % Only use channels with STN activity
    if onlySTN == 1
        for j=1:length(LFPact)
            P_tot_off(:,LFPact{j}<1,j) = NaN;
            P_inc_off(:,LFPact{j}<1,j) = NaN;
            P_coh_off(:,LFPact{j}<1,j) = NaN;
            P_vc_off(:,LFPact{j}<1,j) = NaN;
        end
    end

    % Define frequency bands
    fband = find(f>=13 & f<=20);
%     fband{1} = find(f>=1 & f<=4);
%     fband{2} = find(f>=4 & f<=7);
    

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
    

    
    %% Calculate difference in power spectra OFF->ON
    
    
%     dPabs_tot{i}  = ... squeeze( nanmean( yinc(fband,:)-xinc(fband,:), 1 ) ...
%                         ...    ./ nanmean(totx(fband,:), 1 ) )*100;
    dPabs_inc{i}  = squeeze( nanmean( yinc(fband,:)-xinc(fband,:), 1 ) ...
                            ./ nanmean(totx(fband,:), 1 ) )*100;
    dPabs_coh{i} = squeeze( nanmean( ycoh(fband,:)-xcoh(fband,:), 1 ) ...
                            ./ nanmean(totx(fband,:), 1 ) )*100;
    dPabs_vc{i} = squeeze( nanmean( yvc(fband,:) -xvc(fband,:),  1 ) ...
                            ./ nanmean(totx(fband,:), 1 ) )*100;
    dPabs_acoh{i} = squeeze( nanmean(ycoh(fband,:)+yvc(fband,:) ...
        -xcoh(fband,:)-xvc(fband,:),1) ./ nanmean(totx(fband,:), 1 ) )*100;
    dPabs_nvc{i}  = squeeze( nanmean(ycoh(fband,:)+yinc(fband,:) ...
        -xcoh(fband,:)-xinc(fband,:),1) ./ nanmean(totx(fband,:), 1 ) )*100;
    dPabs_bip{i}  = squeeze( nanmean( ybip(fband,:)-xbip(fband,:), 1 ) ...
                            ./ nanmean(xbip(fband,:), 1 ) )*100;
    clear xinc yinc xcoh ycoh xvc yvc totx toty
    
    
    % Wilcoxon sign-rank test
    sig(1,i)  = signrank ( dPabs_inc{i}  );
    sig(2,i)  = signrank ( dPabs_coh{i} );
    sig(3,i)  = signrank ( dPabs_vc{i} );
%     sig(4,i)  = signrank ( dPabs_nvc{i}  );
    sig(4,i)  = signrank ( dPabs_bip{i}  );
%     sig(6,i)  = signrank ( dPabs_acoh{i} );
    
end



%% Group results
for i=1:3
    X{i} = [dPabs_inc{i} dPabs_coh{i} dPabs_vc{i} dPabs_bip{i}];
    G{i} = zeros(1,length(X{i}));
    n  = length(dPabs_inc{i});
    G{i}(1:n) = 1;
    G{i}(n+1:2*n) = 2;
    G{i}(2*n+1:3*n) = 3;
%     G{i}(3*n+1:4*n) = 4;
    G{i}(3*n+1:end) = 4;

    % Plot
    h=subplot(1,3,i);
    boxplot( X{i}, G{i}, 'labels', textlab, 'whisker', 100)
    ylim(limy), grid on
    if i==1
        ylabel('\DeltaP/P_{OFF} [%]')
    end
    title(exp{i}, 'fontw', 'bold', 'fontsi', 14)
    hold all
    for ii=1:4
        if sig(ii,i)<0.01;        
            plot([ii-0.07 ii+0.07],0.9*limy([2 2]),'*k');
        elseif sig(ii,i)<0.05;
            plot(ii,0.9*limy(2),'*k');
        end
    end
%     a2 = axes('YAxisLocation', 'Right');
%     set(a2, 'color', 'none')
%     set(a2, 'XTick', [])
%     set(a2, 'YLim', 2*limy)

    
    
end
% set(fig1, 'visi', 'on')
% text( -3.8, 0.5, 'Fist','HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'middle', 'Rotation', 90, ...
%     'units', 'normalized', 'Fonts', 14, 'Fontw', 'bold')
% text( -3.8, 1.9, 'Hold','HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'middle', 'Rotation', 90, ...
%     'units', 'normalized', 'Fonts', 14, 'Fontw', 'bold')
% text( -3.8, 3.3, 'Rest','HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'middle', 'Rotation', 90, ...
%     'units', 'normalized', 'Fonts', 14, 'Fontw', 'bold')
% set(fig1, 'visi', 'on')
% saveas(fig1, 'psd_diff_abs_beta_t6w12.pdf')
saveas(fig1, 'psd_diff_abs_beta_ds10_pc15_w12.eps', 'epsc')
close(fig1)