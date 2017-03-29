%%  Plot absolute (percentual) ON-OFF difference of PSDs for rest and
%%  suppression/activation for hold and fist

clear all

nch = 4;
cutoff = 1e-20;
exp = {'Hold', 'Fist'};
textlab = {'\alpha', '\beta1', '\beta2'};
% textlab = {'\delta', '\theta', '\alpha', '\beta', '\beta1', '\beta2', ...
%         '\gamma', 'HF'};
limy = [-80 80];
limy2= [-120 120];
fig1 = figure('Papertype', 'A4', ...
        'Paperunits', 'normalized', 'PaperPosition', [0 0 1 1], ...
        'PaperPositionmode', 'manual', 'Visible', 'off', ...
        'PaperOrientation', 'landscape');

load(['pow_coh_t6_w0_12_Rest.mat'])
totx = P_coh_off;
toty = P_coh_on;
bipx = P_bip_off;
bipy = P_bip_on;
ind{1} = [1 2 3 4 5 7 8];
ind{2} = [1 2 4 5 7 8];

%%%%%%%% vvv ABSOLUTE DIFFERENCES vvv %%%%%%%%%%%%%%%

for ii = 1:2
    
    type = exp{ii};
    load(['pow_coh_t6_w0_12_' type '.mat'])
    
    dP_inc_off = (P_tot_off-P_coh_off-P_vc_off) - totx(:,:,ind{ii});
    dP_inc_on  = (P_tot_on -P_coh_on -P_vc_on)  - toty(:,:,ind{ii});   
    dP_tot_off = P_tot_off - totx(:,:,ind{ii});
    dP_tot_on  = P_tot_on  - toty(:,:,ind{ii});
    dP_coh_off = (P_coh_off - totx(:,:,ind{ii})) ./ totx(:,:,ind{ii});
    dP_coh_on  = (P_coh_on - toty(:,:,ind{ii})) ./ toty(:,:,ind{ii});
    dP_vc_off  = P_vc_off  - totx(:,:,ind{ii});
    dP_vc_on   = P_vc_on   - toty(:,:,ind{ii});
    dP_bip_off = P_bip_off - bipx(:,:,ind{ii});
    dP_bip_on  = P_bip_on  - bipy(:,:,ind{ii});
    
%     dP_inc_off = (P_tot_off-P_coh_off-P_vc_off) ./ totx(:,:,ind{ii});
%     dP_inc_on  = (P_tot_on -P_coh_on -P_vc_on)  ./ toty(:,:,ind{ii});   
%     dP_tot_off = P_tot_off ./ totx(:,:,ind{ii});
%     dP_tot_on  = P_tot_on  ./ toty(:,:,ind{ii});
%     dP_coh_off = P_coh_off ./ totx(:,:,ind{ii});
%     dP_coh_on  = P_coh_on  ./ toty(:,:,ind{ii});
%     dP_vc_off  = P_vc_off  ./ totx(:,:,ind{ii});
%     dP_vc_on   = P_vc_on   ./ toty(:,:,ind{ii});
%     dP_bip_off = P_bip_off ./ bipx(:,:,ind{ii});
%     dP_bip_on  = P_bip_on  ./ bipy(:,:,ind{ii});
    

    % Define frequency bands
    clear fband
%     fband{1} = find(f>=1  & f<=4 );
%     fband{2} = find(f>=4  & f<=7 );
    fband{1} = find(f>=7  & f<=13);
%     fband{4} = find(f>=13 & f<=30);
    fband{2} = find(f>=13 & f<=20);
    fband{3} = find(f>=20 & f<=30);
%     fband{7} = find(f>=60 & f<=90);
%     fband{8} = find(f>=250 & f<=350);
    Nb = length(fband);
    

    % Normalized absolute differences of PSDs
    xinc = dP_inc_off(:,:);
    yinc = dP_inc_on(:,:);
    xcoh = dP_coh_off(:,:);
    ycoh = dP_coh_on(:,:);
    xvc  = dP_vc_off(:,:);
    yvc  = dP_vc_on(:,:);
    xbip = dP_bip_off(:,:);
    ybip = dP_bip_on(:,:);
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
    clear dPabs_inc dPabs_coh dPabs_vc dPabs_bip
    for i = 1:Nb
        dPabs_inc(:,i) = squeeze( nanmean( yinc(fband{i},:)-xinc(fband{i},:), 1 ) ...
                                ./ nanmean(totx(fband{i},:), 1 ) )*100;
        dPabs_coh(:,i) = squeeze( nanmean( ycoh(fband{i},:)-xcoh(fband{i},:), 1 ) ...
                                ./ nanmean(totx(fband{i},:), 1 ) )*100;
        dPabs_vc(:,i)  = squeeze( nanmean( yvc(fband{i},:) -xvc(fband{i},:),  1 ) ...
                                ./ nanmean(totx(fband{i},:), 1 ) )*100;
        dPabs_bip(:,i) = squeeze( nanmean( ybip(fband{i},:)-xbip(fband{i},:), 1 ) ...
                                ./ nanmean(xbip(fband{i},:), 1 ) )*100;
    end
    clear xinc yinc xcoh ycoh xvc yvc totx toty
    
    
    % Wilcoxon sign-rank test
    for i=1:Nb
        sig_inc(i) = signrank ( dPabs_inc(:,i) );
        sig_coh(i) = signrank ( dPabs_coh(:,i) );
        sig_vc(i)  = signrank ( dPabs_vc(:,i)  );
        sig_bip(i) = signrank ( dPabs_bip(:,i) );
    end


    % Monopolar inc
    subplot(3,3,(ii-1)*3+1)
    boxplot( dPabs_inc, ...
        'labels', textlab)
    h = findobj(gca, 'type', 'text');
    set(h, 'Interpreter', 'tex');
    ylim(limy), grid on
    ylabel('\DeltaP/P_{OFF} [%]')
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
    subplot(3,3,(ii-1)*3+2)
    boxplot( dPabs_coh, ...
        'labels', textlab)
    h = findobj(gca, 'type', 'text');
    set(h, 'Interpreter', 'tex');
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
    subplot(3,3,(ii-1)*3+3)
    boxplot( dPabs_vc, ...
        'labels', textlab)
    h = findobj(gca, 'type', 'text');
    set(h, 'Interpreter', 'tex');
    ylim(limy), grid on
    if ii==1
        title('Remote Coherent', 'fontw', 'bold', 'fontsi', 14)
    end
    hold all
    for i=1:Nb
        if sig_vc(i)<0.01;        
            plot([i-0.07 i+0.07],0.9*limy([2 2]),'*k');
        elseif sig_vc(i)<0.05;
            plot(i,0.9*limy(2),'*k');
        end
    end
%     % Bipolar
%     subplot(3,4,(ii-1)*4+4)
%     boxplot( dPabs_bip, ...
%         'labels', textlab)
%     h = findobj(gca, 'type', 'text');
%     set(h, 'Interpreter', 'tex');
%     ylim(limy2), grid on
%     if ii==1
%         title('Bipolar', 'fontw', 'bold', 'fontsi', 14)
%     end
%     hold all
%     for i=1:Nb
%         if sig_bip(i)<0.01;        
%             plot([i-0.07 i+0.07],0.9*limy2([2 2]),'*k');
%         elseif sig_bip(i)<0.05;
%             plot(i,0.9*limy2(2),'*k');
%         end
%     end

end

text( -3.8, 0.5, 'Fist','HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'Rotation', 90, ...
    'units', 'normalized', 'Fonts', 14, 'Fontw', 'bold')
text( -3.8, 1.9, 'Hold','HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'Rotation', 90, ...
    'units', 'normalized', 'Fonts', 14, 'Fontw', 'bold')
text( -3.8, 3.3, 'Rest','HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'Rotation', 90, ...
    'units', 'normalized', 'Fonts', 14, 'Fontw', 'bold')
set(fig1, 'visi', 'on')
% saveas(fig1, 'psd_diff_abs_all.pdf')
% saveas(fig1, 'psd_diff_abs_all.eps', 'epsc')
% close(fig1)