Patients = {'P02_R', 'P03_L', 'P04_L', 'P06_L', ...
    'P07_L', 'P08_L', 'P01_L', 'P05_L'};

cutoff = 1e-10;

for patnum = 1:length(Patients);
    pat    = Patients{patnum};

    load(['pow_coh_' pat '.mat'])
    STNact{patnum} = LFPact{1};
    Nelec(patnum)  = length(LFPact{1});
    rig{patnum} = rigor';
    Nrest(patnum) = length(rigor);
    P_inc(P_inc<=cutoff) = cutoff;
    P_coh(P_coh<=cutoff) = cutoff;
    P_vc(P_vc<=cutoff)   = cutoff;
    P_nvc(P_nvc<=cutoff) = cutoff;
    
    % Define frequency bands
%     fband{1} = find(f>=1  & f<=4 );
%     fband{2} = find(f>=4  & f<=7 );
%     fband{3} = find(f>=7  & f<=13);
%     fband{4} = find(f>=13 & f<=30);
    fband{1} = find(f>=13 & f<=20);
    fband{2} = find(f>=20 & f<=30);
%     fband{7} = find(f>=30 & f<=45);
    fband{3} = find(f>=60 & f<=85);

    
    
    %% Electrode power
    
%     % Average over electrodes
%     totx = squeeze(nanmean(P_tot,2));
%     xinc = squeeze(nanmean(P_inc,2));
%     xcoh = squeeze(nanmean(P_coh,2));
%     xvc  = squeeze(nanmean(P_vc, 2));
%     xnvc = squeeze(nanmean(P_nvc,2));
%     xaco = squeeze(nanmean(P_vc+P_coh,2));
%     xbip = squeeze(nanmean(P_bip,2));
%     Pc{patnum} = P_coh;  
%     Pi{patnum} = P_inc;  
%     Pv{patnum} = P_vc;  
%     Pb{patnum} = P_bip;
    
%     % Patient specific selection
%     if any(patnum==[1 3 7])
%         totx = squeeze(nanmean(P_tot(:,:,2:end),2));
%         xinc = squeeze(nanmean(P_inc(:,:,2:end),2));
%         xcoh = squeeze(nanmean(P_coh(:,:,2:end),2));
%         xvc  = squeeze(nanmean(P_vc(:,:,2:end), 2));
%         xnvc = squeeze(nanmean(P_nvc(:,:,2:end),2));
%         xaco = squeeze(nanmean(P_vc(:,:,2:end)+P_coh(:,:,2:end),2));
%         xbip = squeeze(nanmean(P_bip(:,:,2:end),2));
%         rig{patnum} = rig{patnum}(2:end);
%     end
        
    
    % All electrodes
    totx = P_tot(:,:);
    xinc = P_inc(:,:);
    xcoh = P_coh(:,:);
    xvc  = P_vc(:,:);
%     xaco = P_coh(:,:)+P_vc(:,:);
%     xnvc = P_nvc(:,:);
    xbip = P_bip(:,:);
    tmp          = repmat( rigor, 6, 1 );
    rig2{patnum} = tmp(:);
    tmp          = repmat( rigor, Nelec(patnum), 1 );
    rig{patnum}  = tmp(:);
    clear tmp
    
%     % Specified electrodes
% %     E = [3 3 2 1 3 4 1 1]; % for total
% %     E = [2 3 3 1 1 4 1 1]; % for incoherent
%     E = [1 2 3 2 1 1 1 1]; % for coh signal
% %     E = [1 2 2 1 1 2 1 1]; % for volume conduction
% %     E2 = [2 6 2 2 2 6 1 1]; % for bipolar
% %     E = [1 1 1 1 1 1 0 1]+2;
% %     E(patnum) = Nelec(patnum);
%     totx = squeeze(P_tot(:,E(patnum),:));
%     xinc = squeeze(P_inc(:,E(patnum),:));
%     xcoh = squeeze(P_coh(:,E(patnum),:));
%     xvc  = squeeze(P_vc(:, E(patnum),:));
%     xnvc = squeeze(P_nvc(:,E(patnum),:));
%     xbip = squeeze(P_bip(:,E(patnum),:));

    
    

    Nb = length(fband);
%     i0 = ones(1,length(rig{patnum}));
    i0 = repmat(1:Nelec(patnum), 1, Nrest(patnum));
    i2 = repmat(1:6,             1, Nrest(patnum));
%     if any(patnum==[1 3 7])
%         i0 = i0+1;
%     end
    for i = 1:Nb
% %         dP_tot{patnum}(:,i) = squeeze( nanmean( totx(fband{i},:)-totx(fband{i},i0), 1 ) ...
% %             ./ nanmean(totx(fband{i},i0), 1 ) )*100;
        dP_inc{patnum}(:,i) = squeeze( nanmean( xinc(fband{i},:)-xinc(fband{i},i0), 1 ) ...
            ./ nanmean(totx(fband{i},i0), 1 ) )*100;
        dP_coh{patnum}(:,i) = squeeze( nanmean( xcoh(fband{i},:)-xcoh(fband{i},i0), 1 ) ...
            ./ nanmean(totx(fband{i},i0), 1 ) )*100;
        dP_vc{patnum}(:,i)  = squeeze( nanmean( xvc(fband{i},:) - xvc(fband{i},i0), 1 ) ...
            ./ nanmean(totx(fband{i},i0), 1 ) )*100;
% %         dP_nvc{patnum}(:,i) = squeeze( nanmean( xnvc(fband{i},:)-xnvc(fband{i},i0), 1 ) ...
% %             ./ nanmean(totx(fband{i},i0), 1 ) )*100;
% %         dP_aco{patnum}(:,i) = squeeze( nanmean( xaco(fband{i},:)-xaco(fband{i},i0), 1 ) ...
% %             ./ nanmean(totx(fband{i},i0), 1 ) )*100;
%         dP_bip{patnum}(:,i) = squeeze( nanmean( xbip(fband{i},:)-xbip(fband{i},i0), 1 ) ...
%             ./ nanmean(xbip(fband{i},i0), 1 ) )*100;
        dP_bip{patnum}(:,i) = squeeze( nanmean( xbip(fband{i},:)-xbip(fband{i},i2), 1 ) ...
            ./ nanmean(xbip(fband{i},i2), 1 ) )*100;
    end


    for i=1:Nb %loop over freq.band
% %         [rho_tot(i) p_tot(i)] = corr(rig{patnum}, dP_tot{patnum}(:,i), 'type', 'Pearson');
        [rho_inc(i) p_inc(i)] = corr(rig{patnum}, dP_inc{patnum}(:,i), 'type', 'Pearson');
        [rho_coh(i) p_coh(i)] = corr(rig{patnum}, dP_coh{patnum}(:,i), 'type', 'Pearson');
        [rho_vc(i)  p_vc(i)]  = corr(rig{patnum}, dP_vc{patnum}(:,i),  'type', 'Pearson');
% %         [rho_nvc(i) p_nvc(i)] = corr(rig{patnum}, dP_nvc{patnum}(:,i), 'type', 'Pearson');
%         [rho_bip(i) p_bip(i)] = corr(rig{patnum}, dP_bip{patnum}(:,i), 'type', 'Pearson');
        j = ~isnan(dP_bip{patnum}(:,i));
        [rho_bip(i) p_bip(i)] = corr(rig2{patnum}(j), dP_bip{patnum}(j,i), 'type', 'Pearson');
    end
    Nbib(patnum) = sum(j(1:6));

%     C(:,:,patnum) = [rho_tot; rho_inc; rho_coh; rho_vc; rho_nvc; rho_bip];
%     P(:,:,patnum) = [p_tot;   p_inc;   p_coh;   p_vc;   p_nvc;   p_bip];
    C(:,:,patnum) = [rho_inc; rho_coh; rho_vc; rho_bip];
    P(:,:,patnum) = [p_inc;   p_coh;   p_vc;   p_bip];


    % %% Plot results for one freqeuncy band
    % fig1 = figure('Papersize', [6.5 5], 'PaperPosition', [0.5 0.5 5.5 4], ...
    %         'PaperPositionmode', 'manual', 'Visible', 'off');
    % fb = [5 8]; %freq. bands to plot, 5=lowbeta, 7=highgamma
    % for n=1:2
    %     subplot(1,2,n)    
    %     i    = fb(n);
    %     limx = [0 110];
    %     limy = [-20 5];
    % 
    %     plot(rigor(:), dP_tot(:,i), 'ok', rigor(:), dP_inc(:,i), 'vk', ...
    %         rigor(:), dP_coh(:,i), '^k', rigor(:), dP_vc(:,i),  'xk', ...
    %         rigor(:), dP_nvc(:,i), 'sk', rigor(:), dP_bip(:,i)/4, '*k')
    %     xlim(limx), ylim(limy)
    %     legend(['total, c=' mat2str(rho_tot(i),2)], ...
    %         ['inc., c=' mat2str(rho_inc(i),2)], ...
    %         ['coh., c=' mat2str(rho_coh(i),2)], ...
    %         ['vol.con.,  c=' mat2str(rho_vc(i),2)], ...
    %         ['local, c=' mat2str(rho_nvc(i),2)], ...
    %         ['bipolar, c=' mat2str(rho_bip(i),2)], ...
    %         'location', 'southwest')
    %     xlabel('Rigidity improvement over time [%]')
    %     ylabel('Average \DeltaP/P_{OFF} [%] per patient')
    %     a = polyfit(rigor(:)',dP_coh(:,i),1);
    % 
    %     hold all
    %     plot([15 105], a(1)*[15 105]+a(2), '--k')
    % %     a2 = axes('YAxisLocation', 'Right');
    % %     set(a2, 'color', 'none', 'ycolor', 'r', 'XTick', [], 'YLim', limy*4)
    % %     ylabel('Average \DeltaP/P [%] for bipolar')
    %     grid on
    % end
    % 
    % set(fig1, 'visi', 'on');
    % % saveas(fig1, ['corr_beta_' pat], 'epsc')
    % % close(fig1)

end

% % fprintf('Correlation of rigidity with power decrease, coefficients and p-values for [tot; inc; coh; vc; nvc; bip]')
% 
% j=3; 'loc',
% [c,b]=min(squeeze(C(j,4,1:7)));
% tmp=squeeze(P(j,4:6,1:7));
% for i=1:7
%     p(i)  = tmp(b(i),i);
%     b0(i) = dP_coh{i}(1,3+b(i));
%     b1(i) = dP_coh{i}(end,3+b(i));
% end
% [c mean(c) std(c); p mean(p) std(p); b NaN NaN; ...
%     Nrest(1:7) mean(Nrest(1:7)) std(Nrest(1:7)); ...
%     Nelec(1:7) mean(Nelec(1:7)) std(Nelec(1:7)); ...
%     b0 mean(b0) std(b0); b1-b0 mean(b1-b0) std(b1-b0)]
% 
% 

%% Calculate correlation
R=[]; X=[]; fb=2;
load time_of_rest_measurements

% % for patient-averages
% for i=[2 4 5 6 8];
%     X = [X; dP_bip{i}]; 
%     R = [R T{i}]; 
% end
% j = ~isnan(X(:,fb));
% [rho, p] = corr(R(j)',X(j,fb)),

% for all channels
% for i=[2 4 5 6 8];
for i=1:8;
    for j=1:Nbib(i)
        X = [X; dP_bip{i}(j:6:end,:)];
        R = [R T{i}];
    end
end
j = ~isnan(X(:,fb));
[rho(fb), p(fb)] = corr(R(j)',X(j,fb));
[a(fb,:), S{fb}] = polyfit(R(j)',X(j,fb),1);

%% Plot results
figure
pn = [7 1 2 3 8 4 5 6];
subplot(1,3,1)
% % for patient-averages
% for i=pn
%     plot(T{i}, dP_bip{i}(:,fb), '.-')
%     hold all
% end
% legend('Pat. 1', 'Pat 2', 'Pat 3', 'Pat 4', 'Pat 5', 'Pat 6', 'Pat 7', 'Pat 8', ...
%     'Location', 'SouthEast')
% for all electrodes
for i=pn
    for j=1:Nbib(i)
        plot(T{i}, dP_bip{i}(j:6:end,fb), '.-')
        hold all
    end
end
%
xlim([-2 12])
title('Change in low-beta power over time')
xlabel('mins after apomorphine injection')
ylabel('\DeltaP/P_{OFF}')
% Beta peak
subplot(1,3,2),
[p1, p2, ~, h1, h2] = avg_spectra ( 'pow_coh_ds10_pc15_w0_12_rest', ...
    'P_bip_off', 'P_bip_on' , 'pat', [2 4 5 6 8], 'log', 1); 
xlim([1 45]), %ylim([1e-8 1e-5])
ylabel('$\left<P\right> / (\overline{P/P_\mu})$', 'Interpreter', 'Latex')
xlabel('f [Hz]')
title('PSD for patients (3, 5-8)')
legend([h1.patch h2.patch], 'OFF', 'ON')
% No beta peak
subplot(1,3,3),
[~, ~, ~, h1, h2] = avg_spectra ( 'pow_coh_ds10_pc15_w0_12_rest', ...
    'P_bip_off', 'P_bip_on' , 'pat', [1 3 7], 'log', 1); 
xlim([1 45]), %ylim([6e-9 1e-6])
ylabel('$\left<P\right> / (\overline{P/P_\mu})$', 'Interpreter', 'Latex')
xlabel('f [Hz]')
legend([h1.patch h2.patch], 'OFF', 'ON')
title('PSD for patients (1-2, 4)')