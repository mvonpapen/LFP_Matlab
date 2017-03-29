clear dP_inc dP_coh dP_vc dP_bip fband STNact Nelec rig Nrest rig2 Nbib


Patients = {'BEMI_R', 'BIMA_L', 'GRFR_L', 'REGE_L', ...
    'RIRO_L', 'SUGI_L', 'AUIN_L', 'MEGE_L'};

cutoff = 1e-10;

for patnum = 1:length(Patients);
    pat    = Patients{patnum};

    load(['pow_coh_' pat '.mat'])
    STNact{patnum} = LFPact{1};
    combo = nchoosek(1:length(LFPact{1}),2);
    STNelec{patnum}(1:size(combo,1)) = sum(LFPact{1}(combo)>0,2);
    Nelec(patnum) = length(LFPact{1});
    rig{patnum}   = rigor';
    Nrest(patnum) = length(rigor);
    P_inc(P_inc<=cutoff) = cutoff;
    P_coh(P_coh<=cutoff) = cutoff;
    P_vc(P_vc<=cutoff)   = cutoff;
    P_nvc(P_nvc<=cutoff) = cutoff;
    
    % Define frequency bands
    fband{1} = find(f>=13 & f<=20);
    fband{2} = find(f>=20 & f<=30);
    fband{3} = find(f>=30 & f<=40);
    fband{4} = find(f>=60 & f<=90);
    LabFb    = {'low \beta', 'high \beta', 'low \gamma', 'high \gamma'};
    Nb = length(fband);

    
    
    %% Electrode power 
    
    % All electrodes
    totx = P_tot(:,:);
    xinc = P_inc(:,:);
    xcoh = P_coh(:,:);
    xvc  = P_vc(:,:);
    xbip = P_bip(:,:);
    tmp          = repmat( rigor, 6, 1 );
    rig2{patnum} = tmp(:);
    tmp          = repmat( rigor, Nelec(patnum), 1 );
    rig{patnum}  = tmp(:);
    clear tmp
    

    i0 = repmat(1:Nelec(patnum), 1, Nrest(patnum));
    i2 = repmat(1:6,             1, Nrest(patnum));
    for i = 1:Nb
        dP_inc{patnum}(:,i) = squeeze( nanmean( xinc(fband{i},:)-xinc(fband{i},i0), 1 ) ...
            ./ nanmean(totx(fband{i},i0), 1 ) )*100;
        dP_coh{patnum}(:,i) = squeeze( nanmean( xcoh(fband{i},:)-xcoh(fband{i},i0), 1 ) ...
            ./ nanmean(totx(fband{i},i0), 1 ) )*100;
        dP_vc{patnum}(:,i)  = squeeze( nanmean( xvc(fband{i},:) - xvc(fband{i},i0), 1 ) ...
            ./ nanmean(totx(fband{i},i0), 1 ) )*100;
        dP_aco{patnum}(:,i)  = squeeze( nanmean( xvc(fband{i},:) - xvc(fband{i},i0) ...
            + xcoh(fband{i},:)-xcoh(fband{i},i0), 1 ) ./ nanmean(totx(fband{i},i0), 1 ) )*100;
        dP_bip{patnum}(:,i) = squeeze( nanmean( xbip(fband{i},:)-xbip(fband{i},i2), 1 ) ...
            ./ nanmean(xbip(fband{i},i2), 1 ) )*100;
    end

    clear rho_inc p_inc rho_coh p_coh rho_vc p_vc rho_bip p_bip
    for i=1:Nb %loop over freq.band
        [rho_inc(i), p_inc(i)] = corr(rig{patnum}, dP_inc{patnum}(:,i), 'type', 'Pearson');
        [rho_coh(i), p_coh(i)] = corr(rig{patnum}, dP_coh{patnum}(:,i), 'type', 'Pearson');
        [rho_aco(i), p_aco(i)] = corr(rig{patnum}, dP_aco{patnum}(:,i), 'type', 'Pearson');
        [rho_vc(i),  p_vc(i)]  = corr(rig{patnum}, dP_vc{patnum}(:,i),  'type', 'Pearson');
        j = ~isnan(dP_bip{patnum}(:,i));
        [rho_bip(i), p_bip(i)] = corr(rig2{patnum}(j), dP_bip{patnum}(j,i), 'type', 'Pearson');
    end
    Nbip(patnum) = sum(j(1:6));

    C(:,:,patnum) = [rho_inc; rho_coh; rho_vc; rho_bip];
    P(:,:,patnum) = [p_inc;   p_coh;   p_vc;   p_bip];

end


%% Calculate correlation
% % PCC
% Signal = dP_inc; %dP_inc/coh/vc
% Nch    = Nelec;
% Rig    = rig;
% n      = Nelec;
% Bipolar
Signal = dP_bip;
Nch    = Nbip;
Rig    = rig2;
n      = 6*ones(1,length(Nbip));


R=[]; Y=[]; R2=[]; R3=[];
load time_of_rest_measurements
x = [-2:0.2:12];

for i=1:8;
    for j=1:Nch(i)
        Y =  [Y; Signal{i}(j:n(i):end,:)];
        R =  [R T{i}];
        R2 = [R2; Rig{i}(j:n(i):end)];
%         R3 = [R3; STNelec{i}(j)*ones(length(T{i}),1)];
    end
end

figure
for fb=1:Nb
    j                   = ~isnan(Y(:,fb));
    [rho(fb), p(fb)]    = corr(R(j)',Y(j,fb), 'type', 'Spearman');
    [rho2(fb), p2(fb)]  = corr(R2(j),Y(j,fb), 'type', 'Spearman');
%     [rho3(fb), p3(fb)]  = corr(R3(j),Y(j,fb), 'type', 'Spearman');
    [coef(fb,:), S{fb}] = polyfit(R(j)', Y(j,fb), 1);
    [Y2, dY2]           = polyconf(coef(fb,:), x, S{fb}, 'predopt', 'curve');

    %% Plot results
    pn = [7 1 2 3 8 4 5 6];
    subplot(2,2,fb)
    for i=pn
        for j=1:Nch(i)
%             switch STNelec{i}(j)
%                 case 0
%                     linetyp = '.-';
%                 case 1
%                     linetyp = '--';
%                 case 2
%                     linetyp = '-';
%             end
            plot(T{i}, Signal{i}(j:n(i):end,fb))%, linetyp)
            title(LabFb{fb})
            hold all
        end
    end
    xlim([-2 12])
    xlabel('mins after apomorphine injection')
    ylabel('$\left<P_n-P_0\right>/\left<P_0\right>$', 'Interpreter', 'Latex')
%     plot(x,Y2,'-k', x,Y2-dY2,'--k', x,Y2+dY2,'--k', 'Linewidth', 2)
    plot(x,Y2,'-k', 'Linewidth', 2)
end

fprintf('Temporal correlation:'), [rho;p],
fprintf('Rigor correlation:'), [rho2;p2],

figure,
fb = 4;
j = ~isnan(Y(:,fb));
x=[-10:110];
[coef(fb,:), S{fb}] = polyfit(R2(j), Y(j,fb), 1);
[Y2, dY2] = polyconf(coef(fb,:), x, S{fb}, 'predopt', 'curve');
yr=polyval(coef(fb,:),x);
plot(R2(j),Y(j,fb),'k+')
hold all, plot(x,yr+dY2,'k--', x,yr-dY2,'k--')
rsquare(R2(j),Y(j,fb),1)