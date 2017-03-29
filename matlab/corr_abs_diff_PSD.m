clear all
nch = 4;
cutoff = 1e-10;
updrs_imp = [38 40 55 41 45 47 30 43];
rigor=[100 100 80 30 100 50 20 50];

exp = {'rest', 'fist', 'hold'};

ii = 3;

id{1} = (1:8);
id{2} = [1 2 4 5 7 8];
id{3} = [1 2 3 4 5 7 8];
updrs_imp=updrs_imp(id{ii});
rigor=rigor(id{ii});

task = exp{ii};
load(['LFP_pow_v1_' task '.mat'])
% P_inc_off=P_tot_off-P_coh_off-P_vc_off;
% P_inc_on=P_tot_on-P_coh_on-P_vc_on;
val = cutoff;
P_inc_off(P_inc_off<=cutoff) = val;
P_inc_on(P_inc_on<=cutoff)   = val;
P_coh_off(P_coh_off<=cutoff) = val;
P_coh_on(P_coh_on<=cutoff)   = val;
P_vc_off(P_vc_off<=cutoff)   = val;
P_vc_on(P_vc_on<=cutoff)     = val;
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
% Define frequency bands
% fband{1} = find(f>=1  & f<=4 );
% fband{2} = find(f>=4  & f<=7 );
% fband{3} = find(f>=7  & f<=13);
fband{1} = find(f>=13 & f<=20);
fband{2} = find(f>=20 & f<=30);
fband{3} = find(f>=30 & f<=40);
fband{4} = find(f>=60 & f<=90);% & ~(f>40 & f<60));
Nb = length(fband);
for i = 1:Nb
    j = fband{i};
    dP_inc(:,i) = squeeze( nanmean( yinc(j,:)-xinc(j,:) ) ...
        ./ nanmean( totx(j,:) ) )*100;
    dP_coh(:,i) = squeeze( nanmean( ycoh(j,:)-xcoh(j,:) ) ...
        ./ nanmean( totx(j,:) ) )*100;
    dP_vc(:,i)  = squeeze( nanmean( yvc(j,:) -xvc(j,:) ) ...
        ./ nanmean(totx(j,:) ) )*100;
    dP_bip(:,i) = squeeze( nanmean( ybip(j,:)-xbip(j,:) ) ...
        ./ nanmean(xbip(j,:) ) )*100;
end

X = []; X2 =[];
for i=1:length(rigor); %loop over patients
    dinc(i,:) = nanmean(dP_inc((i-1)*4+1:i*4,:)); 
    dcoh(i,:) = nanmean(dP_coh((i-1)*4+1:i*4,:));
    dvc(i,:) = nanmean(dP_vc((i-1)*4+1:i*4,:));
    dbip(i,:) = nanmean(dP_bip((i-1)*4+1:i*4,:));
    X  = [X rigor([i i i i])]; 
    X2 = [X2 rigor([i i i i i i])];
end
for i=1:Nb %loop over freq.band
    [cr_inc(i)  pr_inc(i)]  = corr(rigor',dinc(:,i), 'type', 'Spearman');
    [cr_coh(i) pr_coh(i)] = corr(rigor',dcoh(:,i), 'type', 'Spearman');
    [cr_vc(i) pr_vc(i)] = corr(rigor',dvc(:,i), 'type', 'Spearman');
    [cr_bip(i) pr_bip(i)] = corr(rigor',dbip(:,i), 'type', 'Spearman');
    [cu_inc(i)  pu_inc(i)] = corr(updrs_imp',dinc(:,i), 'type', 'Spearman');
    [cu_coh(i) pu_coh(i)]  = corr(updrs_imp',dcoh(:,i), 'type', 'Spearman');
    [cu_vc(i) pu_vc(i)] = corr(updrs_imp',dvc(:,i), 'type', 'Spearman'); 
    [cu_bip(i) pu_bip(i)] = corr(updrs_imp',dbip(:,i), 'type', 'Spearman');
    [C(1,i) P(1,i)] = corr(X(~isnan(dP_inc(:,i)))',dP_inc(~isnan(dP_inc(:,i)),i), 'type', 'Spearman');
    [C(2,i) P(2,i)] = corr(X(~isnan(dP_coh(:,i)))',dP_coh(~isnan(dP_coh(:,i)),i), 'type', 'Spearman');
    [C(3,i) P(3,i)] = corr(X(~isnan(dP_vc(:,i)))',dP_vc(~isnan(dP_vc(:,i)),i), 'type', 'Spearman');
    [C(4,i) P(4,i)] = corr(X2(~isnan(dP_bip(:,i)))',dP_bip(~isnan(dP_bip(:,i)),i), 'type', 'Spearman');
end

fprintf('rigor correlation, [corr; p-values]')
[cr_inc;cr_coh;cr_vc;cr_bip]
[pr_inc;pr_coh;pr_vc;pr_bip]
% fprintf('UPDRS correlation, [corr; p-values]')
% [cu_inc;cu_coh;cu_vc;cu_bip]
% [pu_inc;pu_coh;pu_vc;pu_bip]

i=4;        % freq. band
x = rigor; %updrs_imp;  % corr to rigor or updrs_imp
limx = [-5 110]; % [-5 110] for rigor, [25 50] for updrs

fig1 = figure('Papersize', [6.5 5], 'PaperPosition', [0.5 0.5 5.5 4], ...
        'PaperPositionmode', 'manual', 'Visible', 'off');

% plot(updrs_imp',dinc(:,i)','ok', updrs_imp',dcoh(:,i)','+k', ...
%     updrs_imp',dvc(:,i)','xk', updrs_imp',dbip(:,i)','^k')
% xlim([28 47])
% legend(['loc.inc., c=' mat2str(cu_inc(i),2)], ...
%        ['loc.coh., c=' mat2str(cu_coh(i),2)], ...
%        ['vol.con., c=' mat2str(cu_vc(i),2)], ...
%        ['bipolar , c=' mat2str(cu_bip(i),2)], ...
%        'location', 'southwest')
% xlabel('Pre-op UPDRS improvement [%]')
% ylabel('Average \DeltaP/P_{OFF} [%] per patient')
% a = polyfit(x',dvc(:,i),1);

plot(x',dinc(:,i)','ok', x',dcoh(:,i)','+k', ...
    x',dvc(:,i)','xk', x',dbip(:,i)'/4,'^r')
xlim(limx), ylim([-20 10])
legend(['loc.inc., c=' mat2str(cr_inc(i),2)], ...
       ['loc.coh., c=' mat2str(cr_coh(i),2)], ...
       ['vol.con., c=' mat2str(cr_vc(i),2)], ...
       ['bipolar , c=' mat2str(cr_bip(i),2)], ...
       'location', 'southwest')
xlabel('Rigidity improvement [%]')
ylabel('Average \DeltaP/P_{OFF} [%] per patient')
a = polyfit(x',dinc(:,i),1);


hold all
plot(limx, a(1)*limx+a(2), '--k')
a2 = axes('YAxisLocation', 'Right');
set(a2, 'color', 'none', 'ycolor', 'r', 'XTick', [], 'YLim', [-80 30])
ylabel('Average \DeltaP/P [%] for bipolar')
grid on

% % figure
% % plot(updrs_imp',dacoh(:,4)','ok', updrs_imp',dacoh(:,5)','+k', ...
% %     updrs_imp',dacoh(:,6)','xk'), xlim([28 47])
% % legend(['\beta,  c=' mat2str(cu_acoh(4),2)], ...
% %        ['\beta_1, c=' mat2str(cu_acoh(5),2)], ...
% %        ['\beta_2, c=' mat2str(cu_acoh(6),2)])
% % title('Coherent signals, fist: UPDRS improvement-\beta correlation')
% % xlabel('Pre-op UPDRS improvement OFF->ON [%]')
% % ylabel('Change in \beta band OFF->ON [%]')
% saveas(fig1, 'corr_beta1_fist.fig')
% saveas(fig1, ['corr_beta_' task '_ds10_pc15.eps'], 'epsc')
% close(fig1)