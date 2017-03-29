nch = 4;
cutoff = 1e-10;
updrs_imp = [38 40 55 41 45 47 30 43];
rigor=[100 100 80 30 100 50 20 50];
x = rigor;

exp = {'rest', 'hold', 'fist'};
ii = 3;
fist= [1 2 4 5 7 8];
Hold = [1 2 3 4 5 7 8];
x = x(fist);
task = exp{ii};

load(['LFP_pow_v1_' task '.mat'])

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
fband{1} = find(f>=13 & f<=30);
fband{2} = find(f>=13 & f<=20);
fband{3} = find(f>=20 & f<=30);
fband{4} = find(f>=30 & f<=90 & ~(f>40 & f<60));
Nb = length(fband);

for i = 1:Nb
    j = fband{i};
    dP_inc(:,i) = squeeze( trapz( f(j), yinc(j,:)-xinc(j,:) ) ...
    ./ trapz( f(j),totx(j,:) ) )*100;
    dP_coh(:,i) = squeeze( trapz( f(j), ycoh(j,:)-xcoh(j,:) ) ...
    ./ trapz( f(j),totx(j,:) ) )*100;
    dP_vc(:,i)  = squeeze( trapz( f(j), yvc(j,:) -xvc(j,:) ) ...
    ./ trapz( f(j),totx(j,:) ) )*100;
    dP_bip(:,i) = squeeze( trapz( f(j), ybip(j,:)-xbip(j,:) ) ...
    ./ trapz( f(j),xbip(j,:) ) )*100;
end

for i=1:length(x); %loop over patients
    dinc(i,:) = nanmean(dP_inc((i-1)*4+1:i*4,:)); 
    dcoh(i,:) = nanmean(dP_coh((i-1)*4+1:i*4,:));
    dvc(i,:)  = nanmean(dP_vc((i-1)*4+1:i*4,:));
    dbip(i,:) = nanmean(dP_bip((i-1)*4+1:i*4,:));
end
for i=1:Nb %loop over freq.band
    [c_inc(i) p_inc(i)] = corr(x',dinc(:,i), 'type', 'Spearman');
    [c_coh(i) p_coh(i)] = corr(x',dcoh(:,i), 'type', 'Spearman');
    [c_vc(i)  p_vc(i)]  = corr(x',dvc(:,i), 'type', 'Spearman');
    [c_bip(i) p_bip(i)] = corr(x',dbip(:,i), 'type', 'Spearman');
end

fprintf('rigor correlation, [corr; p-values]')
[c_inc;c_coh;c_vc;c_bip]
[p_inc;p_coh;p_vc;p_bip]

i=2;        % freq. band
limx = [-5 110]; % [-5 110] for rigor, [25 50] for updrs

fig1 = figure('Papersize', [6.5 5], 'PaperPosition', [0.5 0.5 5.5 4], ...
        'PaperPositionmode', 'manual', 'Visible', 'off');

plot(x',dinc(:,i)','ok', x',dcoh(:,i)','+k', ...
    x',dvc(:,i)','xk', x',dbip(:,i)'/4,'^r')
xlim(limx), ylim([-20 5])
legend(['loc.inc., c=' mat2str(c_inc(i),2)], ...
       ['loc.coh., c=' mat2str(c_coh(i),2)], ...
       ['vol.con., c=' mat2str(c_vc(i),2)], ...
       ['bipolar , c=' mat2str(c_bip(i),2)], ...
       'location', 'southwest')
xlabel('Rigidity improvement [%]')
ylabel('Average \DeltaP/P_{OFF} [%] per patient')
a = polyfit(x',dcoh(:,i),1);


hold all
plot(limx, a(1)*limx+a(2), '--k')
a2 = axes('YAxisLocation', 'Right');
set(a2, 'color', 'none', 'ycolor', 'r', 'XTick', [], 'YLim', [-80 20])
ylabel('Average \DeltaP/P [%] for bipolar')
grid on