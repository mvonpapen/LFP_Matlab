%% Plots Figure 8 of JoN Methods paper
%% panel a: LFP A-P of BIMA, b: PCC in time-freq., c: PSD of signals


load data_rest_v3

%% Parameters
sf      = 2456;
f       = logspace(log10(1/50),log10(sf/2),200);
Nf      = length(f);
w0      = 12;
nsig    = 6;
usenan  = false;
sig     = sig_coh_thresh(w0, nsig);
phthres = sig_phi_thresh(w0, nsig);
ds      = 20;
ds2     = 20;
scale   = (w0+sqrt(2+w0^2))/4/pi ./ f;
patnum  = 2;
i=2;j=3;



[X,W1,coi1,P_tot_off(:,1:Nch(patnum),patnum)] = procdata(LFP_OFF{patnum}, 'freq', f, 'w0', w0, ...
    'filter', [44 56; 90 110], 'art', ArtLFP_OFF{patnum}); %'filter', [44 56; 90 110],
[Cs, Wxys, W1s] = wave_cohere ( W1, scale, nsig, ds );
coi1s = coi1(1:ds:end,:)/nsig;
[Pcoh, Pinc, Pvc, Pnvc] = psd_acoh ( f, W1s, Cs, coi1s, sig, Wxys, usenan, phthres );
Ph = angle(Wxys)/pi*180;

Pinc(Pinc<1e-10) = NaN;
Pcoh(Pcoh<1e-10) = NaN;
Pvc(Pvc<1e-10) = NaN;
t=(0:length(Wxys)-1)/sf*ds;


% Define W_coh, W_inc and W_vc
coh = Cs(:,:,i,j)>sig & abs(Ph(:,:,i,j))>phthres;
inc = Cs(:,:,i,j)<=sig;
vc = Cs(:,:,i,j)>sig & abs(Ph(:,:,i,j))<=phthres;
Z   = zeros(size(Wxys(:,:,i)));
Z(coh) = 1;
Z(vc ) = 2;
Wcoh = W1(:,1:ds:end,i);
Wcoh(~coh) = 0;
Winc = W1(:,1:ds:end,i);
Winc(~inc) = 0;
Wvc = W1(:,1:ds:end,i);
Wvc(~vc) = 0;


% Reconstruct time series
[ts_coh, P_coh] = wave_recstr ( Wcoh, f, w0 );
[ts_inc, P_inc] = wave_recstr ( Winc, f, w0 );
[ts_vc,  P_vc]  = wave_recstr ( Wvc, f, w0 );
[ts_tot, P_tot, varWT] = wave_recstr ( W1(:,1:ds:end,i), f, w0 );

figure
[~,fi]=min(abs(f-4.7));
subplot(3,2,[1 2])
X  = X(:,i);
t0 = (1:length(X)) / sf;
sigma2 = var(X);
subplot(4,2,[1 2])
[ax,h1,h2] = plotyy(t0,X,t,W1s(fi,:,i));
set(h1,'color', [0.8 0.8 0.8]); 
set(h2,'color', 'k', 'linew', 2)
ax(1).YAxis.Color = 'k';
ax(2).YAxis.Color = 'k';
ax(2).YAxis.Label.String = 'PSD ($f=5\,$Hz)';
ax(1).YAxis.Label.String = 'LFP A [V]';
ax(1).YAxis.Label.Interpreter = 'latex';
ax(2).YAxis.Label.Interpreter = 'latex';
ax(1).XAxis.Limits = [0 48];
ax(1).YAxis.Limits = [-0.019 0.019];
ax(2).XAxis.Limits = [0 48];
set(ax(2),'Ysca','log', 'Ylim', [1e-7 1e-4], 'YTick', [1e-7 1e-6 1e-5 1e-4]);
title('LFP A during tremor episode of patient 3', 'interp', 'latex')

fprintf('variance of ts/WT = %f/%f. Error of %f%%.\n', sigma2, varWT, (sigma2-varWT)/sigma2)
subplot(4,2,[3 4])
plot(t, ts_vc+0.015, 'g-', t, ts_inc-0.015, 'r-', t, ts_coh, 'b-')
xlim([0 48])
ylim([-0.025 0.025]);
xlabel('time [s]', 'interp', 'latex')
ylabel('PCC signals [V]', 'interp', 'latex')

% figure
subplot(4,2,[5 7])
% subplot(1,2,1)
pcolor(t(1:ds2:end),f,Z(:,1:ds2:end)), shading flat
hold all
plot(t(1:ds2:end),1./coi1s(1:ds2:end,i),'-k')
caxis([0 2])
colormap([1 0 0; 0 0 1; 0 1 0])
colorbar('YTickLabel',{'Inc'; 'Coh'; 'VC'}, 'YTick', [1/3 1 5/3])
title('PCC in time-frequency domain', 'interp', 'latex')
set(gca, 'TickDir', 'out', 'Ysca', 'log');
ylim([1 100])
xlabel('time [s]', 'interp', 'latex')
ylabel('$f$ [Hz]', 'interp', 'latex')

subplot(4,2,[6 8])
semilogy(f,P_tot_off(:,i,patnum),'-k', f,Pvc(:,i,j),'-g', f,Pcoh(:,i,j),'b-', f,Pinc(:,i,j),'-r')
xlim([1 30])
xlabel('$f$ [Hz]', 'interp', 'latex')
ylabel('PSD [V$^2$/Hz]', 'interp', 'latex')
title('PSD', 'interp', 'latex')
h = legend('total', 'vc', 'coh', 'inc');
set(h, 'interp', 'latex')

saveas(gcf, 'LFP_BIMA_Ap_subsig.fig')