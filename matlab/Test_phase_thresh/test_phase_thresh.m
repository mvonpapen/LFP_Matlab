%% Parameter
nsig  = 6;
w0    = 12;

alpha = 10;  %coh
beta  = 20;  %vc
gamma = 50;  %coh

phthres = 15;
nl    = 3;
dp    = 30/180*pi;
f     = 1:70;
ds    = 1;
ds2   = 150; %downsampling for plot

dt    = 1/2500;
nt    = 20/dt;
t     = (0:nt-1)*dt;
t     = t(1:ds:end);

A     = 1000; % amplitude correction
limy  = [1e-8 1e-6];
N     = 1000; j=1;

sig   = sig_coh_thresh(w0, nsig);

%% Preallocate
Ptot = zeros(length(f),N);
Pinc = zeros(length(f),N);
Pcoh = zeros(length(f),N);
Pvc  = zeros(length(f),N);

parfor j=1:N
    
    %% Synthetic time series
    [x, y] = synth_ts_PM_paper( t, nl, A, dp );
    
%     %% Real LFP
%     load('data_rest_v3.mat', 'LFP_OFF', 'ArtLFP_OFF');
%     x  = LFP_OFF{2}(:,2);
%     y  = LFP_OFF{2}(:,4);
%     nt = length(x);
%     t     = (0:nt-1)*dt;
%     limy  = [1e-9 1e-5];


    %% Spectral analysis
    scale            = (w0+sqrt(2+w0^2))/4/pi ./ f;
    [~, W, coi, tmp] = procdata([x y], 'freq', f, 'w0', w0, 'dt', dt);
    Ptot(:,j)        = tmp(:,1);
    [C, Wxy, W]      = wave_cohere(W, scale, nsig, ds, dt);
    coi              = coi(1:ds:end,:);
    [a, b, c]        = pcc ( f, W, C, coi, sig, Wxy, 0, phthres );
    Pcoh(:,j)        = a(:,1,2);
    Pinc(:,j)        = b(:,1,2);
    Pvc(:,j)         = c(:,1,2);

end
clear tmp tmp1 tmp2 tmp3 a b c
clear xa xb xg ya yb yg n1 n2


% One more time for data used in figure
x            = (xa + xb + xg + n1) / A;
y            = (ya + yb + yg + n2) / A;
scale        = (w0+sqrt(2+w0^2))/4/pi ./ f;
[ts, W, coi] = procdata([x y], 'freq', f, 'w0', w0, 'dt', dt, 'filt', [0 0.5]);
[C, Wxy, W]  = wave_cohere(W, scale, nsig, ds, dt);
    

%% Plot results
figure
i=1;
subplot(1,2,1)
Ph  = abs(angle(Wxy(:,:,1,2))/pi*180);
coh = C(:,:,1,2)>sig & Ph>phthres & Ph<180-phthres;
vc  = C(:,:,1,2)>sig & (Ph<=phthres | Ph>=180-phthres);
Z   = zeros(size(W(:,:,1)));
Z(coh) = 1;
Z(vc ) = 2;
pcolor(t(1:ds2:end),f,Z(:,1:ds2:end)), shading flat
hold all, plot(t(1:ds2:end),nsig./coi(1:ds2:end,1),'-k')
caxis([0 2])
colormap([1 0 0; 0 0 1; 0 1 0])
colorbar('YTickLabel',{'Inc'; 'Coh'; 'VC'}, 'YTick', [1/3 1 5/3])
title('PCC in time-frequency domain')
ylim([1 max(f)]), xlim([0 max(t)])
xlabel('Time [s]')
ylabel('f [Hz]')
subplot(1,2,2)
semilogy(f, Pinc(:,i), '-r', f, Pcoh(:,i), '-b', f, Pvc(:,i), '-g', f, Ptot(:,i), '--k')
title('PSD of synthetic data')
xlim([1 max(f)]), ylim(limy)
ylabel('PSD [au]')
legend('Incoherent', 'Coherent', 'Volume-conduction', 'Total power')
xlabel('f [Hz]')