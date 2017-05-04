%% Parameter
nsig  = 6;
w0    = 12;

alpha = 10;  %coh
beta  = 20;  %vc
gamma = 50;  %coh

phase_thresh = 15;
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

% Determine TrueMat, the matrix where correct results are stored in order
% to calculate false pos/negs of coherent estimate
TrueMat = zeros(length(f),nt);
switch w0
    case 12
        w0fac=1.14;
end
TrueMat(f>alpha/w0fac & f<alpha*w0fac,nt/5+1:4*nt/5) = 1; %coh
for i=nt/2+1:nt
    TrueMat(f>beta/w0fac+i/nt*10 & f<beta*w0fac+i/nt*10,i) = 1; %coh
end
TrueMat(f>gamma/w0fac & f<gamma*w0fac,1:nt/2) = 1; %coh
true_pos_ind = find(TrueMat==1);
true_pos = sum(TrueMat(:)==1);
true_neg = sum(TrueMat(:)==0);


sig   = sig_coh_thresh(w0, nsig);

%% Preallocate
Ptot = zeros(length(f),N);
Pinc = zeros(length(f),N);
Pcoh = zeros(length(f),N);
Pvc  = zeros(length(f),N);
FPR  = zeros(N,1);
TPR  = zeros(N,1);

fprintf ('\r Processing ')
for j=1:N
    fprintf('%4d/%4d', j, N)
    
    %% Synthetic time series
    [x, y] = synth_ts( t, nl, A, dp );
    
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
    [a, b, c]        = pcc ( f, W, C, coi, sig, Wxy, 0, phase_thresh );
    Pcoh(:,j)        = a(:,1,2);
    Pinc(:,j)        = b(:,1,2);
    Pvc(:,j)         = c(:,1,2);
    
    %% False/true positive rate
    Ph       = abs(angle(Wxy)/pi*180);
    Cthreshd = C(:,:,1,2)>sig & Ph(:,:,1,2)>=phase_thresh ...
                & Ph(:,:,1,2)<=180-phase_thresh;
    FalseMat = TrueMat - Cthreshd;
    FP       = sum(FalseMat(:)==-1);
    TPR(j)   = sum(FalseMat(true_pos_ind)==0)/true_pos;
    FPR(j)   = FP/true_neg;
    
    fprintf('\b\b\b\b\b\b\b\b\b')

end


save('PSD_synth_PCC_ns6_w12_nl3_v2.mat', 'A', 'N', 'phase_thresh', ...
    'Pcoh', 'Pinc', 'Pvc', 'Ptot', 'alpha', 'beta', 'gamma', ...
    'coi', 'f', 'FPR', 'TPR', 'w0', 'nsig', 'nl', 'dp')


% % One more time for data used in figure
% x            = (xa + xb + xg + n1) / A;
% y            = (ya + yb + yg + n2) / A;
% scale        = (w0+sqrt(2+w0^2))/4/pi ./ f;
% [ts, W, coi] = procdata([x y], 'freq', f, 'w0', w0, 'dt', dt, 'filt', [0 0.5]);
% [C, Wxy, W]  = wave_cohere(W, scale, nsig, ds, dt);
%     
% 
% %% Plot results
% figure
% i=1;
% subplot(1,2,1)
% Ph  = abs(angle(Wxy(:,:,1,2))/pi*180);
% coh = C(:,:,1,2)>sig & Ph>phthres & Ph<180-phthres;
% vc  = C(:,:,1,2)>sig & (Ph<=phthres | Ph>=180-phthres);
% Z   = zeros(size(W(:,:,1)));
% Z(coh) = 1;
% Z(vc ) = 2;
% pcolor(t(1:ds2:end),f,Z(:,1:ds2:end)), shading flat
% hold all, plot(t(1:ds2:end),nsig./coi(1:ds2:end,1),'-k')
% caxis([0 2])
% colormap([1 0 0; 0 0 1; 0 1 0])
% colorbar('YTickLabel',{'Inc'; 'Coh'; 'VC'}, 'YTick', [1/3 1 5/3])
% title('PCC in time-frequency domain')
% ylim([1 max(f)]), xlim([0 max(t)])
% xlabel('Time [s]')
% ylabel('f [Hz]')
% subplot(1,2,2)
% semilogy(f, Pinc(:,i), '-r', f, Pcoh(:,i), '-b', f, Pvc(:,i), '-g', f, Ptot(:,i), '--k')
% title('PSD of synthetic data')
% xlim([1 max(f)]), ylim(limy)
% ylabel('PSD [au]')
% legend('Incoherent', 'Coherent', 'Volume-conduction', 'Total power')
% xlabel('f [Hz]')