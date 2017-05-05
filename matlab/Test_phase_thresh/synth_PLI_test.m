%% Parameter
nsig  = 6;
w0    = 12;
sig   = 0.55; %p01=0.70, p05=0.55 (see 'sig_PLI_ns6_w12.mat')
N     = 1000;

% Parameters for synth data
alpha = 10;  %coh
beta  = 20;  %vc
gamma = 50;  %coh
f     = 1:70;
dt    = 1/2500;
nt    = 20/dt;
t     = (0:nt-1)*dt;
nl    = 3;
dp    = 30/180*pi;
A     = 1000; % amplitude correction


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


%% Preallocate
Ptot = zeros(length(f),N);
Pinc = zeros(length(f),N);
Pcoh = zeros(length(f),N);
FPR  = zeros(N,1);
TPR  = zeros(N,1);

fprintf ('\r Processing ')
for j=1:N
    fprintf('%4d/%4d', j, N)
    
    %% Synthetic time series
    [x, y] = synth_ts( t, nl, A, dp );

    %% Spectral analysis
    scale            = (w0+sqrt(2+w0^2))/4/pi ./ f;
    [~, W, coi, tmp] = procdata([x y], 'freq', f, 'w0', w0, 'dt', dt);
    Ptot(:,j)        = tmp(:,1);
    Wxy              = squeeze(W(:,:,1).*conj(W(:,:,2)));
    PLI              = phase_lag_index(Wxy, scale, nsig, dt);
    [a, b]           = pcc ( f, W, repmat(PLI,1,1,2,2), coi, sig );
    Pcoh(:,j)        = a(:,1,2);
    Pinc(:,j)        = b(:,1,2);
    
    %% False/true positive rate
    PLIthreshd = PLI>sig;
    FalseMat = TrueMat - PLIthreshd;
    FP       = sum(FalseMat(:)==-1);
    TPR(j)   = sum(FalseMat(true_pos_ind)==0)/true_pos;
    FPR(j)   = FP/true_neg;

    fprintf('\b\b\b\b\b\b\b\b\b')
end


save('PSD_synth_PLI_p05_ns6_w12_nl3.mat', 'A', 'N', 'sig', ...
    'Pcoh', 'Pinc', 'Ptot', 'alpha', 'beta', 'gamma', 'nt', ...
    'coi', 'f', 'FPR', 'TPR', 'w0', 'nsig', 'nl', 'dp')