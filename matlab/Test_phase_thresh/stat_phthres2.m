function [Prcoh, Prvc, dp] = stat_phthres2(f, phthres, w0, nsig, N)
%STAT_PHTHRES2 statistically estimates coherent and vol.cond. power
%   STAT_PHTHRES2 can be used to estimate resolution of specific phase threshold

%   Date: 20.01.2016
%   Author: M. von Papen

if nargin<5
    N = 500;
end
if nargin<4
    nsig = 3;
end
if nargin<3
    w0 = 12;
end
if nargin<2
    phthres = 15;
end


%% Parameter
nl    = 3;                  %noise level 1 -> powlawnoise ~ 0.1*[1:1e3].^(-1)
dt    = 1/2500;
t     = (1:20/dt)*dt;
nt    = length(t);
dp    = 0:40;               %delta phase
sig   = sig_coh_thresh(w0, nsig);
ds    = 10;                 %downsampling


%% Preallocate variables
Prcoh = zeros(length(dp),N);
Prvc  = zeros(length(dp),N);

pool = ~isempty(gcp('nocreate'));
if ~pool
    mypool = parpool(4);
end

%% Begin with loops
for j = 1:length(dp)
    
    delphase = dp(j)/180*pi;
    scale    = (w0+sqrt(2+w0^2))/4/pi ./ f;
    fprintf('Phase %u/%u.\n', dp(j), max(dp));
    
    parfor n = 1:N

        %% Generate noise
        N1 = powlawnoise(nt,1)';
        N2 = powlawnoise(nt,1)';

        %% Time series [ 1/sqrt(f) leads to SNR=10: P(sin)=10*P(powlawnoise) ]
        x = sin(2*pi*f*t') + nl * N1;
        y = sin(2*pi*f*t'+delphase) + nl * N2;


        %% Spectral analysis
        [~, W, coi]    = procdata([x y],'freq', f, 'w0', w0);
        [C, Wxy, W]    = wave_cohere(W, scale, nsig, ds);
        coi            = coi(1:ds:end,:)/nsig;
        [coh, inc, vc] = psd_acoh ( f, W, C, coi/nsig, sig, Wxy, 0, phthres );
        Prcoh(j,n)     = coh(1,1,2) / (coh(1,1,2)+inc(1,1,2)+vc(1,1,2));
        Prvc(j,n)      = vc(1,1,2) / (coh(1,1,2)+inc(1,1,2)+vc(1,1,2));

    end
    
end

if ~pool
    delete(mypool);
end