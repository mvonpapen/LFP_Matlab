function [Prcoh, Prvc, Princ, dp] = stat_phthres(f, phthres, N)
%STAT_PHTHRES statistically estimates optimal phase threshold
%   STAT_PHTHRES plots the resulting statistics

%   Date: 18.01.2016
%   Author: M. von Papen


%% Parameter
% phthres = 1:20;
% f     = [1 5 20];
nl    = 3;                  %noise level 1 -> powlawnoise ~ 0.1*[1:1e3].^(-1)
dt    = 1/2500;
nt    = 40/dt;
t     = (0:nt-1)*dt;
dp    = [0:40];               %delta phase
sig   = 0.43;               %p=0.01 significance level for w0=12, tsmo=6
tsmo  = 6;
w0    = 12;
dsf   = 1;                  %downsampling factor

for n = N:-1:1

    N1 = powlawnoise(nt,1)';
    N2 = powlawnoise(nt,1)';
%     for i = length(f):-1:1
        parfor j = 1:length(dp)

            %% Time series [ 1/sqrt(f(i)) leads to SNR=10: P(sin)=10*P(powlawnoise) ]
            x = sin(2*pi*f(i)*t') + nl * N1;
            y = sin(2*pi*f(i)*t'-dp(j)/180*pi) + nl * N2;


            %% Spectral analysis
            scale   = (w0+sqrt(2+w0^2))/4/pi ./ f(i);
            [~, W, coi] = procdata([x y],'freq',f(i), 'w0', w0);
            [C, Wxy, W] = wave_cohere(W,scale,tsmo,dsf);
            coi = coi(dsf:dsf:end,:);
            np   = length(phthres);
            for k = np:-1:1
                [coh, inc, vc] = psd_acoh ( f(i), W, C, coi, sig, Wxy, 0, phthres(k) );
                Prcoh(i,j,k,n) = coh(1,1,2)/(coh(1,1,2)+inc(1,1,2)+vc(1,1,2));
                Prvc(i,j,k,n)  = vc(1,1,2)/(coh(1,1,2)+inc(1,1,2)+vc(1,1,2));
                Princ(i,j,k,n) = inc(1,1,2)/(coh(1,1,2)+inc(1,1,2)+vc(1,1,2));
            end

        end
    end
    
end