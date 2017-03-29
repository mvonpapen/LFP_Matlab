function [C, sWxy, sWx, sWy, Wxy] = wave_coh ( Wx, Wy, scale, nsig, dt )
%% WAVE_COH Compute wavelet magnitude squared coherence (!)
% 
%   [C, sWxy, sWx, sWy, Wxy] = WAVE_COH( Wx, Wy, fw, tSmo, scSmo )
% 
%	INPUT:
%           Wx:      Wavelet coefficients of channel x (nf x nt x nchx)
%           Wy:      Wavelet coefficients of channel y (nf x nt x nchy)
%           fw:      Frequency vector in Hz
%           tSmo:    determines size of Gaussian for smoothing [tSmo=6]
%           scSmo:   determines smoothing across how many scales [scSmo=0]
%
%	OUTPUT:
%           C:       Magnitude squared Coherence (nf x nt x nchx x nchy)
%           sWxy:    Smoothed cross-spectra of Wx and Wy (nf x nt x nchx x nchy)
%           sWi:     Smoothed auto-spectra of Wi (nf x nt x nchi)
%           Wxy:     Cross-spectra of Wx and Wy (nf x nt x nchx x nchy), not
%                    smoothed
%
% Author: Michael von Papen
% Date: 13.10.15

%% Parameters
[m,n,nx] = size(Wx);
[m,n,ny] = size(Wy);

%% Setting parameters
fac=1;



%% Set variables
sWx = NaN(m,n,nx); % ??????WARUM NAN???? ZEROS SCHNELLER
sWy = NaN(m,n,ny);
Wxy = NaN(m,n,nx,ny);
sWxy = NaN(m,n,nx,ny);
C = NaN(m,n,nx,ny);

%% Cross-spectra, Smoothing and Coherence (!)
for i = 1:ny;
    sWy(:,:,i) = smooth_mat(fac.*abs(Wy(:,:,i)).^2,dt,nsig*scale);
end

for i = 1:nx;
    sWx(:,:,i) = smooth_mat(fac.*squeeze(abs(Wx(:,:,i))).^2,dt,nsig*scale); 
    for j = 1:ny;
        Wxy(:,:,i,j) = squeeze(Wx(:,:,i).*conj(Wy(:,:,j))); 
        sWxy(:,:,i,j) = smooth_mat(fac.*Wxy(:,:,i,j),dt,nsig*scale);
        C(:,:,i,j) = squeeze(abs(sWxy(:,:,i,j)).^2)./...
                     squeeze(sWx(:,:,i).*sWy(:,:,j)); 
    end
end

% clear dj fac n m kl ke i j s0 j1 ?????WARUM????
end

% %% Coherence can be plotted with
% i = 1;
% j = 1;
% plot_coh(ti,p,C(:,:,i,j),coi)