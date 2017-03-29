function [C, sWxy, sW, Wxy] = wave_cohere ( W, scale, nsig, ds, dt )
%% WAVE_COHERE Compute wavelet magnitude squared coherence (!) from a single wavelet coefficient matrix
% 
%   [C, sWxy, sW, Wxy] = WAVE_COH( W, scale, nsig, ds, dt )
% 
%	INPUT:
%           W:       Wavelet coefficients of channel x (nf x nt x nchx)
%           scale:   Scale vector in s
%           nsig:    determines size of Gaussian for smoothing [nsig=6]
%           dt:      sample interval [dt=1/2456]
%
%	OUTPUT:
%           C:       Magnitude squared Coherence (nf x nt x nchx x nchy)
%           sWxy:    Smoothed cross-spectra of Wx and Wy (nf x nt x nchx x nchy)
%           sW:      Smoothed auto-spectrum of W (nf x nt x nchi)
%           Wxy:     Cross-spectra of W (nf x nt x nchx x nchy), not smoothed
%
% Author: Michael von Papen
% Date: 22.04.16

%% Parameters
if nargin<5; dt   = 1/2456; end
if nargin<4; ds   = 1; end
if nargin<3; nsig = 3; end
[nf, nt, nch] = size ( W );

%% Setting variables
scale = scale(:);
sW    = zeros(nf,nt,nch);
Wxy   = zeros(nf,nt,nch,nch);
sWxy  = zeros(nf,nt,nch,nch);
if ds==1
    C    = zeros(nf,nt,nch,nch);
else
    nt2  = ceil(nt/ds);
    C    = zeros(nf,nt2,nch,nch);
end
    

%% Cross-spectra, Smoothing and Coherence
for i=1:nch;
    sW(:,:,i)=smooth_mat(abs(W(:,:,i)).^2,dt,nsig*scale);
end
for i=1:nch
    for j=i+1:nch;
        Wxy(:,:,i,j)  = squeeze(W(:,:,i).*conj(W(:,:,j)));
        sWxy(:,:,i,j) = smooth_mat(Wxy(:,:,i,j),dt,nsig*scale);
        if nargout>3
            Wxy(:,:,j,i)  = conj(Wxy(:,:,i,j));
        end
        sWxy(:,:,j,i) = conj(sWxy(:,:,i,j));
        if ds==1
            C(:,:,i,j)    = abs(sWxy(:,:,i,j)).^2 ./ (sW(:,:,i).*sW(:,:,j));
            C(:,:,j,i)    = C(:,:,i,j);
        end
    end
end

% Downsample data, take every dec'th value
if ds~=1
    sW   = sW(:,1:ds:end,:); % sW already in units of power [V^2/Hz]
    sWxy = sWxy(:,1:ds:end,:,:);
    for i=1:nch
        for j=i+1:nch
            C(:,:,i,j) = abs(sWxy(:,:,i,j)).^2 ./ (sW(:,:,i).*sW(:,:,j));
            C(:,:,j,i) = C(:,:,i,j);
        end
    end
end
sW = 2*dt*sW; %make sW in units of power

% %% Coherence can be plotted with
% plot_coh_phase(ti, 1./p, C(:,:,1,2), 'coi', coi)