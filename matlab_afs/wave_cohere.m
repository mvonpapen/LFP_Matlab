function [C, sWxy, sW2, Wxy, Cohy] = wave_cohere ( W, scale, nsig, ds, dt )
%% WAVE_COHERE Compute wavelet magnitude squared coherence (!) from a single wavelet coefficient matrix
% 
%   [C, sWxy, sW, Wxy] = WAVE_COHERE( W, scale, nsig, ds, dt )
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
%           ImCoh:   Imaginary part of coherency (nf x nt x nchx x nchy)
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
sW2    = zeros(nf,nt,nch);
Wxy   = zeros(nf,nt,nch,nch);
sWxy  = zeros(nf,nt,nch,nch);
if ds==1
    Cohy    = zeros(nf,nt,nch,nch);
else
    nt2  = ceil(nt/ds);
    Cohy    = zeros(nf,nt2,nch,nch);
end
    

%% Cross-spectra, Smoothing and Coherence
for i=1:nch
    sW2(:,:,i)=smooth_mat(abs(W(:,:,i)).^2,dt,nsig*scale);
end
for i=1:nch
    for j=i+1:nch
        Wxy(:,:,i,j)  = squeeze(W(:,:,i).*conj(W(:,:,j)));
        sWxy(:,:,i,j) = smooth_mat(Wxy(:,:,i,j),dt,nsig*scale);
        if nargout>3
            Wxy(:,:,j,i)  = conj(Wxy(:,:,i,j));
        end
        sWxy(:,:,j,i) = conj(sWxy(:,:,i,j));
        if ds==1
            Cohy(:,:,i,j)    = sWxy(:,:,i,j) ./ sqrt(sW2(:,:,i).*sW2(:,:,j));
            Cohy(:,:,j,i)    = Cohy(:,:,i,j);
        end
    end
end

% Downsample data, take every dec'th value
if ds~=1
    sW2  = sW2(:,1:ds:end,:); % sW2 already in units of power [V^2/Hz]
    sWxy = sWxy(:,1:ds:end,:,:);
    for i=1:nch
        for j=i+1:nch
            Cohy(:,:,i,j) = sWxy(:,:,i,j) ./ sqrt(sW2(:,:,i).*sW2(:,:,j));
            Cohy(:,:,j,i) = Cohy(:,:,i,j);
        end
    end
end
sW2 = 2*dt*sW2; %make sW in units of power

C = abs(Cohy).^2;
% %% Coherence can be plotted with
% plot_coh_phase(ti, 1./p, C(:,:,1,2), 'coi', coi)