%% Compute wavelet magnitude squared coherencE (!)
%%
%%  INPUT:
%%      W       Wavelet coefficients nf x nt x nch
%%      fw      Frequency vector in Hz

function [C, sWxy, sW, Wxy] = wave_acoh3 ( W, scale, t_smo, dt )

%% Parameters
if nargin<4; dt=1/2500; end
if nargin<3; t_smo=6; end
[nf, nt, nch] = size ( W );

%% Setting parameters
scale = scale(:);
% tmp = 1./scale;
fac  = 1; %tmp(:,ones(1,nt)); %factor 1/s to multiply with W (this is needed if smoothing is done over scales)


%% Set variables
sW   = zeros(nf,nt,nch);
Wxy  = zeros(nf,nt,nch,nch);
sWxy = zeros(nf,nt,nch,nch);
C    = zeros(nf,nt,nch,nch);

%% Cross-spectra, Smoothing and Coherence (!)
for i=1:nch;
    sW(:,:,i)=smooth_mat(fac.*abs(W(:,:,i)).^2,dt,t_smo*scale);
end

for i=1:nch;
    for j=i+1:nch;
        Wxy(:,:,i,j)  = squeeze(W(:,:,i).*conj(W(:,:,j)));
        sWxy(:,:,i,j) = smooth_mat(fac.*Wxy(:,:,i,j),dt,t_smo*scale);
        C(:,:,i,j)    = squeeze(abs(sWxy(:,:,i,j)).^2)./...
            squeeze(sW(:,:,i).*sW(:,:,j));
        C(:,:,j,i)    = C(:,:,i,j);
        Wxy(:,:,j,i)  = conj(Wxy(:,:,i,j));
        sWxy(:,:,j,i) = conj(sWxy(:,:,i,j));
    end
end