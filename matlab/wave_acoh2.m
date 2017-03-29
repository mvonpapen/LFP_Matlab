%% Compute wavelet magnitude squared coherencE (!)
%%
%%  INPUT:
%%      W       Wavelet coefficients nf x nt x nch
%%      fw      Frequency vector in Hz

function [C, sWxy, sW, Wxy] = wave_acoh2 ( W, scale, t_smo, sc_smo )

%% Parameters
if nargin<4; sc_smo=0; end
if nargin<3; t_smo=3; end
[nf, nt, nch] = size ( W );

%% Setting parameters
dt  = 1/2500;
scale = scale(:);
tmp = 1./scale;
fac  = tmp(:,ones(1,nt)); %factor 1/s to multiply with W (this is needed if smoothing is done over scales)
clear tmp


%% Set variables
sW   = NaN(nf,nt,nch);
Wxy  = NaN(nf,nt,nch,nch);
sWxy = NaN(nf,nt,nch,nch);
C    = NaN(nf,nt,nch,nch);

%% Cross-spectra, Smoothing and Coherence (!)
for i=1:nch;
    sW(:,:,i)=smooth_mat(fac.*abs(W(:,:,i)).^2,dt,t_smo*scale,...
        'ScaleSmo', sc_smo);
end

for i=1:nch;
    for j=i+1:nch;
        Wxy(:,:,i,j)  = squeeze(W(:,:,i).*conj(W(:,:,j)));
        sWxy(:,:,i,j) = smooth_mat(fac.*Wxy(:,:,i,j),dt,t_smo*scale,...
            'ScaleSmo', sc_smo);
        C(:,:,i,j)    = squeeze(abs(sWxy(:,:,i,j)).^2)./...
            squeeze(sW(:,:,i).*sW(:,:,j));
        C(:,:,j,i)    = C(:,:,i,j);
        Wxy(:,:,j,i)  = conj(Wxy(:,:,i,j));
        sWxy(:,:,j,i) = conj(sWxy(:,:,i,j));
    end
end

clear dj fac n nf kl ke i j s0 j1

% %% Coherence can be plotted with
% plot_coh_phase(ti, 1./p, C(:,:,1,2), 'coi', coi)
