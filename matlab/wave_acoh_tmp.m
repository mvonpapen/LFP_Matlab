%% Compute wavelet magnitude squared coherencE (!)
%%
%%  INPUT:
%%      W       Wavelet coefficients nf x nt x nch
%%      fw      Frequency vector in Hz

function [C, Wxy, sW, sWxy] = wave_acoh_tmp ( W, fw, t_smo, sc_smo )

%% Parameters
if nargin<4; sc_smo=0; end
if nargin<3; t_smo=3; end
if nargin<2; fw=logspace(-1,2,30); end
[nf, nt, nch] = size ( W );

%% Setting parameters
dt  = 1/2500;


%% Set variables
sW   = NaN(nf,nt,nch);
Wxy  = NaN(nf,nt,nch,nch);
sWxy = NaN(nf,nt,nch,nch);
C    = NaN(nf,nt,nch,nch);

%% Cross-spectra, Smoothing and Coherence (!)
for i=1:nch;
    sW(:,:,i)=smooth_mat(abs(W(:,:,i)).^2,dt,t_smo./fw,...
        'ScaleSmo', sc_smo);
end

for i=1:nch;
    for j=i+1:nch;
        Wxy(:,:,i,j)  = squeeze(W(:,:,i).*conj(W(:,:,j)));
        Wxy(:,:,j,i)  = conj(Wxy(:,:,i,j));
        sWxy(:,:,i,j) = smooth_mat(Wxy(:,:,i,j),dt,t_smo./fw,...
            'ScaleSmo', sc_smo);
        sWxy(:,:,j,i) = conj(sWxy(:,:,i,j));
        C(:,:,i,j)    = squeeze(abs(sWxy(:,:,i,j)).^2)./...
            squeeze(sW(:,:,i).*sW(:,:,j));
        C(:,:,j,i)    = C(:,:,i,j);
    end
end

clear dj fac n nf kl ke i j s0 j1

% %% Coherence can be plotted with
% i=1;
% j=1;
% plot_coh(ti,p,C(:,:,i,j),coi)