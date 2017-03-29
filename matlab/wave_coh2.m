%% Compute wavelet coherencY (!)
%%
%%  INPUT:
%%      Wx   Wavelet coefficients of x (nf x nt x nch)
%%      Wy   Wavelet coefficients of y (nf x nt x nch)
%%      fw      Frequency vector in Hz
%%      dt      Time increment

function [C, Wxy, sWx, sWy, sWxy] = wave_coh2 ( Wx, Wy, fw, dt, t_smo, sc_smo )

%% Parameters
if nargin<6; sc_smo=0; end
if nargin<5; t_smo=1; end
if nargin<4; dt=1/250; end
if nargin<3; fw=logspace(-1,2,30); end
[m,n,nx]=size(Wx);
[m,n,ny]=size(Wy);

%% Setting parameters
tmp=fw(:)/1.033043648; %freq -> 1/s
fac=tmp(:,ones(1,n)); %factor 1/s to multiply with W (this is needed if smoothing is done in freq.domain as well)
clear tmp


%% Set variables
sWE=NaN(m,n,ny);
sWL=NaN(m,n,nx);
WLE=NaN(m,n,nx,ny);
sWLE=NaN(m,n,nx,ny);
C=NaN(m,n,nx,ny);

%% Cross-spectra, Smoothing and Coherency (!)
for i=1:ny;
    sWy(:,:,i)=sqrt(smooth_mat(fac.*abs(Wy(:,:,i)).^2,dt,t_smo./fw,...
        'ScaleSmo', sc_smo));
end

for i=1:nx;
    sWx(:,:,i)=sqrt(smooth_mat(fac.*squeeze(abs(Wx(:,:,i))).^2,dt,t_smo./fw,...
        'ScaleSmo', sc_smo)); 
    for j=1:ny;
        Wxy(:,:,i,j)=squeeze(Wx(:,:,i).*conj(Wy(:,:,j))); 
        sWxy(:,:,i,j)=smooth_mat(fac.*Wxy(:,:,i,j),dt,t_smo./fw,...
        'ScaleSmo', sc_smo);
        C(:,:,i,j)=squeeze(sWxy(:,:,i,j))./...
            squeeze(sWx(:,:,i).*sWy(:,:,j)); 
    end
end

clear dj fac n m kl ke i j s0 j1

% %% Coherence can be plotted with
% i=1;
% j=1;
% plot_coh(ti,p,C(:,:,i,j),coi)