%% Compute wavelet magnitude squared coherencE (!)
%%
%%  INPUT:
%%      Wx      Wavelet coefficients of channel x (nf x nt x nchx)
%%      Wy      Wavelet coefficients of channel y (nf x nt x nchy)
%%      fw      Frequency vector in Hz
%%      dt      Time increment
%%
%%  OUTPUT:
%%      C       Magnitude squared Coherence (nf x nt x nchx x nchy)
%%      sWxy    Smoothed cross-spectra of Wx and Wy (nf x nt x nchx x nchy)
%%      sWi     Smoothed auto-spectra of Wi (nf x nt x nchi)
%%      Wxy     Cross-spectra of Wx and Wy (nf x nt x nchx x nchy), not smoothed


function [C, sWxy, sWx, sWy, Wxy] = wave_coh3 ( Wx, Wy, scale, nsig, ds )

%% Parameters
if nargin<5; ds=1; end
if nargin<4; nsig=3; end
[~,~,nx]=size(Wx);
[m,n,ny]=size(Wy);

%% Setting parameters
dt=1/2500;
scale = scale(:);


%% Set variables
sWx=NaN(m,n,nx);
sWy=NaN(m,n,ny);
Wxy=NaN(m,n,nx,ny);
sWxy=NaN(m,n,nx,ny);
if ds==1
    C    = NaN(m,n,nx,ny);
else
    n2  = ceil(n/ds);
    C   = NaN(m,n2,nx,ny);
end

%% Cross-spectra, Smoothing and Coherence (!)
for i=1:ny;
    sWy(:,:,i)=smooth_mat(abs(Wy(:,:,i)).^2,dt,nsig*scale);
end

for i=1:nx;
    sWx(:,:,i)=smooth_mat(squeeze(abs(Wx(:,:,i))).^2,dt,nsig*scale); 
    for j=1:ny;
        Wxy(:,:,i,j)=squeeze(Wx(:,:,i).*conj(Wy(:,:,j))); 
        sWxy(:,:,i,j)=smooth_mat(Wxy(:,:,i,j),dt,nsig*scale);
        if ds==1
            C(:,:,i,j)    = abs(sWxy(:,:,i,j)).^2 ./ (sWx(:,:,i).*sWy(:,:,j));
        end
    end
end

% Downsample data, take every dec'th value
if ds~=1
    sWx  = sWx(:,1:ds:end,:);
    sWy  = sWy(:,1:ds:end,:);
    sWxy = sWxy(:,1:ds:end,:,:);
    for i=1:nx
        sX = sWx(:,:,i);
        for j=1:ny
            C(:,:,i,j) = abs(sWxy(:,:,i,j)).^2 ./ (sX.*sWy(:,:,j));
        end
    end
end


clear dj fac n m kl ke i j s0 j1

% %% Coherence can be plotted with
% i=1;
% j=1;
% plot_coh(ti,p,C(:,:,i,j),coi)
