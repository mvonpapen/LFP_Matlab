%% Calculate power spectral density

function [Pw, fw, E]=mygws(W,coi,f,dt,dj)

if nargin<5; dj=0; end
if nargin<4; dt=1/2500; end


%% Calculate energy
if dj ~= 0
    [m n]=size(W);
    S=repmat(1./f,n,1)'; %periods -> scales
    E=sum(sum(abs(W).^2./S,2)*dj*dt/n/0.776);
end

%% Respect COI
if nargin>=3
    W=coi2nan(f,W,coi);
end


%% Calculate PSD
Pw=nanmean(2*dt*abs(W).^2,2);
