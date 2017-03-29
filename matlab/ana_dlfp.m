%% Calculate coherence and phase locking value between cannels of x
%% and conduct principal value decomposition


function [Coh, Wxy, plv, pca] = ana_dlfp ( x, f, arts, nondiag )

if nargin<4
    nondiag=0;
end



%% Size of dLFP (nch x nt x nexp)
[nt, nch] = size ( x );
if nch>nt
    x=x';
    [nt, nch] = size ( x );
end

if nargin<3
    for i=1:nch
        arts{i}=[];
    end
end


%% Set variables
nf=length(f);
t_smo=3; % number of scales used to smooth
sc_smo=0; % number of neighboring scales to smooth in freq-domain
W=NaN(nf,nt,nch);
coi=NaN(nt,nch);
% PLFP=NaN(nf,5,ns);
% f=logspace(0,2.7,nf);
[dj,s0,j1] = scale4wavelet(f);    % Compute scales according to f
dt=1/2500;


%% Wavelet Transform
for i=1:nch
    [W(:,:,i),p,s,coi(:,i)] = wavelet(x(:,i),dt,0,dj,s0,j1);
    coi(:,i) = artcoi(coi(:,i),arts{i});
%     Pw(:,i) = mygws(abs(squeeze(W(:,:,i))).^2,dt,coi(:,i)',f);    
end
    
%% Coherence
% Coh=NaN(nf,nt,nch,nch);
if nondiag==1
    [Coh, Wxy] = wave_coh(W(:,:,1), W(:,:,2:end), f, dt, t_smo, sc_smo);
% for i=1:nch
%     Coh(:,:,i,i) = 1;
%     Coh(:,:,i,1) = tmp;
%     Coh(:,:,1,i) = tmp;
% end
% clear tmp
else
    [Coh, Wxy] = wave_coh(W, W, f, dt, t_smo, sc_smo);
end

%% Phase lockin value
if nargout>2
    plv = NaN(nf,nch,nch);
    for ii=1:nch
        for j=1:nch
            if ii==j
                plv(:,ii,j)=1;
                continue
            end
            for k=1:nf
                phase=phases(f(k),Wxy(k,:,ii,j),minmax(f),...
                    'coi',min([coi(:,ii)';coi(:,j)']));
    %             [h, bin]=hist(phase/pi*180,[-170:20:170]);
    %             h=h/sum(h);
                jj=~isnan(phase);
                plv(k,ii,j)=abs(sum(exp(sqrt(-1)*phase(jj))))/length(phase(jj));
            end
        end
    end
end

%% Principal component analysis
if nargout>3
    pca=princomp(x);
end