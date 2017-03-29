function [W, coi, P] = wave ( f, data, arts )

if nargin<3; arts=[]; end

m = length ( f );
[n,nch] = size ( data );
if nch>n
    data=data';
    [nch,n] = size ( data );
end
    
dt=1/2500;

Pw=NaN(m,nch);
coi=NaN(n,nch);
W=NaN(m,n,nch);

[dj,s0,j1] = scale4wavelet(1./f);    % Compute scales according to f
for i=1:nch
    [W(:,:,i),p,s,coi(:,i)] = wavelet(data(:,i),dt,0,dj,s0,j1);
    if ~isempty(arts)
        coi(:,i) = artcoi(coi(:,i),arts(i,:));
    end
    P(:,i) = mygws(abs(squeeze(W(:,:,i))).^2,dt,coi(:,i)',f);
end

if isempty(arts)
    coi=coi(:,1);
end