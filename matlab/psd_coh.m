%% Calculate PSD for coherent and incoherent parts
%
%
% Output:
%   Pcoh    PSD corresponding to coherent parts of C
%   Pinc    PSD corresponding to incoherent parts of C
%   Pvc     PSD corresponding to coherent parts of C with |dPhi|<phase_thresh
%
% Ptotal = Pcoh + Pinc + Pvc

function [Pcoh, Pinc, Pvc] = psd_coh ( f, W, C, coi, sig, Wxy )

if nargin<6 && nargout>2
    error('Cannot compute Pvc without Phase matrix!')
end
if nargin<6
    warning('Pcoh also contains volume conducted signals!')
end
if nargin<4
    sig=0.495;
end

[nf, nt, nx] = size ( W );
[nf2, nt2, nxw, ny] = size ( C );

if nf~=nf2 || nt~=nt2 || nx~=nxw
    error('Sizes of W and C do not match!')
end
clear nf2 nt2 nxw

% Set cone of influence
if nargin<3
    coi=ones(nt,nx)*9e9;
end

% Set phase threshold to determine volume conduction PSD
phase_thresh = 10;

% Pre-allocate variables
Pcoh=NaN(nf,nx,ny);
Pinc=NaN(nf,nx,ny);
if nargout>2
    Pvc=NaN(nf,nx,ny);
end

if nargin<6
    Ph = phase_thresh+1;
else
    Ph = abs(angle(Wxy)/pi*180);
end

%% Begin loop
for i=1:nx
    
    for j=1:ny
        
        % find indices for (in)coherent coefficients
        coh = C(:,:,i,j)>sig & Ph(:,:,i,j)>=phase_thresh ...
            & Ph(:,:,i,j)<=180-phase_thresh;
        inc = C(:,:,i,j)<=sig;
               
        % calculate PSD for coherent coefficients
        Wi = zeros(nf,nt);
        tmp = W(:,:,i);
        Wi(coh) = tmp(coh);
        Pcoh(:,i,j) = mygws( Wi, min([coi(:,i)'; coi(:,j)']), f );
        
        % calculate PSD for incoherent coefficients
        Wi = zeros(nf,nt);
        tmp = W(:,:,i);
        Wi(inc) = tmp(inc);
        Pinc(:,i,j) = mygws( Wi, min([coi(:,i)'; coi(:,j)']), f );
        
        % Optional: calculate PSD for coherent coefficients with dPhi~0
        % corresponding to volume conduction
        if nargout>2
            vc = find( C(:,:,i,j)>sig & ... % without C>sig very conservative
                ( Ph(:,:,i,j)<phase_thresh | Ph(:,:,i,j)>180-phase_thresh ) );
            Wi = zeros(nf,nt);
            tmp = W(:,:,i);
            Wi(vc) = tmp(vc);
            Pvc(:,i,j) = mygws( Wi, min([coi(:,i)'; coi(:,j)']), f );
            
        end
        
    end
    
end
        