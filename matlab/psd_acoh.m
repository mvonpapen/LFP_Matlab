%% Calculate PSD for coherent and incoherent parts
%
%
% Output:
%   Pcoh    PSD corresponding to coherent parts of C
%   Pinc    PSD corresponding to incoherent parts of C
%   Pvc     PSD corresponding to coherent parts of C with |dPhi|<phase_thresh

function [Pcoh, Pinc, Pvc, Pnvc] = psd_acoh ( f, W, C, coi, sig, Wxy, usenan, phase_thresh, zero_only )

if nargin<9; zero_only = false; end
if nargin<7; usenan    = false; end
if nargin<6 && nargout>2
    error('Cannot compute Pvc without Phase matrix!')
end
if nargin<6
    warning('Pcoh also contains volume conducted signals!')
end
if nargin<4
    sig=0.43;
end

[nf, nt, nx] = size ( W );
[nf2, nt2, nxw, ny] = size ( C );

if nf~=nf2 || nt~=nt2 || nx~=nxw || nx~=ny
    error('Sizes of W and C do not match!')
end
clear nf2 nt2 nxw ny

nout = nargout;

% Set cone of influence
if nargin<3
    coi=ones(nt,nx)*9e9;
end

% Set phase threshold to determine volume conduction PSD
if nargin<8
    phase_thresh = 15;
end

% Pre-allocate variables
Pcoh=NaN(nf,nx,nx);
Pinc=NaN(nf,nx,nx);
if nargout>2
    Pvc=NaN(nf,nx,nx);
end
if nout>3
    Pnvc=NaN(nf,nx,nx);
end

if nargin<6
    Ph = phase_thresh+1;
else
    Ph = abs(angle(Wxy)/pi*180);
end

%% Begin loop
for i=1:nx
    
    Xw = W(:,:,i);
    Xc = coi(:,i)';
    
    for j=1:nx
        
        if i==j
            continue
        end
        
        % find indices for (in)coherent coefficients
        if ~zero_only
            coh = C(:,:,i,j)>sig & Ph(:,:,i,j)>=phase_thresh ...
                & Ph(:,:,i,j)<=180-phase_thresh;
        else
            coh = C(:,:,i,j)>sig & Ph(:,:,i,j)>=phase_thresh;
        end
        inc = C(:,:,i,j)<=sig;
               
        % calculate PSD for coherent coefficients
        if ~usenan
            Wi = zeros(nf,nt);
        else
            Wi = NaN(nf,nt);
        end
        tmp = Xw;
        Wi(coh) = tmp(coh);
        Pcoh(:,i,j) = psdw( Wi, min([Xc; coi(:,j)']), f );
        
        % calculate PSD for incoherent coefficients
        if ~usenan
            Wi = zeros(nf,nt);
        else
            Wi = NaN(nf,nt);
        end
        tmp = Xw;
        Wi(inc) = tmp(inc);
        Pinc(:,i,j) = psdw( Wi, min([Xc; coi(:,j)']), f );
        
        % Optional: calculate PSD for coherent coefficients with dPhi~0
        % corresponding to volume conduction
        if nout>2
            vc = ~coh & ~inc;
            if ~usenan
                Wi = zeros(nf,nt);
            else
                Wi = NaN(nf,nt);
            end
            tmp = Xw;
            Wi(vc) = tmp(vc);
            Pvc(:,i,j) = psdw( Wi, min([Xc; coi(:,j)']), f ); 
        end
        if nout>3
            if ~zero_only
                vc = C(:,:,i,j)>sig & (Ph(:,:,i,j)<phase_thresh ...
                    | Ph(:,:,i,j)>180-phase_thresh);
            else
                vc = C(:,:,i,j)>sig & Ph(:,:,i,j)<phase_thresh;
            end
            Wi = zeros(nf,nt);
            tmp = Xw;
            Wi(~vc) = tmp(~vc);
            Pnvc(:,i,j) = psdw( Wi, min([Xc; coi(:,j)']), f );
        end            
        
    end
    
end
        