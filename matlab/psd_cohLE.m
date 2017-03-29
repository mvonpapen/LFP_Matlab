%% Calculate PSD for coherent and incoherent parts of LFP-EMG data
%
%
% Output:
%   Pcoh    PSD corresponding to coherent parts of CLE
%   Pinc    PSD corresponding to incoherent parts of CLE

function [PLcoh, PEcoh, PLinc, PEinc] = psd_cohLE ( f, WL, WE, CLE, coiL, coiE, sig, WLE, usenan )

if nargin<9; usenan=false; end
if nargin<8
    error('Cannot compute Pvc without Phase matrix!')
end
if nargin<7
    sig=0.43;
end

[~, ~, nx] = size ( WL );
[nf, nt, ny] = size ( WE );
[nf2, nt2, nxw, nyw] = size ( CLE );

if nf~=nf2 || nt~=nt2 || nx~=nxw || ny~=nyw
    error('Sizes of W and C do not match!')
end
clear nf2 nt2 nxw nyw

% Pre-allocate variables
PLcoh=NaN(nf,nx,ny);
PLinc=NaN(nf,nx,ny);
PEcoh=NaN(nf,ny,nx);
PEinc=NaN(nf,ny,nx);

Ph = angle(WLE)/pi*180;


%% Begin loop for LFP
for i=1:nx
    
    wl = WL(:,:,i);
    cl = coiL(:,i)';
    
    for j=1:ny
        
        % find indices for (in)coherent coefficients
        coh = CLE(:,:,i,j)>sig;
               
        % calculate PSD for coherent coefficients
        if ~usenan
            Wi = zeros(nf,nt);
        else
            Wi = NaN(nf,nt);
        end
        tmp = wl;
        Wi(coh) = tmp(coh);
        PLcoh(:,i,j) = psdw( Wi, min([cl; coiE(:,j)']), f );
        
        % calculate PSD for incoherent coefficients
        if ~usenan
            Wi = zeros(nf,nt);
        else
            Wi = NaN(nf,nt);
        end
        tmp = wl;
        Wi(~coh) = tmp(~coh);
        PLinc(:,i,j) = psdw( Wi, min([cl; coiE(:,j)']), f );           
        
    end
    
end


%% Begin loop for EMG
for i=1:ny
    
    we = WE(:,:,i);
    ce = coiE(:,i)';
    
    for j=1:nx
        
        % find indices for (in)coherent coefficients
        coh = CLE(:,:,j,i)>sig;
               
        % calculate PSD for coherent coefficients
        if ~usenan
            Wi = zeros(nf,nt);
        else
            Wi = NaN(nf,nt);
        end
        tmp = we;
        Wi(coh) = tmp(coh);
        PEcoh(:,i,j) = psdw( Wi, min([ce; coiL(:,j)']), f );
        
        % calculate PSD for incoherent coefficients
        if ~usenan
            Wi = zeros(nf,nt);
        else
            Wi = NaN(nf,nt);
        end
        tmp = we;
        Wi(~coh) = tmp(~coh);
        PEinc(:,i,j) = psdw( Wi, min([ce; coiL(:,j)']), f );           
        
    end
    
end
        