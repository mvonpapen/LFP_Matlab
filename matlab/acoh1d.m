%% Calculates time-averaged coherence from matrix C, where diagonal elements
%% correspond to autocoherences and C is symmetric
%% alternatively calculates coherence for non-volume-conducted signals only

function [c1d, c1d_novc] = acoh1d ( f, C, coi, Wxy, phthres, ref )

[nf, nt, nx, ny] = size ( C );
if nx~=ny
    error('Coherence matrix C is not symmetric!')
end

if nargin<5
    phthres     = 10;
end
c1d         = NaN (nf, nx, nx);
c1d_novc    = NaN (nf, nx, nx);


for i=1:nx
    
    for j=i+1:nx
        
        if nargin < 6
            ref = j;
        end
        
        c1d(:,i,j) = nanmean(coi2nan(f,C(:,:,i,j),min([coi(:,i)'; coi(:,j)'])),2);
        c1d(:,j,i) = c1d(:,i,j);
        
        if nargout>1
            x  = zeros(nf,nt);
            c  = C(:,:,i,j);
            ph = abs(angle(Wxy(:,:,i,ref))/pi*180);
            vc = (ph<phthres | ph>180-phthres) & c>0.495;
            x(~vc) = c(~vc);
            c1d_novc(:,i,j) = nanmean(coi2nan(f,x,min([coi(:,i)'; coi(:,j)'])),2);
            c1d_novc(:,j,i) = c1d_novc(:,i,j);
        end
        
    end
    
end