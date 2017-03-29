%% Calculates time-averaged coherence from multi-group matrix C without
%% taking into account volume-conducted signals with dPhi~0 and C>sig

function c1d_novc = coh_novc ( f, C, coi1, coi2, Wxy, phthres, ref )

[nf, nt, nx, ny] = size ( C );

c1d_novc    = NaN (nf, nx, ny);

if nargin<6
    phthres = 10;
end

for i=1:nx
    
    for j=1:ny
        
        if nargin < 7
            ref = j;
        end
        
        if nargout>1
            x  = zeros(nf,nt);
            c  = C(:,:,i,j);
            ph = abs(angle(Wxy(:,:,i,ref))/pi*180);
            vc = (ph<phthres | ph>180-phthres) & c>0.495;
            x(~vc) = c(~vc);
            c1d_novc(:,i,j) = nanmean(coi2nan(f,x,min([coi1(:,i)'; coi2(:,j)'])),2);
        end
        
    end
    
end