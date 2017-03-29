%% Calculates time-averaged coherence from single-group matrix C without
%% taking into account volume-conducted signals with dPhi~0 and C>sig.

function c1d_novc = acoh_novc ( f, C, coi, Wxy, phthres, ref )

[nf, nt, nx, ny] = size ( C );
if nx~=ny
    error('Coherence matrix C is not symmetric!')
end

if nargin<5
    phthres     = 10;
end
c1d_novc    = NaN (nf, nx, nx);


for i=1:nx
    
    for j=i+1:nx
        
        if nargin < 6
            ref = j;
        end
        
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