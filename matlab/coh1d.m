%% Calculates time-averaged coherence from matrix C
%% alternatively calculates coherence for non-volume-conducted signals only

function [c1d, c1d_novc] = coh1d ( f, C, coi1, coi2, Wxy, phthres, sig )

nout = nargout;

[nf, nt, nx, ny] = size ( C );

c1d         = NaN (nf, nx, ny);
c1d_novc    = NaN (nf, nx, ny);

if nargin<7
    sig = 0.43;
end
if nargin<6
    phthres = 10;
end

for i=1:nx
    
    for j=1:ny
        
        c1d(:,i,j) = nanmean(coi2nan(f,C(:,:,i,j),min([coi1(:,i)'; coi2(:,j)'])),2);
        
        if nout>1
            x  = zeros(nf,nt);
            c  = C(:,:,i,j);
            ph = abs(angle(Wxy(:,:,i,j))/pi*180);
            vc = (ph<phthres | ph>180-phthres) & c>sig;
            x(~vc) = c(~vc);
            c1d_novc(:,i,j) = nanmean(coi2nan(f,x,min([coi1(:,i)'; coi2(:,j)'])),2);
        end
        
    end
    
end