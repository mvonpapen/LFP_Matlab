%% Calculates time-averaged coherence from matrix C, where diagonal elements
%% correspond to autocoherences and C is symmetric
%% alternatively calculates coherence for non-volume-conducted signals only

function [c1d, c1d_novc] = acoh1d_uptri ( f, C, coi, Wxy, phthres, sig )

[nf, nt, nx, ny] = size ( C );
if nx~=ny
    error('Coherence matrix C is not symmetric!')
end

if nargin<6; sig=0.43; end
if nargin<5
    phthres     = 10;
end
c1d         = NaN (nf, nx, nx);
c1d_novc    = NaN (nf, nx, nx);

nout = nargout;

for i=1:nx
    
    for j=i+1:nx
        
        c1d(:,i,j) = nanmean(coi2nan(f,C(:,:,i,j),min([coi(:,i)'; coi(:,j)'])),2);
        
        if nout>1
            x  = zeros(nf,nt);
            c  = C(:,:,i,j);
            ph = abs(angle(Wxy(:,:,i,j))/pi*180);
            vc = (ph<phthres | ph>180-phthres) & c>sig;
            x(~vc) = c(~vc);
            c1d_novc(:,i,j) = nanmean(coi2nan(f,x,min([coi(:,i)'; coi(:,j)'])),2);
        end
        
    end
    
end