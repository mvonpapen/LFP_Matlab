%% Finds coefficients of significant coherence and non-zero phase
% difference. These can be interpreted as being caused by local coherent
% subpopulations in contrast to volume conduction.

% Input:
%   Coh     Matrix of coefficients (nf x nt x nch x nch)
%   Ph      Matrix of phase differences (nf x nt x nch x nch)
%
% Output:
%   iVC     Indices for coefficients in a (!) nf x nt (!) matrix


function iLC = indloco ( Coh, Ph )

%% Set threshold
sig_thresh  =0.5;
phase_thresh=10;

%% Check sizes
[a1,b1,c1,d1]=size(Coh);
[a2,b2,c2,d2]=size(Ph);
if a1~=a2 || b1~=b2 || c1~=c2 || d1~=d2
    error('Matrices not of equal size!')
end

%% Find indices
iLC  = 1:a1*b1;

for i=1:c1
    for j=1:c1
        if i==j
            continue
        end
        sig  = find( Coh(:,:,i,j)>=sig_thresh );
        phas = find( abs(Ph(:,:,i,j))*180/pi>=phase_thresh );
        lc   = intersect( sig, phas );
        iLC  = intersect( iLC, lc );
    end
end