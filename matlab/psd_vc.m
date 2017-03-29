%% This function calculates the total power and the volume conducted power
%% for a set of electrodes. Volume conduction is defined as near-zero phase
%% for all electrode combinations.

function [Pvc, Ptot] = psd_vc ( W, Wxy, phase_thresh, coi, f, C )

[nf, nt, nch]       = size ( W );
[nf2, nt2, nx, ny]  = size ( Wxy );
Ptot                = NaN(nf, nch);
Pvc                 = NaN(nf, nch);


if (nf~=nf2) || (nt~=nt2) || (nch~=nx) || (nch~=ny)
    error ('Matrices have different sizes!')
end

Ph = abs(angle(Wxy)/pi*180);


for i=1:nch
    %% Determine coefficients with near-zero phase
    ni = setdiff(1:nch, i);
    ind = false( nch-1,numel(Ph(:,:,1,1)) );
    for j=1:nch-1
        if nargin>5
            ph      = Ph(:,:,i,ni(j));
            c       = C(:,:,i,ni(j));
            x       = (ph<phase_thresh | ph>180-phase_thresh) & c>0.495;
        else
            ph      = Ph(:,:,i,ni(j));
            x       = ph<phase_thresh | ph>180-phase_thresh;
        end
        ind(j,:)= x(:);
    end
    if nch>2
        VC = all(ind);
    else
        VC = ind;
    end
    
    %% Calculate PSD
    Wvc         = W(:,:,i);
    Wvc(~VC)     = 0;
    Pvc(:,i)    = mygws ( Wvc,      coi(:,i), f);
    if nargout>1
        Ptot(:,i)   = mygws ( W(:,:,i), coi(:,i), f);
    end
end