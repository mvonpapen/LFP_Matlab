%% Calculate pairwise wavelet coherence and saves as mat-file
%
%
% INPUT:
%   X       Nt x Nch - matrix, containing Nch time series of length Nt
%   f       frequency vector for which coherence should be estimated
%   w0      wavelet paremeter
%   nsig    averaging width to determine coherence
%   ds      downsampling factor (can usually be applied as coherence
%           estimate involves time averaging)
%   dt      sampling period of time series
%   tag     string for filename 'rs_[tag].mat'
%
%
% OUTPUT (contents of saved mat-file)
%   W       wavelet transform


function wave_coh_ts ( X, f, w0, nsig, ds, dt, tag, avg )

    if nargin<8
        avg = false;
    end
    
    %% Parameters
    Nf          = length(f);
    sig         = sig_coh_thresh(w0, nsig);
    [Nt, Nch]   = size(X);
    scale       = (w0+sqrt(2+w0^2))/4/pi ./ f;



    %% Wavelet transform
    W   = zeros (Nf, Nt, Nch);
    for i = 1:Nch
        [W(:,:,i),~,~,coi] = waveletlin(X(:,i),dt,f,1,'MORLET',w0);
    end
    if avg
        C = wave_cohere ( W, scale, nsig, ds );
        C = squeeze( nanmean(C, 2) );
    else
        [C, Wxy, PW] = wave_cohere ( W, scale, nsig, ds );
    end
    coi = coi(1:ds:end)/nsig;
    


    %% Save Variables
    if avg
        save(['rs_' tag '.mat'], 'C', 'f', 'w0', 'nsig', 'ds', ...
            'dt', 'tag', 'sig', 'Nf', 'scale', 'Nch', 'Nt')
    else
        save(['rs_' tag '.mat'], 'C', 'Wxy', 'W', 'coi', 'f', 'w0', 'nsig', 'ds', ...
            'dt', 'tag', 'sig', 'Nf', 'scale', 'Nch', 'Nt', 'PW')
    end
  