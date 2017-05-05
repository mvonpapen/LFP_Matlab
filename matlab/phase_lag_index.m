%% Function to compute PLI from wavelet cross spectrum

function PLI = phase_lag_index ( Im, scale, nsig, dt, varargin )

% input imaginary part of coherency matrix for wPLI [Vinck et al., 2011]
args = struct('wPLI', false);
        
args = parseArgs(varargin,args);
wPLI = args.wPLI;

if nargin<4
    dt = 1/2456;
end
if nargin<3
    nsig = 6;
end

% Average sign of phases over time
if ~wPLI
    Ph_avg = smooth_mat( sign(Im), dt, nsig*scale );
else
    % biased est
    Ph_avg_num   = smooth_mat( Im, dt, nsig*scale );
    Ph_avg_denom = smooth_mat( abs(Im), dt, nsig*scale );
    Ph_avg       = Ph_avg_num./Ph_avg_denom;
%     % debiased est (taken from fieldtrip
%     Im      = imag(Wxy);
%     outsum  = smooth_mat( Im, dt, nsig*scale );
%     outsumW = smooth_mat( abs(Im), dt, nsig*scale );
%     outssq  = smooth_mat( Im.^2, dt, nsig*scale, 'smooth','gauss_square' );
%     Ph_avg  = (outsum.^2 - outssq)./(outsumW.^2 - outssq);
end
  
PLI = abs(Ph_avg);