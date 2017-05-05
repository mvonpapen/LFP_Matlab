%% Function to compute PLI from wavelet cross spectrum

function PLI = phase_lag_index ( Wxy, scale, nsig, dt, varargin )

% input imaginary part of coherency matrix for wPLI [Vinck et al., 2011]
args = struct('wPLI', []);
        
args = parseArgs(varargin,args);
ImCo = args.wPLI;

if nargin<4
    dt = 1/2456;
end
if nargin<3
    nsig = 6;
end


% Extract phases
Ph = angle(Wxy)/pi*180;

% Average sign of phases over time
if isempty(ImCo)
    Ph_avg = smooth_mat( sign(Ph), dt, nsig*scale );
else
    Ph_avg_num   = smooth_mat( abs(ImCo).*sign(ImCo), dt, nsig*scale );
    Ph_avg_denum = smooth_mat( abs(ImCo), dt, nsig*scale );
    Ph_avg       = Ph_avg_num./Ph_avg_denum;
end

PLI = abs(Ph_avg);
