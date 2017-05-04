%% Function to compute PLI from wavelet cross spectrum

function PLI = phase_lag_index ( Wxy, scale, nsig, dt )

if nargin<4
    dt = 1/2456;
end
if nargin<3
    nsig = 6;
end


% Extract phases
Ph = angle(Wxy)/pi*180;

% Average sign of phases over time
Ph_avg = smooth_mat( sign(Ph), dt, nsig*scale );

PLI = abs(Ph_avg);