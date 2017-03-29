function PLV = wave_plv ( amp, phase, N )
%% WAVE_PLI Compute wavelet phase locking value from a single wavelet coefficient matrix
% 
%   PLV = WAVE_PLV( amp, phase, N )
% 
%	INPUT:
%           amp:     Amplitude of signal at frequency f1
%           phase:   Phase of signal in degrees(!) at frequency f2
%           N:       Number of data points used for average
%
%	OUTPUT:
%           PLV:     Phase locking value as a function of time
%
% Author: Michael von Papen
% Date: 17.05.16

%% Parameters
if nargin<3; N = 2457; end % 1s averaging interval
[namp, nt] = size ( amp );

%% Analytic signal
clear i
Z = zeros(namp,nt);
for j=1:namp
    Z(j,:) = amp(j,:).*exp(1i*phase);
end

%% Calculate PLV
As  = smooth_mat(amp,1,repmat(N,namp,1));
Z   = smooth_mat(Z  ,1,repmat(N,namp,1));
PLV = abs(Z)./As;