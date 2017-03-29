function Pw = psdw(W,coi,f,dt)
%% PSDW Calculate power spectral density from Wavelet coefficients
% 
%   Pw = PSDW(W,coi,f,dt,dj)
% 
%   INPUT:
%           W:   Wavelet coefficient matrix
%           coi: Cone of influence vector
%           f:   Wavelet frequency vector
%           dt:  Time increment
% 
%   OUTPUT:
%           Pw: Power spectral density vector of length of f
% 
% Author: Michael von Papen
% 
% Date: 15.10.15

if nargin<4; dt = 1/2456; end

if ~isreal(W)
    W = 2*dt*abs(W).^2;
end


%% Respect COI
if nargin >= 3
    W = coi2nan(f,W,coi);
end


%% Calculate PSD
Pw = nanmean(W,2);

end