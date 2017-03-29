function coi = gencoi ( n, dt, tSmo )
%% GENCOI Generates Cone of influence after Torrence & Compo for 
% a Morlet wavelet with w0=6
%
%   IN:
%       n:      length of t or x
%       dt:     time step (default=1/2500)
%       tSmo:   ????TODO?????
% 
% Author: Michael von Papen
% 
% Date: 14.10.15
if nargin < 3; tSmo = 1; end
if nargin < 2; dt = 1/2456; end

coi = 1.03/sqrt(2)/tSmo*dt*...
      [1e-5,1:(n+1)/2-1,fliplr(1:n/2-1),1e-5];
  
end