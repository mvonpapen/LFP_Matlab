
function [xf, poles] = bfilt (x,pass,N,type)
%% BFILT Filter data x(t) with two-side Butterworth high pass filter according to SPM8
%   [xf, poles] = BFILT(x,pass,N,type)
%  OUTPUT:
%      xf:      filtered signal
%
%  INPUT:
%      x:       unfiltered signal (nchannels x ndata)
%      pass:    highpass normalized by Nyquist frequency (pass = f/Fn)
%      N:       filter order (default N=5)
%      type:    type of filter (high, stop, low, passband)
%
%  Author: Michael von Papen
%  Date:   16.09.2014

%% Set to defaults
if nargin<4 && numel(pass)==1; type='high'; end
if nargin<4 && numel(pass)==2 && pass(1)~=0; type='stop'; end
if nargin<4 && numel(pass)==2 && pass(1)==0;
    type='high';
    pass=pass(2);
end
if nargin<4 && numel(pass)==2 && pass(2)==0;
    type='low';
    pass=pass(1);
end
if nargin<3; N=5; end

%% Compute filter coefficients
switch numel(pass)
    case 1
        [B, A] = butter(N, pass, type);
    case 2
        [B, A] = butter(N, [min(pass) max(pass)], type);
    otherwise
        error(['Too much filter frequencies! NUMEL(pass)=' num2str(numel(pass))]);
end

%% Check for poles
poles = roots(A);
if any(abs(poles) >= 1)
  disp(poles);
  error(['Calculated filter coefficients have poles on or outside ' ...
        'the unit circle and will not be stable. Try a higher cutoff ' ...
        'frequency or a different type/order of filter.']);
end

%% Apply filter with automatic correction
xf = filtfilt(B, A, x')';