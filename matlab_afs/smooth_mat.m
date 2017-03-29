function smoWave = smooth_mat(wave,dt,speriod,varargin)
%% SMOOTH_MAT Function used to smooth wavelets, input parameters are:
% 
%   INPUT:
%           wave:   (matrix from wavelet transform)
%           dt:     (time step)
%           speriod:(vector of periods used for smoothing, e.g. speriod=1./f. 
%                   Must have the same number of elements as the frequency 
%                   vector used for wavelet transformation. speriod is in
%                   unit of seconds!)
%
% 	OPTIONAL INPUT:
%           'scaleSmooth':  (1 on, or 0 off. Uses additional smoothing for 
%                           scales (box window after grinsted)
%           'smooth':       ('gauss' or 'box', default 'gauss', type of 
%                           smoothing window used for smoothing in time)
% 
%   OUTPUT:
%           smoWave: Smoothed wavelet coefficient matrix
% 
% Author: Michael von Papen
% 
% Date: 14.10.15


%% parameters
speriod=speriod(:);
[m, n] = size(wave);
smoWave = zeros(m,n);

args = struct('scaleSmooth',0,...
            'smooth','gauss');
args = parseArgs(varargin,args);


%% creating frequency vector for fft
f = [0:fix(n/2), -(fix((n-1)/2):-1:1)]'/(n*dt);

%% Take out zero frequency so that ifft does not produce imaginary values
tag = 0;
if speriod(1) == Inf
    wave0 = wave(1,:);
    wave  = wave(2:end,:);
    speriod = speriod(2:end);
    tag     = 1;
end


%% actual smoothing for each type by convolution using fft
switch upper(args.smooth)
    case 'GAUSS'
        F = exp( - 2 * pi^2 * f.^2 * speriod'.^2  ); % = fft of normalized gaussian with std = speriod
        smoWave = ifft(F.*fft(wave'))';
    case 'BOX'
        F = sinc( f * speriod' ); % = fft of normalized boxcar with width = speriod
        smoWave = ifft(F.*fft(real(wave')))'+1i*ifft(F.*fft(imag(wave')))';
end

%% scale smoothing after grinsted et al. 
ns = args.scaleSmooth;
if ns ~= 0
    for i = 1+ns:length(speriod)-ns
        smoWave(i,:) = nanmean(smoWave(i-ns:i+ns,:));
    end
end

%% Insert zero frequency value
if tag == 1
    smoWave = [wave0; smoWave];
end