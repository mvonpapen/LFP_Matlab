%% WAVELET_FREQ  1D Wavelet transform for specific frequency vector freq
%
%   [wave,period,scale,coi] = WAVELET_FREQ(Y,dt,freq,pad,mother,param)
%
%   Computes the wavelet transform of the vector Y (length N),
%   with sampling rate DT at given frequencies FREQ.
%
%   By default, the Morlet wavelet (k0=6) is used.
%   The wavelet basis is normalized to have total energy=1 at all scales.
%
%
% INPUTS:
%
%    Y = the time series of length N.
%    DT = amount of time between each Y value, i.e. the sampling time.
%    FREQ = frequencies for wavelet analysis
%
% OUTPUTS:
%
%    WAVE is the WAVELET transform of Y. This is a complex array
%    of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
%    ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
%    The WAVELET power spectrum is ABS(WAVE)^2.
%    Its units are sigma^2 (the time series variance).
%
%
% OPTIONAL INPUTS:
% 
% *** Note *** setting any of the following to -1 will cause the default
%               value to be used.
%
%    PAD = if set to 1 (default is 0), pad time series with enough zeroes to get
%         N up to the next higher power of 2. This prevents wraparound
%         from the end of the time series to the beginning, and also
%         speeds up the FFT's used to do the wavelet transform.
%         This will not eliminate all edge effects (see COI below).
%
%    FREQ = Frequency vector.
%
%    MOTHER = the mother wavelet function.
%             The choices are 'MORLET', 'PAUL', or 'DOG'
%
%    PARAM = the mother wavelet parameter.
%            For 'MORLET' this is k0 (wavenumber), default is 6.
%            For 'PAUL' this is m (order), default is 4.
%            For 'DOG' this is m (m-th derivative), default is 2.
%
%
% OPTIONAL OUTPUTS:
%
%    PERIOD = the vector of "Fourier" periods (in time units) that corresponds
%           to the SCALEs.
%
%    SCALE = the vector of scale indices, given by S0*2^(j*DJ), j=0...J1
%            where J1+1 is the total # of scales.
%
%    COI = if specified, then return the Cone-of-Influence, which is a vector
%        of N points that contains the maximum period of useful information
%        at that particular time.
%        Periods greater than this are subject to edge effects.
%        This can be used to plot COI lines on a contour plot by doing:
%
%              contour(time,log(period),log(power))
%              plot(time,log(coi),'k')
%
%----------------------------------------------------------------------------
%   This script is a modified version of Torrence & Compo's wavelet.m. It
%   has been modified by Felix Gerick and Michael von Papen
%   Institute of Geophysics & Meteorology, University of Cologne, Germany
%   Date: 7/24/15
%
%
%   Copyright (C) 1995-1998, Christopher Torrence and Gilbert P. Compo
%   University of Colorado, Program in Atmospheric and Oceanic Sciences.
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.
%
% Notice: Please acknowledge the use of this program in any publications:
%   ``Wavelet software was provided by C. Torrence and G. Compo,
%     and is available at URL: http://paos.colorado.edu/research/wavelets/''.
%
% Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
%            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
%
% Please send a copy of such publications to either C. Torrence or G. Compo:
%  Dr. Christopher Torrence               Dr. Gilbert P. Compo
%  Advanced Study Program                 NOAA/CIRES Climate Diagnostics Center
%  National Center for Atmos. Research    Campus Box 216
%  P.O. Box 3000                          University of Colorado at Boulder
%  Boulder CO 80307--3000, USA.           Boulder CO 80309-0216, USA.
%  E-mail: torrence@ucar.edu              E-mail: gpc@cdc.noaa.gov
%----------------------------------------------------------------------------

function [wave,period,scale,coi,k] = waveletlin(Y,dt,freq,pad,mother,param)


if (nargin < 6), param = -1; end
if (nargin < 5), mother = -1; end
if (nargin < 4), pad = 0; end
if (nargin < 3), freq = -1; end
if (nargin < 2)
	error('Must input a vector Y and sampling time DT')
end

n1 = length(Y);

if (mother == -1), mother = 'MORLET'; end

%....construct SCALE array.
if (freq == -1)
    s0=2*dt;
    dj = 1./4.;
    J1=fix((log(n1*dt/s0)/log(2))/dj);
    scale = s0*2.^((0:J1)*dj);
else
    switch upper(mother)
        case 'MORLET'
            if (param == -1), param = 6.; end
            fourier_factor = (4*pi)/(param + sqrt(2 + param^2));
        case 'PAUL'
            if (param == -1), param = 4.; end
            fourier_factor = 4*pi/(2*param+1);
        case 'DOG'
            if (param == -1), param = 2.; end
            fourier_factor = 2*pi*sqrt(2./(2*param+1));
    end
    scale=1./(fourier_factor*freq);
end


%....construct time series to analyze, pad if necessary
x(1:n1) = Y - mean(Y);
if (pad == 1)
	base2 = fix(log(n1)/log(2) + 0.4999);   % power of 2 nearest to N
	x = [x,zeros(1,2^(base2+1)-n1)];
end
n = length(x);

%....construct wavenumber array used in transform [Eqn(5)]
k = [1:fix(n/2)];
k = k.*((2.*pi)/(n*dt));
k = [0., k, -k(fix((n-1)/2):-1:1)];

%....compute FFT of the (padded) time series
f = fft(x);    % [Eqn(3)]

%....construct WAVE arrays
wave = zeros(length(scale),n);  % define the wavelet array
wave = wave + 1i*wave;  % make it complex

% loop through all scales and compute transform
for a1 = 1:length(scale)
	[daughter,fourier_factor,coi]=wave_bases(mother,k,scale(a1),param);	
	wave(a1,:) = ifft(f.*daughter);  % wavelet transform[Eqn(4)]
end

period = 1./freq;
coi = coi*dt*[1E-5,1:((n1+1)/2-1),fliplr((1:(n1/2-1))),1E-5];  % COI [Sec.3g]
wave = wave(:,1:n1);  % get rid of padding before returning

return