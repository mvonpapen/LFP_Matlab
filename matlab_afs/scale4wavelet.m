function [dj, s0, j1] = scale4wavelet ( f, fourierFactor )
%% SCALE4WAVELET Creates variables needed for Wavelet transformation (Morlet,k0=6)
% from given frequency vector
% 
%   INPUT:
%           f: Frequencies for wavelet transformation
%
%   OUTPUT: (needed for Torrence & Compo wavelet function)
%           DJ: the spacing between discrete scales. Default is 0.25.
%               A smaller # will give better scale resolution, 
%               but be slower to plot.
%
%           S0: the smallest scale of the wavelet.  Default is 2*DT.
%
%           J1: the # of scales minus one. Scales range from S0 up to 
%               S0*2^(J1*DJ),to give a total of (J1+1) scales.
%               Default is J1 = (LOG2(N DT/S0))/DJ.
%
%    Here, direction is reversed, so that wavelet transform is according to
%       increasing frequencies, not increasing periods p

f = f(:);
s0 = 1./min(f)/fourierFactor;
j1 = length(f)-1;
dj = log2(1/max(f)/fourierFactor/s0)/j1;