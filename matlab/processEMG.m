%%  Process EMG data for spectral analysis:
%%  1. 10 Hz High-pass filter
%%  2. Rectify signal
%%  3. Smooth signal
%%  4. Downsample signal
%%  5. Compute wavelet transform
%%  6. Compute Global wavelet spectra
%%  7. Compute Fourier spectra
%%
%%  Input:
%%      t       Equidistant time vector in s (n elements)
%%      EMG     EMG signal of k channels (n x k elements)
%%      hp      High pass filter frequency
%%      win     Window to use for smoothing
%%      ndt     New dt for downsampling
%%      fw      Log-equidistant frequency vector in Hz (m elements) to evaluate
%%
%%  Output:
%%      nt      New time vector with ndt in s (n2 elements)
%%      EMGout  Processed EMG data (n2 elements)
%%      W       Wavelet transform of EMG data (m x n2 x k elements)
%%      Pw      Global wavelet spectrum (m x k elements)
%%      ff      Frequency vector for Fourier spectrum
%%      Pf      Fourier-based power spectrum (m x k elements)
%%
%%  Date:   05.11.2014
%%  Author: M. von Papen
%%
%%

function [ nt, EMGout, W, Pw, ff, Pf, coi ] = processEMG ( t, EMG, hp, win, ndt, fw )

if nargin<6; fw=logspace(-1,2,100); end
if nargin<5; ndt=1/250; end
if nargin<4; win=gausswin(11)/sum(gausswin(11)); end
if nargin<3; hp=10; end


%% Parameters
[n, k] = size(EMG);
m=length(fw);
dt=t(2)-t(1);


%% Preset variables
EMGf=NaN(n,k);
W=NaN(m,n,k);
Pw=NaN(m,k);


%% 1. Filter time series with two-side butterworth filter order 5
for i=1:k
    EMGf (:,i) = bfilt ( EMG(:,i), hp*2*dt );
end


%% 2. & 3. Rectify and smooth signal
EMGo = abs(EMGf);
EMGo = convn(EMGo, win,'same');
clear EMGf


%% 4. Downsample signal with linear interpolation
nt=[min(t):ndt:max(t)];
EMGo = interp1(t, EMGo, nt);


%% 5. & 6. Wavelet Transform and global wavelet spectra
if nargout>2
    [dj,s0,j1] = scale4wavelet(1./fw);    % Compute scales according to fw
    for i=1:k
        [W(:,:,i),p,s,coi] = wavelet(EMGo(:,i),ndt,0,dj,s0,j1);
        Pw(:,i) = mygws(squeeze(W(:,:,i)),ndt,coi,p);
    end
end
clear dj s0 j1 i

%% 7. Fourier Spectra
if nargout>4
    [ff,Pf] = mypsd ( EMGo, ndt );
end