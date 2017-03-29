%% Reconstruct time series from wavelet transform
%
%
%
% Input:
%
%   W       wavelet coefficient matrix obtained from wavelet.m with
%           size(W) = [nf, nt, nch], nf=#freq, nt=#time, nch=#channels
%   f       vector of frequencies
%   ind     vector of indices defining part of W to use for reconstruction

function [ts, P] = wave_recstr ( W, f, ind )

if nargin<3;
    ind=[];
end


%% Basic parameters
mother      = 'morlet';
foufac      = 1.033043648;
param       = 6;
dt          = 1/2500;
[nf, nt, nch]= size (W);


%% Define wavelet matrix
if ~isempty(ind)
    Wtmp        = zeros (nf, nt, nch);
    Wtmp(ind)   = W(ind);
    W           = Wtmp;
    clear Wtmp
end


%% Reconstruct time series
if strcmp(mother, 'morlet')
    
    % Simple reconstruction according to equation (11) in TC1998
    dj  = scale4wavelet ( f );
    s   = repmat( 1./f(:)/foufac, 1, nt);
    for i = 1:nch
        ts(:,i) = dj * sqrt(dt) / 0.776 / pi^(1/4) ...
                    * sum( real( W(:,:,i) ) ./ sqrt(s) );
    end
    
else
    
    %....construct wavenumber array used in transform [Eqn(5)] (see wavelet.m)
    k = [1:fix(nt/2)];
    k = k.*((2.*pi)/(nt*dt));
    k = [0., k, -k(fix((nt-1)/2):-1:1)];
    for i = 1:nch
        ts(:,i) = invcwt (W(:,:,i), mother, foufac./f, param, k);
    end
    
end



%% Calculate PSD
if nargout>1
    coi = gencoi (nt, dt);
    for i = 1:nch
        P(:,i) = mygws (abs( W(:,:,i) ).^2, dt, coi, f);
    end
end