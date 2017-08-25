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

function [ts, P, varWT] = wave_recstr ( W, f, w0, dt, ind )

if nargin<5
    ind=[];
end
if nargin<4
    dt = 1/2456;
end


%% Basic parameters
mother      = 'morlet';
foufac      = 4*pi/(w0+sqrt(2+w0^2));
switch w0
    case 6
        Cdelta = 0.776;
    case 12
        Cdelta = 0.3804;
end
if w0~=6 && w0~=12
    warning('Cdelta not known for this value of w0')
end
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
    dj  = scale4wavelet ( f, foufac );
    s   = repmat( 1./f(:)/foufac, 1, nt);
    for i = 1:nch
        ts(:,i)  = -dj * sqrt(dt) / Cdelta / pi^(-1/4) ...
                    * sum( real( W(:,:,i) ) ./ sqrt(s) );
        varWT(i) = -dj * dt / Cdelta / nt ...
                    * sum(sum( abs( W(:,:,i) ).^2 ./ s ));
    end
    
else
    
    %....construct wavenumber array used in transform [Eqn(5)] (see wavelet.m)
    k = [1:fix(nt/2)];
    k = k.*((2.*pi)/(nt*dt));
    k = [0., k, -k(fix((nt-1)/2):-1:1)];
    for i = 1:nch
        ts(:,i) = invcwt (W(:,:,i), mother, foufac./f, w0, k);
    end
    
end



%% Calculate PSD
if nargout>1
    coi = gencoi (nt, dt);
    for i = 1:nch
        P(:,i) = mygws (W(:,:,i), coi, f, dt);
    end
end