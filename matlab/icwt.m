%% Recontruction of time series from wavelet transform using inverse filter
%% method as described in Torrence & Compo (1998), equation (11), for a
%% Morlet wavelet

function x = icwt ( Wreal, f, dt )

if nargin<3
    dt=1/2500;
end

dj = scale4wavelet( f );

Cdelta = 0.776;
Psi0   = pi^(1/4);

s = repmat ( 1.03*f(:), 1, length(Wreal) );

x = dj * sqrt(dt) / Cdelta / Psi0 * sum( Wreal./sqrt(s) );