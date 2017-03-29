%% Morlet wavelet after Torrence & Compo, 1998, equation (1)

function [Psi, smooth_win] = wavelet_morlet( t, t0, w0, f, dt, nsig )

[n1, n2] = size(f);
if n1==1 && n2==1
    scale   = (w0+sqrt(2+w0^2))/4/pi ./ f;
    eta     = (t-t0) ./ scale;
else
    if n1==1 && n2>1
        f       = f';
        n1      = n2;
    end
    f       = repmat(f, 1, length(t));
    t       = repmat(t, n1, 1);
    scale   = (w0+sqrt(2+w0^2))/4/pi ./ f;
    eta     = (t-t0) ./ scale;
end



Psi0 = pi^(-1/4) * exp(1i*w0*eta-eta.^2/2);

Psi  = sqrt(dt./scale) .* Psi0;

if nargout>1
    smooth_win = abs( pi^(-1/4) * exp(1i*w0*eta/nsig-(eta/nsig).^2/2) );
    smooth_win = sqrt(dt./scale) .* smooth_win;
end