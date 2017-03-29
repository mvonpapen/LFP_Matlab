%% x=lsqcurvefit(@mygauss,[5e3 0 3],bin,h_lfp)

% Output: x(1) = Amplitude
%         x(2) = Mean
%         x(3) = Std

function F = mygauss(x,xdata)
F = x(1)*exp(-(xdata-x(2)).^2/2/x(3)^2);