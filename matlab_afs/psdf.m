function [P,f] = psdf(x, dt, inc, y)
%%PSDF Calculate power spectral density using fast fourier transform
% 
%   [P,f] = PSDF(x, dt, inc, y)
%
%   INPUT:
%           x:   Data time series
%           dt:  Time increment
%           inc: Pre-whitening (incrementing) of time series on/off 
%                (default 0)
%           y:   second data time series. then a cross power spectral
%                density is computed
% 
%   OUTPUT:
%           P:   Power spectral density vector of length of frequencies 
%           f:   Frequency vector
% 
% Author: Michael von Papen
% 
% Date: 14.10.15

if nargin<4; y = 0; end
if nargin<3; inc = 0; end

%% Check sizes, transpone if necessary
[n, m] = size(x);
if m>n
    x = x';
    if y ~= 0
        y = y';
    end
    [n, m] = size(x);
end

%% Pre-Whitening
if inc == 1
    x = x(2:end,:)-x(1:end-1,:);
    if y ~= 0
        y = y(2:end,:)-y(1:end-1,:);
    end
    [n, m] = size(x);
end


%% Calculate PSD
% for i = 1:m
%     x(:,i) = x(:,i)-mean(x(:,i));
%     if y ~= 0
%         y(:,i) = y(:,i)-mean(y(:,i));
%     end
% end
xf = fft(x);
if y ~= 0
    yf = fft(y);
end

if mod(n,2) == 0
    if y ~= 0
        P = 2/n*dt*xf(1:n/2+1,:).*conj(yf(1:n/2+1,:));
    else
        P = 2/n*dt*abs(xf(1:n/2+1,:)).^2;
    end
    P(n/2+1,:) = P(n/2+1,:)/2; %This is the Nyquist frequency
else
    if y ~= 0
        P = 2/n*dt*xf(1:(n+1)/2,:).*conj(yf(1:(n+1)/2,:));
    else
        P = 2/n*dt*abs(xf(1:(n+1)/2,:)).^2;
    end
end
P(1,:) = P(1,:)/2; %This is the zero-frequency


%% Check Parseval
if y == 0
    for i = 1:m
        if y~= 0
            pt = sum(x(:,i))*sum(y(:,i))/n;
        else
            pt = sum(x(:,i).^2)/n;
        end
        pf = sum(abs(P(:,i)))/n/dt;
        if abs(pt-pf)>1e-9
            fprintf('Parseval check failed!\n');
        end
    end
end
    
%% Output is one-sided spectrum
if mod(n,2) == 0
    f = (0:n/2)/(n*dt);
    P = P(1:n/2+1,:);
else
    f = (0:(n+1)/2-1)/(n*dt);
    P = P(1:(n+1)/2,:);
end


%% Post-Darkening
if inc == 1
    for i = 1:m
        P(:,i) = P(:,i)./(4*sin(pi*dt*f).^2)';
    end
end