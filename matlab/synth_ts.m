function [x, y] = synth_ts( t, nl, A, dp )

if nargin<4
    dp = 30;
end
if nargin<3
    A = 1e3;
end
if nargin<2
    nl = 3;
end


%% Parameters
alpha = 10;
beta  = 20;
gamma = 50;
dt    = 1/2500;
nt    = length(t);

%% Synthetic time series
n1 = nl * powlawnoise(nt,1)';
n2 = nl * powlawnoise(nt,1)';
xa = [zeros(nt/5,1); sin(2*pi*alpha*(1:3*nt/5)'*dt); zeros(nt/5,1)];
xb = [sin(2*pi*beta*(1+(1:nt/2)'/nt*0.25).*(1:nt/2)'*dt); ...
      sin(2*pi*beta*(1+(nt/2+1:nt)'/nt*0.25).*(nt/2+1:nt)'*dt-dp)];
xg = [sin(2*pi*gamma*t(1:nt/2)'+dp); sin(2*pi*gamma*t(nt/2+1:end)')];
x  = xa + xb + xg + n1;
ya = [zeros(nt/5,1); sin(2*pi*alpha*(1:3*nt/5)'*dt+dp); zeros(nt/5,1)];
yb = sin(2*pi*beta*(1+t'/(nt*dt)*0.25).*t');
yg = [sin(2*pi*gamma*t(1:nt/2)'); sin(2*pi*gamma*t(nt/2+1:end)')];
y  = ya + yb + yg + n2;
x  = x / A;
y  = y / A;