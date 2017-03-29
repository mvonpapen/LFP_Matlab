%% Calculate power spectral density

function [P,f]=psdtf(x,dt,inc,y,win)

if nargin<5; win=0; end
if nargin<4 | y==0; flag=0; y=0; else flag=1; end
if nargin<3; inc=0; end

%% Pre-Whitening
if inc==1
    x=x(2:end)-x(1:end-1);
    if flag~=0
        y=y(2:end)-y(1:end-1);
    end
end

[m,n]=size(x);
if m>n
    x=x';
    y=y';
    [m,n]=size(x);
end

%% Construct frequency vector
f=[0:n-1]/(n*dt);

%% Calculate PSD
x=x-mean(x);
y=y-mean(y);
if win==1
    x=x.*hamming(n)';
    y=y.*hamming(n)';
end
xf=fft(x);
if flag~=0
    yf=fft(y);
end

for i=1:m;
    if mod(n,2)==0
        if flag~=0
            P(i,1)      =   1/n*dt*xf(:,1).*conj(yf(:,1));
            P(i,2:n/2)  =   2/n*dt*xf(:,2:n/2).*conj(yf(:,2:n/2));
        else
            P(i,1)      =   1/n*dt*abs(xf(:,1)).^2;
            P(i,2:n/2)  =   2/n*dt*abs(xf(:,2:n/2)).^2;
        end
    else
        if flag~=0
            P(i,1)          =   1/n*dt*xf(:,1).*conj(yf(:,1));
            P(i,2:(n+1)/2)  =   2/n*dt*xf(:,2:(n+1)/2).*conj(yf(:,2:(n+1)/2));
        else
            P(i,1)          =   1/n*dt*abs(xf(:,1)).^2;
            P(i,2:(n+1)/2)  =   2/n*dt*abs(xf(:,2:(n+1)/2)).^2;
        end
    end
end

%% Check Parseval
for i=1:m
    if flag~=0
        pt=cov(x,y);
        pt=abs(pt(2,1));
        pf=abs(sum(xf(i,:).*conj(yf(i,:)))/n)/n;
    else
        pt=sum(x(i,:).^2)/n;
        pf=sum(P)/n/dt;
    end
    
    if abs(pt-pf)/pt>0.01
        fprintf('Parseval check failed!\n');
    end
end

%% Use only frequencies up to Nyquist and take out f=0;
if mod(n,2)==0
    f=f(2:n/2);
    P=P(:,2:n/2);
else
    f=f(2:(n+1)/2);
    P=P(:,2:(n+1)/2);
end


%% Post-Darkening
if inc==1
    P=P./(4*sin(pi*dt*f).^2)';
end