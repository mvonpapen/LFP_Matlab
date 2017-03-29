function smoWave=smoothed(wave,dt,scale,dj, varargin)

% function used to smooth wavelets, input parameters are:
%   wave (coefficient matrix from wavelet transform)
%   dt (time step)
%   scale (vector of scales used for wavelet transform)
%   dj (scale step)
% optional parameters are:
%   'scaleSmooth' (1 on, or 0 off. Uses additional smoothing for scales (box 
%   window after grinsted)
%   'smooth' ('gauss' or 'box', default 'gauss', type of smoothingwindow used
%   for smoothing in time)



%% parameters

n=size(wave,2);
m=size(wave,1);
smoWave=zeros(size(wave));

args=struct('scaleSmooth',1,'smooth','gauss');
args=parseArgs(varargin,args);

%% creating w_k vektor after torrence/compo 1998 eq. (5)
k = 1:fix(n/2);
k = k.*((2.*pi)/(n*dt));
k = [0., k, -k(fix((n-1)/2):-1:1)];

%% actual smoothing for each type by convolution using fft

switch upper(args.smooth)
    case 'GAUSS'
        for i=1:m
            F=scale(i)*exp(-scale(i)^2.*k.^2/2);
            smoWave(i,:)=ifft(F.*fft(wave(i,:)));
        end
        smoWave=real(smoWave);
    case 'BOX'
        for i=1:m   
            k(1)=0.001;
            F=sin(4*k.*scale(i))./(2*k.*scale(i));
%           F=F./abs(F);
            smoWave(i,:)=ifft(F.*fft(wave(i,:)));
        end
end

% %% scale smoothing after grinsted et al. 
% 
% if args.scaleSmooth==1
%     dj0=0.6;
%     dj0steps=dj0/(dj*2);
%     boxCar=[mod(dj0steps,1); ones(2 * round(dj0steps)-1,1); ...
%         mod(dj0steps,1)] / (2*round(dj0steps)-1+2*mod(dj0steps,1));
%     smoWave=conv2(smoWave,boxCar,'same');
% end

end
