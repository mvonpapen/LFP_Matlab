% Construct four dimensional spectral matrix from pair of time series x and y.

function S = spectral_matrix ( x, y, fnorm, w0, nsig )

N  = length(x);
Nf = length(fnorm);
scale_norm = (w0+sqrt(2+w0^2))/4/pi ./ fnorm;

% Wavelet transformation
W        = zeros(Nf, N, 2);
W(:,:,1) = waveletlin(x(:), 1, fnorm, 1, 'MORLET', w0);
W(:,:,2) = waveletlin(y(:), 1, fnorm, 1, 'MORLET', w0);


% Construct 4D spectral matrix
S = zeros(2,2,Nf,N);
S(1,1,:,:) = 2*abs(W(:,:,1)).^2;
S(2,2,:,:) = 2*abs(W(:,:,2)).^2;
S(1,2,:,:) = 2*W(:,:,1).*conj(W(:,:,2));
S(2,1,:,:) = conj(S(1,2,:,:));


% Temporal averaging
for i=1:2
    for j=1:2
        S(i,j,:,:) = smooth_mat(squeeze(S(i,j,:,:)), 1, scale_norm*nsig);
    end
end

% Extrapolate power at f=0
for j=1:2;
    for k=1:2;
        Sav = squeeze(mean(S(j,k,2:end,:),4));
        Sf0 = real(interp1(fnorm(2:end), Sav, 0, 'pchip'));
        if j==k && Sf0<=0
            Sf0 = real(interp1(fnorm(2:end), Sav, 0, 'linear'));
            Sf0 = max([0 Sf0]);
        end
        S(j,k,1,:) = Sf0;
    end
end