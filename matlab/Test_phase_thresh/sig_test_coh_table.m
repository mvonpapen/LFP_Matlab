%% Test significance of coherence with power law noise


%% Parameter
w0    = [2 4 6 8 10 12];
nsig  = [1:10];
kappa = 1;
nt    = 2.^15;
dt    = 1/2500;
f     = 20;
ntrial= 500;
bin   = 0:0.005:1;

thres = ones( length(w0), length(nsig) );

%% Compute coherence
for j=1:length(w0)
    W0    = w0(j);
    fourier_factor = 4*pi/(W0+sqrt(2+W0^2));
    scale = 1 ./ (fourier_factor * f);
    coi   = fourier_factor/sqrt(2)*dt*...
        [1E-5,1:((nt+1)/2-1),fliplr((1:(nt/2-1))),1E-5];
    for k=1:length(nsig)
        NS   = nsig(k);
        tmp2 = zeros( nt, ntrial );
        tmp3 = zeros( nt, ntrial );
        parfor i=1:ntrial
            % noise
            x = powlawnoise(nt,kappa);
            y = powlawnoise(nt,kappa);
            % Wavelet transform
            W         = zeros(1, nt, 2);
            W(1,:,1)  = waveletlin( x, dt, f, 0, 'MORLET', W0 );
            W(1,:,2)  = waveletlin( y, dt, f, 0, 'MORLET', W0 );
            tmp       = wave_cohere( W, scale, NS, 1, dt );
            tmp2(:,i) = coi2nan(f, tmp(1,:,1,2), coi/NS);
        end
        c   = tmp2(:);
        H   = hist(c, bin);
        H   = H/sum(H);
        id = find(cumsum(H)>0.99,1,'first');
        if ~isempty(id)
            thres(j,k) = bin(id);
        end
        clear tmp2 tmp W c id x y NS
    end
end
clear i j k W0