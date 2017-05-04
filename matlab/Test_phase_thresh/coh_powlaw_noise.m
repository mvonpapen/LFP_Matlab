%% Calculate coherence for noise with spectral index kappa at frequency f
%% Determine p.d.f. and 1% threshold

function coh_powlaw_noise ( f, kappa )


f     = str2double(f);
kappa = str2double(kappa);

fdir = '/home/vpapenm/matlab/';
mypool = parpool(8);


%% Parameter
w0    = 6:18;
nsig  = 1:10;
nt    = 2.^16;
dt    = 1/2456;
ntrial= 1e3;
bin   = 0:0.005:1;

H     = zeros( length(bin), length(w0), length(nsig) );
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
        
        H(:,j,k) = hist(tmp2(:), bin);
        H(:,j,k) = H(:,j,k)/sum(H(:,j,k));
        id       = find(cumsum(H(:,j,k))>0.99,1,'first');
        if ~isempty(id)
            thres(j,k) = bin(id);
        end
        
        clear tmp2 tmp W id x y NS
        
    end
    
end

save([fdir 'coh_noise_f' mat2str(f) '_si' mat2str(kappa) '.mat'], ...
    'H', 'f', 'nsig', 'w0', 'dt', 'ntrial', 'coi', 'scale', 'nt', 'kappa', ...
    'thres', 'bin');

delete(mypool)