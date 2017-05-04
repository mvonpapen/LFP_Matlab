%% Test significance of coherence with power law noise for single frequency

% sig = 0.555 for two-sided test of ImCoh (P=0.995)

%% Parameter
w0    = 12;
nsig  = 6;
nt    = 2.^15;
dt    = 1/1000;
f     = 20;
ntrial= 1000;
bin   = 0:0.005:1;

H   = zeros( length(bin), length(nsig) );
p01 = zeros( length(nsig), 1 );
p05 = zeros( length(nsig), 1 );

%% Compute coherence
for j=1:length(nsig)
    fourier_factor = 4*pi/(w0+sqrt(2+w0^2));
    scale = 1 / (fourier_factor * f);
    tmp  = zeros( nt, ntrial );
    NT   = nt;
    coi  = fourier_factor/sqrt(2)*dt*...
        [1E-5,1:((NT+1)/2-1),fliplr((1:(NT/2-1))),1E-5];
    fprintf('Begin with loop %u/%u.\n', j, length(nsig));
    
    for i=1:ntrial
        % White noise
        x = powlawnoise(NT,1);
        y = powlawnoise(NT,1);
        % Wavelet transform
        W           = waveletlin( x, dt, f, 0, 'MORLET', w0 );
        W(1,:,2)    = waveletlin( y, dt, f, 0, 'MORLET', w0 );
        Wxy         = squeeze(W(1,:,1).*conj(W(1,:,2)));
        PLI         = phase_lag_index( Wxy, scale, nsig );
        tmp(:,i)    = squeeze( coi2nan(f, PLI, coi) );
    end
    % Histogram
    H(:,j) = hist(tmp(:), bin);
    H(:,j) = H(:,j)/sum(H(:,j));
    id1    = find(cumsum(H(:,j))>0.99,1,'first');
    id5    = find(cumsum(H(:,j))>0.95,1,'first');
    p01(j) = bin(id1);
    p05(j) = bin(id5);
end

save('sig_PLI_ns6_w12_nl3.mat', 'H', 'bin', 'p01', 'p05', ...
    'id1', 'id5', 'PLI', 'f', 'dt', 'w0', 'nsig', 'ntrial', 'p01', 'nt')