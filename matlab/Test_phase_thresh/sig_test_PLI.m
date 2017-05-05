%% Test significance of coherence with power law noise for single frequency

%% Parameter
w0    = 12;
nsig  = 6;
nt    = 2.^15;
dt    = 1/2500;
f     = 20;
ntrial= 1000;
bin   = 0:0.005:1;
wPLI  = false;
fnam = ['sig_PLI_ns' mat2str(nsig) '_w' mat2str(w0) '_v2.mat'];


%% Compute coherence
fourier_factor = 4*pi/(w0+sqrt(2+w0^2));
scale = 1 / (fourier_factor * f);
tmp  = zeros( nt, ntrial );
NT   = nt;
coi  = fourier_factor/sqrt(2)*dt*...
    [1E-5,1:((NT+1)/2-1),fliplr((1:(NT/2-1))),1E-5];


for i=1:ntrial
    % White noise
    x = powlawnoise(NT,1);
    y = powlawnoise(NT,1);
    % Wavelet transform
    W           = waveletlin( x, dt, f, 1, 'MORLET', w0 );
    W(1,:,2)    = waveletlin( y, dt, f, 1, 'MORLET', w0 );
    Wxy         = squeeze(W(1,:,1).*conj(W(1,:,2)));
    PLI = phase_lag_index( imag(Wxy), scale, nsig, dt, 'wPLI', wPLI );
    tmp(:,i) = squeeze( coi2nan(f, PLI, coi) );
end
% Histogram
H   = hist(tmp(:), bin);
H   = H/sum(H);
id1 = find(cumsum(H)>0.99,1,'first');
id5 = find(cumsum(H)>0.95,1,'first');
p01 = bin(id1);
p05 = bin(id5);

save(fnam, 'H', 'bin', 'p01', 'p05', ...
    'id1', 'id5', 'PLI', 'f', 'dt', 'w0', 'nsig', 'ntrial', 'nt')