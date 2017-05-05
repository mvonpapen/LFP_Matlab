%% Test significance of coherence with power law noise

%% Parameter
w0    = 12;
nsig  = 6;
nt    = 2.^15;
dt    = 1/1000;
f     = 20;
ntrial= 1000;
bin   = 0:0.005:1;

% Preallocate
H    = zeros( length(bin), length(nsig) );
p01  = zeros( length(nsig), 1 );
p05  = zeros( length(nsig), 1 );
tmp2 = zeros( 1, nt, ntrial );


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
    W           = waveletlin( x, dt, f, 0, 'MORLET', w0 );
    W(1,:,2)    = waveletlin( y, dt, f, 0, 'MORLET', w0 );
    tmp         = wave_cohere( W, scale, nsig, 1, dt );
    tmp2(:,:,i) = coi2nan(f, tmp(:,:,1,2), coi/nsig);
end

% Histogram
H   = hist(tmp2(:), bin);
H   = H/sum(H);
id1 = find(cumsum(H)>0.99,1,'first');
id5 = find(cumsum(H)>0.95,1,'first');
p01 = bin(id1);
p05 = bin(id5);


%% Save
save('sig_PCC_ns6_w12_nl3.mat', 'H', 'bin', 'p01', 'p05', ...
    'id1', 'id5', 'f', 'dt', 'w0', 'nsig', 'ntrial', 'nt')