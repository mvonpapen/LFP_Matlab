%% Test significance of coherence with power law noise for single frequency

% sig = 0.555 for two-sided test of ImCoh (P=0.995)

%% Parameter
w0    = 12;
nsig  = 6;
nt    = 2.^15;
dt    = 1/1000;
f     = 20;
ntrial= 1000;
bin   = -1:0.005:1;

H = zeros( length(bin), length(nsig) );
thres = ones( length(nsig), 1 );

%% Compute coherence
for j=1:length(nsig)
    fourier_factor = 4*pi/(w0+sqrt(2+w0^2));
    scale = 1 / (fourier_factor * f);
    tmp2 = zeros( nt, ntrial );
    tmp3 = zeros( ntrial, 1  );
    NT   = nt;
    coi  = fourier_factor/sqrt(2)*dt*...
        [1E-5,1:((NT+1)/2-1),fliplr((1:(NT/2-1))),1E-5];
    fprintf('Begin with loop %u/%u.\n', j, length(nsig));
    parfor i=1:ntrial
        % White noise
        x = powlawnoise(NT,1);
        y = powlawnoise(NT,1);
        % Wavelet transform
        W           = waveletlin( x, dt, f, 0, 'MORLET', w0 );
        W(1,:,2)    = waveletlin( y, dt, f, 0, 'MORLET', w0 );
        [~,~,~,~,tmp] = wave_cohere( W, scale, nsig(j), 1, dt );
        tmp2(:,i)   = squeeze( coi2nan(f, tmp(:,:,1,2), coi) );
        tmp3(i)     = nanmean(tmp2(:,i));
    end
    % For single frequency
    c      = tmp2;
    c      = c(:);
    H(:,j) = hist(c, bin);
    H(:,j) = H(:,j)/sum(H(:,j));
    id     = find(cumsum(H(:,j))>0.99,1,'first');
    if ~isempty(id)
        thres(j) = bin(id);
    end
    clear tmp2 W c id
end
clear tmp tmp2 W x y