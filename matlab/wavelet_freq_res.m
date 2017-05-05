%% Analyze frequency resolution of Morlet wavelet

% Parameters
dt = 1/2500;
Nt = 1e4;
f0 = 1:100;
w0 = 6:20;

% Wave number range
wk = 1:fix(Nt/2);
wk = wk.*((2.*pi)/(Nt*dt));
wk = [0., wk, -wk(fix((Nt-1)/2):-1:1)];
f  = wk/2/pi;

% Preallocate
sigma   = zeros(length(w0),length(f0));
dfexp   = zeros(length(w0),length(f0));
rel_sig = zeros(1, length(w0));

% Run loop
for j = 1:length(w0)
    
    for i = 1:length(f0)
        
        % Wavelet in freq. domain
        scale = (w0(j)+sqrt(2+w0(j)^2))/4/pi ./ f0(i);
        Psi   = exp( -(scale'.*wk - w0(j)).^2/2 ); % see Torrence & Compo, 1998
        Psi2  = Psi.^2; %relates to power

        % Gauss fit
        [sigma(j,i), mu] = gaussfit(f,Psi2,f0(i)/10,f0(i));
        
        % Look for 1/e
        e1 = find(Psi2>exp(-1),1,'first')-1;
        e2 = find(Psi2>exp(-1),1,'last')-1;
        dfexp(j,i) = mean( [f0(i)-f(e1) f(e2)-f0(i)] ) / f0(i);
        
    end
    
    rel_sig(j) = mean(sigma(j,:)./f0);
    
end