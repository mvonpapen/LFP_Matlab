%% Test phase distribution

%% Parameter
f     = 20;  %sets coherent frequency of phase variation
dp    = 0 / 180*pi;
nsig  = 2:10;
w0    = 6:15;
nl    = 1:5; %noise level

beta  = 20;

bin   = -180:180;
dt    = 1/2500;
nt    = 2^16;
t     = (0:nt-1)*dt;
T     = t(end);
ds    = 1;
A     = 1e-3; % amplitude correction
N     = 1000;


%% Preallocate
n1 = length(bin);
n2 = length(nsig);
n3 = length(w0);
n4 = length(nl);
H = zeros(n1, n2, n3, n4, N);

for i=1:n2
    
    ns = nsig(i);
    
    for j=1:n3
        
        scale = (w0(j)+sqrt(2+w0(j)^2))/4/pi ./ f;
        w     = w0(j);
        
        for k=1:n4
            
            Pn = nl(k);
            
            [i, j, k],

            for l=1:N

                %% Synthetic time series
                Pn1 = Pn * powlawnoise(nt,1)';
                Pn2 = Pn * powlawnoise(nt,1)';
                x  = sin(2*pi*beta*t') + Pn1;
                y  = sin(2*pi*beta*t') + Pn2;
                x  = x / A;
                y  = y / A;


                %% Spectral analysis
                [~,W,coi]   = procdata([x y], 'freq', f, 'w0', w);
                Wxy         = squeeze(W(:,:,1).*conj(W(:,:,2)));
                Wxy         = smooth_mat(Wxy, dt, ns*scale);
                Wxy         = Wxy(:,1:ds:end);
                coi         = coi(1:ds:end,1)/ns;
                Wxy         = coi2nan(f, Wxy, coi);

                %% Phase distribution
                Ph           = angle(Wxy)/pi*180;
                H(:,i,j,k,l) = hist(Ph, bin);
                H(:,i,j,k,l) = H(:,i,j,k,l)/trapz(bin,H(:,i,j,k,l));

            end
            
        end
        
    end
    
end

% % Plot results via
% % mseb(bin, mean(H,3)', std(H,0,3)');
% 
% %% Test for 1/e to find Phi_c
p10 = zeros(n2,n3,n4);
pexp = zeros(n2,n3,n4);
for i=1:n2
    for j=1:n3
        for k=1:n4
            X   = mean(squeeze(H(:,i,j,k,:)),2);
            tmp = bin(cumsum(X)>0.1 & cumsum(X)<0.9);
            p10(i,j,k)  = (tmp(end)-tmp(1))/2;
            pexp(i,j,k) = (bin(find(X/max(X)>exp(-1),1,'last')) ...
                           -bin(find(X/max(X)>exp(-1),1,'first')))/2;
        end
    end
end
clear i j k X tmp
save 'phase_dist_ns_w0_nl_v3.mat'
% db = bin(2)-bin(1);
% X = mean(H,3);
% for i=1:length(nsig)
%     tmp     = bin(cumsum(X(:,i))*db>0.1 & cumsum(X(:,i))*db<0.9);
%     Xi(i,1) = (tmp(end)-tmp(1))/2;
%     Xi(i,2) = bin(find(X(:,i)/max(X(:,i))>exp(-1),1,'first'));
% end
% Xi
% 
% clear i j db