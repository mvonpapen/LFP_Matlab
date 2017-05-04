%% Test significance of coherence with power law noise


%% Parameter
w0    = 12;
nsig  = 6;
kappa = 1;
nt    = 2.^14;
dt    = 1/2500;
f     = [1 5 10 20 50 100];
ntrial= 1000;
bin   = 0:0.005:1;

H = zeros( length(bin), length(f) );
thres = ones( length(f), 1 );
% H = zeros( length(bin), length(f), length(nsig) );
% thres = ones( length(f), length(nsig) );

%% Compute coherence
% for j=1:length(nsig)
    fourier_factor = 4*pi/(w0+sqrt(2+w0^2));
    scale = 1 ./ (fourier_factor * f);
    tmp2 = zeros( length(f), nt, ntrial );
    tmp3 = zeros( length(f), nt, ntrial );
    NT   = nt;
    coi  = fourier_factor/sqrt(2)*dt*...
        [1E-5,1:((NT+1)/2-1),fliplr((1:(NT/2-1))),1E-5];
%     fprintf('Begin with loop %u/%u.\n', j, length(nsig));
    parfor i=1:ntrial
        % noise
        x = powlawnoise(NT,kappa);
        y = powlawnoise(NT,kappa);
        % Wavelet transform
        W           = waveletlin( x, dt, f, 0, 'MORLET', w0 );
        W(:,:,2)    = waveletlin( y, dt, f, 0, 'MORLET', w0 );
        tmp         = wave_cohere( W, scale, nsig, 1, dt );
        tmp2(:,:,i) = coi2nan(f, tmp(:,:,1,2), coi/nsig);
%         fprintf('\b\b\b\b\b\b\b%3u/%3u. ', i, ntrial);
    end
    for i=1:length(f)
        c        = tmp2(i,:,:);
        c        = c(:);
        H(:,i)   = hist(c, bin);
        H(:,i)   = H(:,i)/sum(H(:,i));
        id = find(cumsum(H(:,i))>0.99,1,'first');
        if ~isempty(id)
            thres(i) = bin(id);
        end
%         H(:,i,j)   = hist(c, bin);
%         H(:,i,j)   = H(:,i,j)/sum(H(:,i,j));
%         id = find(cumsum(H(:,i,j))>0.99,1,'first');
%         if ~isempty(id)
%             thres(i,j) = bin(id);
%         end
    end
    clear tmp2 W c id
% end
clear tmp tmp2 W x y


% %% Calculate pdf
% for j=1:length(nt)
% end
% clear c i j


% %% Find 1% and 5% significances
% sig = NaN( 2, length(f), length(t_smo) );
% for j=1:length(nt)
%     for i=1:length(f)
%         p1 = find( cumsum((bin(2)-bin(1))*H(:,i,j))>=0.99, 1, 'first' );
%         p5 = find( cumsum((bin(2)-bin(1))*H(:,i,j))>=0.95, 1, 'first' );
%         sig(:,i,j) = bin([p1 p5]);
%     end
% end
% clear p1 p5 i j
        

% %% Plot results
% ylimits = [4 8 16];
% for i=1:length(nt)
%     subplot(1,length(t_smo),i)
%     plot(bin, H(:,:,i))
%     hold all
%     plot(sig(1,1,i)*[1 1], [0 ylimits(i)], '--k')
%     ylim([0 ylimits(i)])
%     xlabel('Coherence')
%     if i==1
%         ylabel('p.d.f.')
%     end
%     title(['n_\sigma=' num2str(t_smo(i))])
%     if i==3
%         legend(['f = ' num2str(f(1)) ' Hz'], ['f = ' num2str(f(2)) ' Hz'], ...
%             ['f = ' num2str(f(3)) ' Hz'])
%     end
% end
% clear ylimits i