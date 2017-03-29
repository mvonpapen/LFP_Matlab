%% Phase amplitude coupling

function [M_off, M_on, M_raw_off, M_raw_on] = pac ( experiment, patnum )

%% Setup


% % Start parallel pool
% mypool = parpool(8);


% Convert input
m = str2double(patnum);


% Parameters
w0   = 12;
Nf   = 60;
nsig = 6;
phthres = 15;
% fnam= ['/home/vpapenm/DATA/pac_m' patnum '_' experiment '.mat'];
fnam= ['pac_pat' patnum '_' experiment '.mat'];
f   = logspace(0,2.7,Nf); % 1-500 Hz


% Load and initialize
% load(['/home/vpapenm/DATA/' experiment '.mat'])
load(['data_' experiment '_v3.mat'])
M_raw_off = zeros(Nf,Nf,4);
M_off     = zeros(Nf,Nf,4);
M_raw_on  = zeros(Nf,Nf,4);
M_inc_off = zeros(Nf,Nf,4,4);
M_inc_on  = zeros(Nf,Nf,4,4);
M_coh_off = zeros(Nf,Nf,4,4);
M_coh_on  = zeros(Nf,Nf,4,4);
M_vc_off  = zeros(Nf,Nf,4,4);
M_vc_on   = zeros(Nf,Nf,4,4);




%% OFF

% Calculate wavelet transform
[~,b] = size(LFP_OFF{m});
[~, Wx, cx] = procdata(LFP_OFF{m}, 'freq', f, 'w0', w0, 'art', ArtLFP_OFF{m});
scale   = (w0+sqrt(2+w0^2))/4/pi ./ f;

% Coherences
sig=sig_coh_thresh(w0,nsig);
[C, Wxy] = wave_cohere ( Wx, scale, nsig, 1 );
inc = C<sig;
coh = C>sig & angle(Wxy)/pi*180>phthres;
vc  = C>sig & angle(Wxy)/pi*180<=phthres;

% Run loop over channels
for n=1:b
    
    W = coi2nan(f, Wx(:,:,n), cx(:,n));

    % Run loop over frequencies
    for i=1:Nf
        
        phase = repmat( angle(W(i,:)), Nf, 1 );
        PW  = abs(W).^2;
        amp = (PW-nanmean(PW,2))./nanstd(PW,[],2);
        M_raw_off(i,:,n) = abs( nanmean(amp.*exp(1i*phase),2) );

        M_sur = zeros(Nf,200);
        for k=1:200
            ind = randi(length(amp));
            A = amp(:,[ind:end, 1:ind-1]);
            M_sur(:,k) = abs( nanmean(A.*exp(1i*phase),2) );
        end
        M_off(i,:,n) = ( M_raw_off(i,:,n)-mean(M_sur,2) ) ./ std(M_sur,[],2);
        
        % Loop for inc, coh, vc
        for j=setxor(n,1:length(LFP_OFF))

            % Inc
            ti = inc(i,:,n,j);
            phase = repmat( angle(W(i,ti)), Nf, 1 );
            PW    = abs(W(:,ti)).^2;
            amp = (PW-nanmean(PW,2))./nanstd(PW,[],2);
            M_raw = abs( nanmean(amp.*exp(1i*phase),2) );
            M_sur = zeros(Nf,200);
            for k=1:200
                ind = randi(length(amp));
                A = amp(:,[ind:end, 1:ind-1]);
                M_sur(:,k) = abs( nanmean(A.*exp(1i*phase),2) );
            end
            M_inc_off(i,:,n,j) = ( M_raw-mean(M_sur,2) ) ./ std(M_sur,[],2);

            % Coh
            ti = coh(i,:,n,j);
            phase = repmat( angle(W(i,ti)), Nf, 1 );
            PW    = abs(W(:,ti)).^2;
            amp = (PW-nanmean(PW,2))./nanstd(PW,[],2);
            M_raw = abs( nanmean(amp.*exp(1i*phase),2) );
            M_sur = zeros(Nf,200);
            for k=1:200
                ind = randi(length(amp));
                A = amp(:,[ind:end, 1:ind-1]);
                M_sur(:,k) = abs( nanmean(A.*exp(1i*phase),2) );
            end
            M_coh_off(i,:,n,j) = ( M_raw-mean(M_sur,2) ) ./ std(M_sur,[],2);

            % VC
            ti = vc(i,:,n,j);
            phase = repmat( angle(W(i,ti)), Nf, 1 );
            PW    = abs(W(:,ti)).^2;
            amp = (PW-nanmean(PW,2))./nanstd(PW,[],2);
            M_raw = abs( nanmean(amp.*exp(1i*phase),2) );
            M_sur = zeros(Nf,200);
            for k=1:200
                ind = randi(length(amp));
                A = amp(:,[ind:end, 1:ind-1]);
                M_sur(:,k) = abs( nanmean(A.*exp(1i*phase),2) );
            end
            M_vc_off(i,:,n,j) = ( M_raw-mean(M_sur,2) ) ./ std(M_sur,[],2);
            
        end

    end
    
end




%% ON

% Calculate wavelet transform
[~,b] = size(LFP_ON{m});
[~, Wx, cx] = procdata(LFP_ON{m}, 'freq', f, 'w0', w0, 'art', ArtLFP_ON{m});

% Coherences
sig=sig_coh_thresh(w0,nsig);
[C, Wxy] = wave_cohere ( Wx, scale, nsig, 1 );
inc = C<sig;
coh = C>sig & angle(Wxy)/pi*180>phthres;
vc  = C>sig & angle(Wxy)/pi*180<=phthres;


% Run loop over channels
for n=1:b
    
    W = coi2nan(f, Wx(:,:,n), cx(:,n));

    % Run loop over frequencies
    for i=1:Nf
        
        phase = repmat( angle(W(i,:)), Nf, 1 );
        PW  = abs(W).^2;
        amp = (PW-nanmean(PW,2))./nanstd(PW,[],2);
        M_raw_on(i,:,n) = abs( nanmean(amp.*exp(1i*phase),2) );

        M_sur = zeros(Nf,200);
        for k=1:200
            ind = randi(length(amp));
            A = amp(:,[ind:end, 1:ind-1]);
            M_sur(:,k) = abs( nanmean(A.*exp(1i*phase),2) );
        end
        M_on(i,:,n) = ( M_raw_on(i,:,n)-mean(M_sur,2) ) ./ std(M_sur,[],2);
        
        % Loop for inc, coh, vc
        for j=setxor(n,1:length(LFP_OFF))

            % Inc
            ti = inc(i,:,n,j);
            phase = repmat( angle(W(i,ti)), Nf, 1 );
            PW    = abs(W(:,ti)).^2;
            amp = (PW-nanmean(PW,2))./nanstd(PW,[],2);
            M_raw = abs( nanmean(amp.*exp(1i*phase),2) );
            M_sur = zeros(Nf,200);
            for k=1:200
                ind = randi(length(amp));
                A = amp(:,[ind:end, 1:ind-1]);
                M_sur(:,k) = abs( nanmean(A.*exp(1i*phase),2) );
            end
            M_inc_on(i,:,n,j) = ( M_raw-mean(M_sur,2) ) ./ std(M_sur,[],2);

            % Coh
            ti = coh(i,:,n,j);
            phase = repmat( angle(W(i,ti)), Nf, 1 );
            PW    = abs(W(:,ti)).^2;
            amp = (PW-nanmean(PW,2))./nanstd(PW,[],2);
            M_raw = abs( nanmean(amp.*exp(1i*phase),2) );
            M_sur = zeros(Nf,200);
            for k=1:200
                ind = randi(length(amp));
                A = amp(:,[ind:end, 1:ind-1]);
                M_sur(:,k) = abs( nanmean(A.*exp(1i*phase),2) );
            end
            M_coh_on(i,:,n,j) = ( M_raw-mean(M_sur,2) ) ./ std(M_sur,[],2);

            % VC
            ti = vc(i,:,n,j);
            phase = repmat( angle(W(i,ti)), Nf, 1 );
            PW    = abs(W(:,ti)).^2;
            amp = (PW-nanmean(PW,2))./nanstd(PW,[],2);
            M_raw = abs( nanmean(amp.*exp(1i*phase),2) );
            M_sur = zeros(Nf,200);
            for k=1:200
                ind = randi(length(amp));
                A = amp(:,[ind:end, 1:ind-1]);
                M_sur(:,k) = abs( nanmean(A.*exp(1i*phase),2) );
            end
            M_vc_on(i,:,n,j) = ( M_raw-mean(M_sur,2) ) ./ std(M_sur,[],2);
            
        end

    end
    
end



%% Save and end

save(fnam, 'f', 'w0', 'M_raw_off', 'M_off', 'M_raw_on', 'M_on')
% delete(mypool);