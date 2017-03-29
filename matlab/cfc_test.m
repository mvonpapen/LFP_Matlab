%% Phase amplitude coupling

function [M, M_raw] = pac ( experiment, patnum )


% Convert input
m = str2double(patnum);


% Parameters
w0  = 6;
Nf  = 60;
fnam= ['pac_m' patnum '_' experiment '.mat'];
f   = logspace(0,2.7,Nf); % 1-500 Hz


% Load and initialize
load(['/home/vpapenm/DATA/' experiment '.mat'])
M_raw = zeros(Nf,Nf,4,8);
M_sur = zeros(Nf,4,8);
M     = zeros(Nf,Nf,4,8);


% Calculate wavelet transform
[~,b] = size(LFP_OFF{m});
[~, Wx, cx] = procdata(LFP_OFF{m}, 'freq', f, 'w0', w0, 'art', ArtLFP_OFF{m});


% Start parallel pool
mypool = parpool(8);


for n=1:b
    
    W     = coi2nan(f, Wx(:,:,n), cx(:,n));
    phase = angle(W);

    for i=1:Nf
        
        PW  = abs(W(i,:)).^2;
        amp = repmat(zscore(PW),Nf,1);
        M_raw(:,i,n,m) = abs( nanmean(amp.*exp(1i*phase),2) );

        for k=1:200
            ind = randi(length(amp));
            A = amp(:,[ind:end, 1:ind-1]);
            M_sur(:,k,i) = abs( nanmean(A.*exp(1i*phase),2) );
        end

        M(:,i,n,m) = ( M_raw(:,i,n,m)-mean(M_sur(:,:,i),3) ) ...
            ./ std(M_sur(:,:,i),[],3);

    end

end

save(fnam, 'cfc', 'f', 'w0', 'phi', 'M')
delete(mypool);