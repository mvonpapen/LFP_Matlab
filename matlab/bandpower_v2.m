%% Check changes on power for rest OFF/ON

%% Load LFP data during rest
clear all
[t, dat, rigor, CH, art, ndat] = load_data ( 'experiment', 'Ruhe', ...
    'activity', 1, 'commontime', 1, 'noEMG', 1);

% Sites with rest ON OFF
i0   = [2  6 11 18 26 32 36 41 47 65];
iend = [5 10 17 25 31 67 40 46 55 66];

% % Sites with fist ON OFF
% i0   = [1 3 6 8 10 11 14];
% iend = [2 4 7 9 19 12 15];

% % Sites with hold ON OFF
% i0   = [1 3 5 7  9 11 12 15 19 21];
% iend = [2 4 6 8 10 23 13 16 20 22];

nex  = length(i0);

% Reference (stimulating?) electrodes: P02-A2, P03-C1, P04-C1, PECH-C1, 
% P06-C1, P07-C1, SCUD-M1, P08-P2, VAJU-C1, P05-P2
% channel numbers disregarding channels with activity=0
re = [2 1 1 1 1 1 1 2 1 2];

%% Parameters
nf=50;
f=logspace(0,2.6,nf);
delta=find(f>=1 & f<4);
theta=find(f>=4 & f<7);
alpha=find(f>=7 & f<13);
beta=find(f>=13 & f<30);
beta1=find(f>=13 & f<20);
beta2=find(f>=20 & f<30);
gamma=find(f>=30 & f<70);
gamma2=find(f>=60 & f<80);
HF=find(f>=250 & f<350);


%% Calculate PSD during ON and OFF
nch=0; for i=1:length(CH), nch=max([nch length(CH{i})]); end
Poff = NaN(nf,nch,length(i0));
Pon  = NaN(nf,nch,length(iend));
Pcoh_off = NaN(nf,nch,length(i0));
Pcoh_on  = NaN(nf,nch,length(iend));
Pinc_off = NaN(nf,nch,length(i0));
Pinc_on  = NaN(nf,nch,length(iend));
Pvc_off  = NaN(nf,nch,length(i0));
Pvc_on   = NaN(nf,nch,length(iend));
% P=NaN(nf,nch,nex);
for n=1:length(i0) %nex
%     nLFP=length(CH{n});
    nLFP(n)=length(CH{i0(n)});
    [x,Woff,coi_off,Poff(:,1:nLFP(n),n)]=procdata(dat{i0(n)}, 'freq', f, ...
        'filter', [45 55; 90 110], 'art', art{i0(n)});
    [x,Won,coi_on,Pon(:,1:nLFP(n),n)]=procdata(dat{iend(n)}, 'freq', f, ...
        'filter', [45 55; 90 110], 'art', art{iend(n)});
%     [x,x,x,P(:,1:nLFP,n)]=procdata(dat{n}, 'freq', f, 'filter', [45 55; 90 110], ...
%         'art', art{n});

    %% Calculate coherent spectra of stim-elec w.r.t. other STN-LFP
    if nLFP(n)>1
        [C_off,Wxy_off] = wave_acoh ( Woff,f );
        [C_on ,Wxy_on ] = wave_acoh ( Won ,f );
        [Pcoh_on(:,1:nLFP(n),1:nLFP(n),n) , Pinc_on(:,1:nLFP(n),1:nLFP(n),n) , ...
            Pvc_on(:,1:nLFP(n),1:nLFP(n),n) ] = psd_coh( f, Won , C_on , ...
            coi_on , 0.5, angle(Wxy_on) );
        [Pcoh_off(:,1:nLFP(n),1:nLFP(n),n), Pinc_off(:,1:nLFP(n),1:nLFP(n),n), ...
            Pvc_off(:,1:nLFP(n),1:nLFP(n),n)] = psd_coh( f, Woff, C_off, ...
            coi_off, 0.5, angle(Wxy_off));
    end
end

%% Determine ratio
Prel = NaN (nf, nch-1, nex);
for n=1:length(i0)
    x = setxor ( 1:nLFP(n), re(n) );
    Prel(:,1:length(x),n) = squeeze( Pcoh_off(:,re(n),x,n)./Pcoh_on(:,re(n),x,n) );    
end
Prel(isinf(Prel))=NaN;
Prel(Prel==0)=NaN;

%% Determine averages in freq. bands
tmp = nanmean(Prel(HF,:,:)); Prelband(:,9) = tmp(:);
tmp = nanmean(Prel(gamma2,:,:)); Prelband(:,8) = tmp(:);
tmp = nanmean(Prel(gamma,:,:)); Prelband(:,7) = tmp(:);
tmp = nanmean(Prel(beta2,:,:)); Prelband(:,6) = tmp(:);
tmp = nanmean(Prel(beta1,:,:)); Prelband(:,5) = tmp(:);
tmp = nanmean(Prel(beta,:,:)); Prelband(:,4) = tmp(:);
tmp = nanmean(Prel(alpha,:,:)); Prelband(:,3) = tmp(:);
tmp = nanmean(Prel(theta,:,:)); Prelband(:,2) = tmp(:);
tmp = nanmean(Prel(delta,:,:)); Prelband(:,1) = tmp(:);
clear tmp

%% Plot the result
subplot(3,1,1)
boxplot(Prelband,'labels', {'1-4Hz', '4-7Hz', '7-13Hz', '13-30Hz', ...
    '13-20Hz', '20-30Hz', '30-70Hz', '60-80Hz', '250-350Hz'})
ylim([0 5]), grid on
ylabel('P_{OFF}/P_{ON}')