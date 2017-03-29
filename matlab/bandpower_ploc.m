%% Check changes in power OFF/ON
clear all

%% Load LFP data
[t, dat, CH, art, rigor, ndat, pat, LFPactivity] = ...
    load_dat_from_DATA ( 'experiment', 'Faust', 'activity', 0, 'commontime', 1, 'noEMG', 1);
for i=1:length(dat)
    nLFP(i)=size(dat{i},2);
end
load /afs/.geo.uni-koeln.de/usr/neuro/DATA/DATA_last.mat

%% Parameters
Nf = 30;
f = logspace(0, log10(70), Nf);
delta=find(f>=1 & f<4);
theta=find(f>=4 & f<7);
alpha=find(f>=7 & f<13);
beta=find(f>=13 & f<30);
beta1=find(f>=13 & f<20);
beta2=find(f>=20 & f<30);
gamma=find(f>=30 & f<70);

%% Determine sites and analysis
LFPtypes={'raw'; 'Ploc'};
% % Sites with rest ON OFF and nLFP>1
% i0   = [1  6 11 26 32 41 65];
% iend = [5 10 17 31 68 46 67];
% ref = [2 3 3 3 3 1 1];
% % Sites with hold ON OFF and nLFP>1
% i0   = [1 3 5  9 11 21];
% iend = [2 4 6 10 23 22];
% ref = [2 3 3 3 3 1];
% Sites with fist ON OFF and nLFP>1
i0   = [1 3 8 10 18];
iend = [2 4 9 20 19];
ref = [2 3 3 3 1];
Nex = length(i0);

for i=1:Nex
    DataOFF{i}=dat{i0(i)};
    RigorOFF(i)=rigor(i0(i));
    ChannelOFF{i}=CH{i0(i)};
    NdatOFF(i)=ndat(i0(i));
    ArtOFF{i}=art{i0(i)};
    Patient{i}=DATA(ndat(i0(i))).patient;
    DataON{i}=dat{iend(i)};
    RigorON(i)=rigor(iend(i));
    ChannelON{i}=CH{iend(i)};
    NdatON(i)=ndat(iend(i));
    ArtON{i}=art{iend(i)};
    Patient2{i}=DATA(ndat(iend(i))).patient;
    NLFP(i)=size(dat{iend(i)},2);
    LFPact{i}=LFPactivity{iend(i)};
end

clear i i0 iend



%% Calculate phases relative to non-STN channel
Nch = NLFP;
for i=1:Nex
    D{i,1} = DataOFF{i};
    D{i,2} = DataON {i};
    A{i,1} = ArtOFF{i};
    A{i,2} = ArtON {i};
end
   

%% Calculate wavelet transformation
Poff = NaN(Nf,max(Nch),Nex);
Pon  = NaN(Nf,max(Nch),Nex);

Pcoh_off = NaN(Nf,max(Nch),Nex);
Pcoh_on  = NaN(Nf,max(Nch),Nex);
Pinc_off = NaN(Nf,max(Nch),Nex);
Pinc_on  = NaN(Nf,max(Nch),Nex);
Pvc_off  = NaN(Nf,max(Nch),Nex);
Pvc_on   = NaN(Nf,max(Nch),Nex);

for i=1:Nex
    [x,Woff,coi_off,Poff(:,1:Nch(i),i)] = procdata(D{i,1}, 'freq', f, ...
        'filter', [45 55; 90 110], 'art', A{i,1});
    [x,Won, coi_on, Pon(:,1:Nch(i),i)]  = procdata(D{i,2}, 'freq', f, ...
        'filter', [45 55; 90 110], 'art', A{i,2});
    clear x
    
    % Calculate coherent spectra of all channels w.r.t. non-STN
    [C_off,Wxy_off] = wave_coh ( Woff, Woff(:,:,ref(i)), f );
    [C_on ,Wxy_on ] = wave_coh ( Won , Won (:,:,ref(i)), f );
    [Pcoh_off(:,1:Nch(i),i), Pinc_off(:,1:Nch(i),i), ...
        Pvc_off(:,1:Nch(i),i)] = psd_coh( f, Woff, C_off, ...
        coi_off, 0.5, Wxy_off);
    [Pcoh_on(:,1:Nch(i),i) , Pinc_on(:,1:Nch(i),i) , ...
        Pvc_on(:,1:Nch(i),i) ] = psd_coh( f, Won , C_on , ...
        coi_on , 0.5, Wxy_on );
                  
end


%% Determine ratio

for typ=1:2


    LFPtype=LFPtypes{typ};


    Prel=NaN(Nf,max(Nch),Nex);
    for i=1:Nex

        switch LFPtype

            case 'Ploc'
                x = setxor (1:Nch(i), ref(i));
                Prel(:,1:Nch(i)-1,i) = ...
                    ( Poff(:,x,i)-Pvc_off(:,x,i) ) ./ ( Pon(:,x,i)-Pvc_on(:,x,i) );

            case 'raw'
                x = setxor (1:Nch(i), ref(i));
                Prel(:,1:Nch(i)-1,i) = Poff(:,x,i) ./ Pon(:,x,i) ;
        end

    end
    Prel(isinf(Prel))=9e9;
    Prel(Prel==0)=9e-9;
    clear i ind

    %% Determine averages in freq. bands
    clear Prelband
    tmp = nanmean(Prel(gamma,:,:)); Prelband(:,7) = tmp(:);
    tmp = nanmean(Prel(beta2,:,:)); Prelband(:,6) = tmp(:);
    tmp = nanmean(Prel(beta1,:,:)); Prelband(:,5) = tmp(:);
    tmp = nanmean(Prel(beta,:,:));  Prelband(:,4) = tmp(:);
    tmp = nanmean(Prel(alpha,:,:)); Prelband(:,3) = tmp(:);
    tmp = nanmean(Prel(theta,:,:)); Prelband(:,2) = tmp(:);
    tmp = nanmean(Prel(delta,:,:)); Prelband(:,1) = tmp(:);
    clear tmp

    %% Plot the result
    subplot(2,1,typ)
    % figure
    boxplot(log10(Prelband),'labels', {'\delta', '\theta', '\alpha', '\beta', ...
        '\beta1', '\beta2', '\gamma'})
    h = findobj(gca, 'type', 'text');
    set(h, 'Interpreter', 'tex');
    ylim([-1 1]), grid on
    ylabel('log_{10}(P_{OFF}/P_{ON})')
    hold all
    for i=1:7
        p(i)=signrank(log(Prelband(:,i)));
        if p(i)<0.01;        
            plot([i-0.05 i+0.05],[0.9 0.9],'*k');
        elseif p(i)<0.05;
            plot(i,0.9,'*k');
        end
    end
    prob(:,typ)=p;
end

save 'bandpower_fist.mat' Pon Poff Pvc_off Pvc_on Nch Nex ref f alpha beta beta1 beta2 delta theta gamma prob Pcoh_on Pcoh_off Pinc_on Pinc_off