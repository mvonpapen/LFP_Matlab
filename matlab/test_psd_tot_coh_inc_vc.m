%% Test PSD (total, coh, inc, vc) changes from OFF to ON
clear all

%% Load LFP data
[t, dat, CH, art, rigor, ndat, pat, LFPactivity] = ...
    load_dat_from_DATA ( 'experiment', 'Ruhe', 'activity', 0, 'commontime', 1, 'noEMG', 1);

load /afs/.geo.uni-koeln.de/usr/neuro/DATA/DATA_last.mat

%% Parameters
Nf  = 60;
f   = logspace(0, log10(350), Nf);
sig = 0.495;

%% Determine sites and analysis
% Sites with rest ON OFF and rigor>=30
i0   = [1  6 11 18 26 32 36 41 47 65];
iend = [5 10 17 25 31 68 40 46 55 67];
ref  = [1  3  1  3  0  3  0  2  0  1]; %hand-picked reference electrodes (act=0) 
% % Sites with hold ON OFF and nLFP>1
% i0   = [1 3 5  9 11 21];
% iend = [2 4 6 10 23 22];
% % Sites with fist ON OFF and rigor>=30
% i0   = [1 3 6 8 10 11 14 18];
% iend = [2 4 7 9 20 12 15 19];
% i0=6; iend=10;

Nex = length(i0);

for i=1:Nex
    DataOFF{i}=dat{i0(i)};
    ArtOFF{i}=art{i0(i)};
    
    DataON{i}=dat{iend(i)};
    ArtON{i}=art{iend(i)};
    
    Channel{i}=CH{i0(i)};
    Patient{i}=DATA(ndat(i0(i))).patient;
    Patient2{i}=DATA(ndat(iend(i))).patient;
    Nch(i)=size(dat{iend(i)},2);
    LFPact{i}=LFPactivity{iend(i)};
    
    if ~strcmp(Patient{i},Patient2{i})
        error('ERROR: i0 and iend refer to different patients for i=%i!', i)
    end
end
rigor=rigor(iend);
LFPact{5}(2) = NaN; %PSD sieht hier komisch aus, manual out

clear i i0 iend DATA dat art CH LFPactivity t ndat pat Patient2

   

%% Calculate spectral powers: total, coherent, incoherent and volume
%% conducted
Ptot_off    = NaN(Nf,max(Nch),Nex);
Pcoh_off    = NaN(Nf,max(Nch),max(Nch),Nex);
Pinc_off    = NaN(Nf,max(Nch),max(Nch),Nex);
Pvc_off     = NaN(Nf,max(Nch),max(Nch),Nex);
Ptot_on     = NaN(Nf,max(Nch),Nex);
Pcoh_on     = NaN(Nf,max(Nch),max(Nch),Nex);
Pinc_on     = NaN(Nf,max(Nch),max(Nch),Nex);
Pvc_on      = NaN(Nf,max(Nch),max(Nch),Nex);

for i=1:Nex
    if Nch(i)>1
        [x,W,coi,Ptot_off(:,1:Nch(i),i)] = procdata(DataOFF{i}, 'freq', f, 'filter', [], 'art', ArtOFF{i});
        [C,Wxy]   = wave_acoh ( W, f );
        [Pcoh_off(:,1:Nch(i),1:Nch(i),i), Pinc_off(:,1:Nch(i),1:Nch(i),i), ...
            Pvc_off(:,1:Nch(i),1:Nch(i),i)] = psd_coh ( f, W, C, coi, sig, Wxy );
                
        [x,W,coi,Ptot_on(:,1:Nch(i),i) ] = procdata(DataON{i},  'freq', f, 'filter', [], 'art', ArtON{i} );
        [C,Wxy]   = wave_acoh ( W, f );
        [Pcoh_on(:,1:Nch(i),1:Nch(i),i),  Pinc_on(:,1:Nch(i),1:Nch(i),i),  ...
            Pvc_on(:,1:Nch(i),1:Nch(i),i) ] = psd_coh ( f, W, C, coi, sig, Wxy );
    end
                  
end

clear DataOFF DataON ArtOFF ArtON x W Wxy C coi





%% Determine if electrodes are within STN (2), STN<->non-STN (1), or
%% outside STN (0) and plot results
STN = NaN(Nex,max(Nch));

%% nonSTN power relative to nonSTN-electrodes
k=1;
for n=1:Nex
    STN(n,1:Nch(n))=LFPact{n};
    i1 = find(STN(n,:)==0);
    i2 = find(STN(n,:)==0);
    for i=1:length(i1)
        for j=1:length(i2)
            if i1(i)==i2(j)
                continue
            end
            mpsd_coh_on_0(:,k)  = Pcoh_on(:,i1(i),i2(j),n);
            mpsd_inc_on_0(:,k)  = Pinc_on(:,i1(i),i2(j),n);
            mpsd_vc_on_0(:,k)   = Pvc_on(:,i1(i),i2(j),n);
            mpsd_coh_off_0(:,k) = Pcoh_off(:,i1(i),i2(j),n);
            mpsd_inc_off_0(:,k) = Pinc_off(:,i1(i),i2(j),n);
            mpsd_vc_off_0(:,k)  = Pvc_off(:,i1(i),i2(j),n);
            k=k+1;
        end
    end
end

%% STN power relative to reference electrodes
k=1;
for n=1:Nex
    i1 = find(STN(n,:)>=1);
%     i2 = find(STN(n,:)==0);
    i2 = ref(n);
    if i2 == 0
        continue
    end
    for i=1:length(i1)
        for j=1:length(i2)
            mpsd_coh_on_1(:,k)  = Pcoh_on(:,i1(i),i2(j),n);
            mpsd_inc_on_1(:,k)  = Pinc_on(:,i1(i),i2(j),n);
            mpsd_vc_on_1(:,k)   = Pvc_on(:,i1(i),i2(j),n);
            mpsd_coh_off_1(:,k) = Pcoh_off(:,i1(i),i2(j),n);
            mpsd_inc_off_1(:,k) = Pinc_off(:,i1(i),i2(j),n);
            mpsd_vc_off_1(:,k)  = Pvc_off(:,i1(i),i2(j),n);
            k=k+1;
        end
    end
end

%% STN power relative to (non)STN-electrodes
k=1;
for n=1:Nex
    i1 = find(STN(n,:)>=1);
    i2 = find(STN(n,:)>=1);
    for i=1:length(i1)
        for j=1:length(i2)
            if i1(i)==i2(j)
                continue
            end
            mpsd_coh_on_2(:,k)  = Pcoh_on(:,i1(i),i2(j),n);
            mpsd_inc_on_2(:,k)  = Pinc_on(:,i1(i),i2(j),n);
            mpsd_vc_on_2(:,k)   = Pvc_on(:,i1(i),i2(j),n);
            mpsd_coh_off_2(:,k) = Pcoh_off(:,i1(i),i2(j),n);
            mpsd_inc_off_2(:,k) = Pinc_off(:,i1(i),i2(j),n);
            mpsd_vc_off_2(:,k)  = Pvc_off(:,i1(i),i2(j),n);
            k=k+1;
        end
    end
end


for i=1:30;
    j = find( STN' == 0);
    sig0(1,i)=signrank( Ptot_off(i,j) ./ Ptot_on(i,j), 1 );
    sig0(2,i)=signrank((mpsd_coh_off_0(i,:)+mpsd_inc_off_0(i,:)-mpsd_vc_off_0(i,:)) ./ ...
        (mpsd_coh_on_0(i,:)+mpsd_inc_on_0(i,:)-mpsd_vc_on_0(i,:)),1);
    
    j = find( STN' > 0);
    sig1(1,i)=signrank( Ptot_off(i,j) ./ Ptot_on(i,j), 1 );
    sig1(2,i)=signrank((mpsd_coh_off_1(i,:)+mpsd_inc_off_1(i,:)-mpsd_vc_off_1(i,:)) ./ ...
        (mpsd_coh_on_1(i,:)+mpsd_inc_on_1(i,:)-mpsd_vc_on_1(i,:)),1);
    
    sig2(1,i)= sig2(1,i); %signrank((mpsd_coh_off_2(i,:)+mpsd_inc_off_2(i,:)) ./ ...
%         (mpsd_coh_on_2(i,:)+mpsd_inc_on_2(i,:)),1);
    sig2(2,i)=signrank((mpsd_coh_off_2(i,:)+mpsd_inc_off_2(i,:)-mpsd_vc_off_2(i,:)) ./ ...
        (mpsd_coh_on_2(i,:)+mpsd_inc_on_2(i,:)-mpsd_vc_on_2(i,:)),1);
end

figure
subplot(2,3,1)
j = find( STN' == 0);
loglog(f,Ptot_off(:,j) ./ Ptot_on(:,j))
title('total OFF/ON nonSTN<-nonSTN'), ylim([0.05 20])
hold all
plot(f, nanmean(Ptot_off(:,j) ./ Ptot_on(:,j),2), '--', 'linew', 2)
i= sig0(1,:)<0.05; if any(i); plot(f(i),10,'*k','linew',2); end
i= sig0(1,:)<0.01; if any(i); plot(f(i),13,'+k','linew',2); end
subplot(2,3,4)
loglog(f,(mpsd_coh_off_0+mpsd_inc_off_0-mpsd_vc_off_0) ./ ...
    (mpsd_coh_on_0+mpsd_inc_on_0-mpsd_vc_on_0))
title('total-vc OFF/ON nonSTN<-nonSTN'), ylim([0.05 20])
hold all
plot(f, nanmean((mpsd_coh_off_0+mpsd_inc_off_0-mpsd_vc_off_0) ./ ...
    (mpsd_coh_on_0+mpsd_inc_on_0-mpsd_vc_on_0),2), '--', 'linew', 2)
i= sig0(2,:)<0.05; if any(i); plot(f(i),10,'*k','linew',2); end
i= sig0(2,:)<0.01; if any(i); plot(f(i),13,'+k','linew',2); end

subplot(2,3,2)
loglog(f,(mpsd_coh_off_1+mpsd_inc_off_1)./(mpsd_coh_on_1+mpsd_inc_on_1))
title('total OFF/ON STN<-nonSTN'), ylim([0.05 20])
hold all
plot(f, nanmean((mpsd_coh_off_1+mpsd_inc_off_1)./(mpsd_coh_on_1+mpsd_inc_on_1), 2), '--', 'linew', 2)
i= sig1(1,:)<0.05; if any(i); plot(f(i),10,'*k','linew',2); end
i= sig1(1,:)<0.01; if any(i); plot(f(i),13,'+k','linew',2); end
subplot(2,3,5)
loglog(f,(mpsd_coh_off_1+mpsd_inc_off_1-mpsd_vc_off_1) ./ ...
    (mpsd_coh_on_1+mpsd_inc_on_1-mpsd_vc_on_1))
title('total-vc OFF/ON STN<-nonSTN'), ylim([0.05 20])
hold all
plot(f, nanmean((mpsd_coh_off_1+mpsd_inc_off_1-mpsd_vc_off_1) ./ ...
    (mpsd_coh_on_1+mpsd_inc_on_1-mpsd_vc_on_1),2), '--', 'linew', 2)
i= sig1(2,:)<0.05; if any(i); plot(f(i),10,'*k','linew',2); end
i= sig1(2,:)<0.01; if any(i); plot(f(i),13,'+k','linew',2); end

subplot(2,3,3)
loglog(f,(mpsd_coh_off_2+mpsd_inc_off_2)./(mpsd_coh_on_2+mpsd_inc_on_2))
title('total OFF/ON STN<-STN'), ylim([0.05 20])
hold all
plot(f, nanmean((mpsd_coh_off_2+mpsd_inc_off_2)./(mpsd_coh_on_2+mpsd_inc_on_2),2), '--', 'linew', 2)
i= sig2(1,:)<0.05; if any(i); plot(f(i),10,'*k','linew',2); end
i= sig2(1,:)<0.01; if any(i); plot(f(i),13,'+k','linew',2); end
subplot(2,3,6)
loglog(f,(mpsd_coh_off_2+mpsd_inc_off_2-mpsd_vc_off_2) ./ ...
    (mpsd_coh_on_2+mpsd_inc_on_2-mpsd_vc_on_2))
title('total-vc OFF/ON STN<-STN'), ylim([0.05 20])
hold all
plot(f, nanmean((mpsd_coh_off_2+mpsd_inc_off_2-mpsd_vc_off_2) ./ ...
    (mpsd_coh_on_2+mpsd_inc_on_2-mpsd_vc_on_2),2), '--', 'linew', 2)
i= sig2(2,:)<0.05; if any(i); plot(f(i),10,'*k','linew',2); end
i= sig2(2,:)<0.01; if any(i); plot(f(i),13,'+k','linew',2); end