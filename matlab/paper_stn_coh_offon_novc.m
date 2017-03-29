%% Test PSD (total, coh, inc, vc) changes from OFF to ON


%% Load LFP data
[~, dat, CH, art, rigor, ndat, ~, LFPactivity] = ...
    load_dat_from_DATA_nb ( 'experiment', 'Ruhe', 'activity', 0, 'commontime', 1, 'noEMG', 1);

load /home/mvonpapen/Neuro/DATA/DATA_last.mat

%% Parameters
Nf  = 30;
f   = logspace(0, log10(90), Nf);
sig = 0.495;
phthres = 10;

%% Determine sites and analysis
% Sites with rest ON OFF and rigor>=30
i0   = [2  6 11 18 26 32 36 41 47 65];
iend = [5 10 17 25 31 68 40 46 55 67];
ref  = [1  3  1  0  2  3  0  2  0  1]; %hand-picked reference electrodes (act=0) 
% % Sites with hold ON OFF and nLFP>1
% i0   = [1 3 5  9 11 21];
% iend = [2 4 6 10 23 22];
% ref  = [1 3 1  2  3  1]; %hand-picked reference electrodes (act=0) 
% % Sites with fist ON OFF and rigor>=30
% i0   = [1 3 6 8 10 11 14 18];
% iend = [2 4 7 9 20 12 15 19];
% ref  = [1 3 0 2  3  0  0  1]; %hand-picked reference electrodes (act=0)

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
        [~,W,coi,Ptot_off(:,1:Nch(i),i)] = procdata(DataOFF{i}, 'freq', f, 'filter', [], 'art', ArtOFF{i});
        [C,Wxy]   = wave_acoh ( W, f );
        [c1d_off(:,1:Nch(i),1:Nch(i),i), c1d_novc_off(:,1:Nch(i),1:Nch(i),i)] ...
            = acoh1d( f, C, coi, Wxy, phthres );
                
        [~,W,coi,Ptot_on(:,1:Nch(i),i) ] = procdata(DataON{i},  'freq', f, 'filter', [], 'art', ArtON{i} );
        [C,Wxy]   = wave_acoh ( W, f );
        [c1d_on(:,1:Nch(i) ,1:Nch(i),i), c1d_novc_on(:,1:Nch(i) ,1:Nch(i),i)] ...
            = acoh1d( f, C , coi,  Wxy, phthres  );
    end
                  
end

clear DataOFF DataON ArtOFF ArtON x W Wxy C coi



%% Determine if coherences are within STN (2), STN<->non-STN (1), or
%% outside STN (0)
flag=NaN(max(Nch),max(Nch),Nex);
for k=1:Nex
    for j=1:Nch(k)
        for i=j+1:Nch(k)
            if  all([LFPact{k}(j) LFPact{k}(i)])
                flag(i,j,k)=2;
            else
                if xor(LFPact{k}(j), LFPact{k}(i))
                    flag(i,j,k)=1;
                else
                    flag(i,j,k)=0;
                end
            end
        end
    end
end


for i=1:30
    j = flag==0;
    coh0_off(i,:)      = c1d_off(i,j);
    coh0_on(i,:)       = c1d_on(i,j);
    coh0_novc_off(i,:) = c1d_novc_off(i,j);
    coh0_novc_on(i,:)  = c1d_novc_on(i,j);
    j = flag==1;
    coh1_off(i,:)      = c1d_off(i,j);
    coh1_on(i,:)       = c1d_on(i,j);
    coh1_novc_off(i,:) = c1d_novc_off(i,j);
    coh1_novc_on(i,:)  = c1d_novc_on(i,j);
    j = flag==2;
    coh2_off(i,:)      = c1d_off(i,j);
    coh2_on(i,:)       = c1d_on(i,j);
    coh2_novc_off(i,:) = c1d_novc_off(i,j);
    coh2_novc_on(i,:)  = c1d_novc_on(i,j);
    j = flag>=1;
    coh3_off(i,:)      = c1d_off(i,j);
    coh3_on(i,:)       = c1d_on(i,j);
    coh3_novc_off(i,:) = c1d_novc_off(i,j);
    coh3_novc_on(i,:)  = c1d_novc_on(i,j);
end
for i=1:30;
    sig1(i,1)= signrank( coh1_off(i,:)-coh1_on(i,:) ); 
    sig1(i,2)= signrank( coh1_novc_off(i,:)-coh1_novc_on(i,:) );
    sig2(i,1)= signrank( coh2_off(i,:)-coh2_on(i,:) ); 
    sig2(i,2)= signrank( coh2_novc_off(i,:)-coh2_novc_on(i,:) );
    sig3(i,1)= signrank( coh3_off(i,:)-coh3_on(i,:) );
    sig3(i,2)= signrank( coh3_novc_off(i,:)-coh3_novc_on(i,:) ); 
end

subplot(1,2,1), plot(f, nanmean(coh3_novc_off,2), '--k', 'linew', 3)
hold all, plot(f, nanmean(coh3_off,2), '-k', 'linew', 2)
subplot(1,2,2), plot(f, nanmean(coh3_novc_on,2), '--r', 'linew', 3)
hold all, plot(f, nanmean(coh3_on,2), '-r', 'linew', 2)