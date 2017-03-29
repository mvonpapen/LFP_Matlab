%% Test PSD (total, coh, inc, vc) changes from OFF to ON


%% Load LFP data
[x, dat, CH, art, rigor, ndat, x, LFPactivity] = ...
    load_dat_from_DATA ( 'experiment', 'Faust', 'activity', 0, 'commontime', 1, 'noEMG', 1);

load DATA_last.mat

%% Parameters
Nf  = 60;
f   = logspace(0, log10(350), Nf);
sig = 0.495;

%% Determine sites and analysis
% % Sites with rest ON OFF and rigor>=30
% i0   = [1  6 11 18 26 32 36 41 47 65];
% iend = [5 10 17 25 31 68 40 46 55 67];
% ref  = [1  3  1  0  2  3  0  2  0  1]; %hand-picked reference electrodes (act=0) 
% % Sites with hold ON OFF and nLFP>1
% i0   = [1 3 5  9 11 21];
% iend = [2 4 6 10 23 22];
% ref  = [1 3 1  2  3  1]; %hand-picked reference electrodes (act=0) 
% Sites with fist ON OFF and rigor>=30
i0   = [1 3 6 8 10 11 14 18];
iend = [2 4 7 9 20 12 15 19];
ref  = [1 3 0 2  3  0  0  1]; %hand-picked reference electrodes (act=0)

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


%% STN power relative to nonSTN-reference-electrode
k=1;
for n=1:Nex
    STN(n,1:Nch(n))=LFPact{n};
    i1 = find(STN(n,:)>=1);
    i2 = ref(n);
    if i2 == 0
        continue
    end
    for i=1:length(i1)
        if i1(i) == i2
            continue
        end
        mpsd_tot_off(:,k) = Ptot_off(:,i1(i),n);
        mpsd_vc_off(:,k)  = Pvc_off(:,i1(i),i2,n);
        mpsd_tot_on(:,k)  = Ptot_on(:,i1(i),n);
        mpsd_vc_on(:,k)   = Pvc_on(:,i1(i),i2,n);
        k=k+1;
    end
end


mpsd_vc_off(mpsd_vc_off==mpsd_tot_off) = NaN;
mpsd_vc_on( mpsd_vc_on ==mpsd_tot_on ) = NaN;

for i=1:30;
    sig_tot(i)=signrank( mpsd_tot_off(i,:) ./ mpsd_tot_on(i,:), 1 );
    sig_novc(i)=signrank((mpsd_tot_off(i,:)-mpsd_vc_off(i,:)) ./ ...
        (mpsd_tot_on(i,:)-mpsd_vc_on(i,:)),1);
end


%% Plot results
figure
subplot(2,1,1)
loglog( f, mpsd_tot_off ./ mpsd_tot_on )
title('Total power change during rest'), ylim([0.05 20])
ylabel('P_{OFF}/P_{ON}'), xlabel('f [Hz]')
hold all
plot( f, nanmean( mpsd_tot_off ./ mpsd_tot_on, 2 ), '--k', 'linew', 2 )
i= sig_tot<0.05; if any(i); plot(f(i),10,'*k','markers',5); end
i= sig_tot<0.01; if any(i); plot(f(i),13,'+k','markers',5); end

subplot(2,1,2)
loglog(f,(mpsd_tot_off-mpsd_vc_off) ./ ...
    (mpsd_tot_on-mpsd_vc_on))
title('Change of power w/o vol. cond. during rest'), ylim([0.05 20])
ylabel('P_{OFF}/P_{ON}'), xlabel('f [Hz]')
hold all
plot(f, nanmean((mpsd_tot_off-mpsd_vc_off) ./ ...
    (mpsd_tot_on-mpsd_vc_on),2), '--k', 'linew', 2)
i= sig_novc<0.05; if any(i); plot(f(i),10,'*k','markers',5); end
i= sig_novc<0.01; if any(i); plot(f(i),13,'+k','markers',5); end