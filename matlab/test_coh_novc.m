%% Test coherences with and without vc-correction
clear all

%% Load LFP data
[t, dat, CH, art, rigor, ndat, pat, LFPactivity] = ...
    load_dat_from_DATA ( 'experiment', 'Ruhe', 'activity', 0, 'commontime', 1, 'noEMG', 1);

load /afs/.geo.uni-koeln.de/usr/neuro/DATA/DATA_last.mat

%% Parameters
Nf = 30;
f = logspace(0, log10(90), Nf);

%% Determine sites and analysis
% Sites with rest ON OFF and rigor>=30
i0   = [2  6 11 18 26 32 36 41 47 65];
iend = [5 10 17 25 31 68 40 46 55 67];
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

clear i i0 iend DATA dat art CH LFPactivity t ndat pat Patient2

   

%% Calculate wavelet transformation
c1d_off         = NaN(Nf,max(Nch),max(Nch),Nex);
c1d_on          = NaN(Nf,max(Nch),max(Nch),Nex);
c1d_novc_off    = NaN(Nf,max(Nch),max(Nch),Nex);
c1d_novc_on     = NaN(Nf,max(Nch),max(Nch),Nex);

for i=1:Nex
    
    [x,Woff,coi_off] = procdata(DataOFF{i}, 'freq', f, 'filter', [], 'art', ArtOFF{i});
    [x,Won, coi_on ] = procdata(DataON{i},  'freq', f, 'filter', [], 'art', ArtON{i} );
    
    if Nch(i)>1
        % Calculate coherent spectra of all channels w.r.t. non-STN
        [C_off,Wxy_off]  = wave_acoh ( Woff, f );
        [C_on, Wxy_on ]  = wave_acoh ( Won,  f );
        [c1d_off(:,1:Nch(i),1:Nch(i),i), c1d_novc_off(:,1:Nch(i),1:Nch(i),i)] ...
            = acoh1d( f, C_off, coi_off, Wxy_off );
        [c1d_on(:,1:Nch(i) ,1:Nch(i),i), c1d_novc_on(:,1:Nch(i) ,1:Nch(i),i)] ...
            = acoh1d( f, C_on , coi_on,  Wxy_on  );
    end
                  
end

%% Determine if coherences are within STN (2), STN<->non-STN (1), or
%% outside STN (0)
flag=NaN(max(Nch),max(Nch),Nex);
for k=1:Nex
    for j=1:Nch(k)
        for i=j+1:Nch(k) %setdiff(1:Nch(k),j)
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


%% Plot results in three figures according to flag
col = jet(20);

%% flag=0
[i,j]=find(flag==0);
k=ceil(j/max(Nch));
j=mod(j,max(Nch));
j(j==0) = max(Nch);

figure
for n=1:length(i)
    subplot(2,1,1), title('Non-STN <-> non-STN coherence (no correction)')
    semilogx(f,c1d_off(:,i(n),j(n),k(n)) ,'-', 'color', col(n,:))
    hold all, ylim([0 1])    
    semilogx(f,c1d_on(:,i(n),j(n),k(n)) ,'--','color', col(n,:))
    
    subplot(2,1,2), title('Non-STN <-> non-STN coherence (w/o volume conduction)')
    semilogx(f,c1d_novc_off(:,i(n),j(n),k(n)),'-', 'color', col(n,:))
    hold all, ylim([0 1])
    semilogx(f,c1d_novc_on(:,i(n),j(n),k(n)),'--', 'color', col(n,:))
    
    mc1d_off_0(:,n) = c1d_novc_off(:,i(n),j(n),k(n));
    mc1d_on_0(:,n)  = c1d_novc_on(:,i(n),j(n),k(n)) ;
end

%% flag=1
[i,j]=find(flag==1);
k=ceil(j/max(Nch));
j=mod(j,max(Nch));
j(j==0) = max(Nch);

figure
for n=1:length(i)
    subplot(2,1,1), title('STN <-> non-STN coherence (no correction)')
    semilogx(f,c1d_off(:,i(n),j(n),k(n)) ,'-', 'color', col(n,:))
    hold all, ylim([0 1])    
    semilogx(f,c1d_on(:,i(n),j(n),k(n)) ,'--','color', col(n,:))
    
    subplot(2,1,2), title('STN <-> non-STN coherence (w/o volume conduction)')
    semilogx(f,c1d_novc_off(:,i(n),j(n),k(n)),'-', 'color', col(n,:))
    hold all, ylim([0 1])
    semilogx(f,c1d_novc_on(:,i(n),j(n),k(n)),'--', 'color', col(n,:))
    
    mc1d_off_1(:,n) = c1d_novc_off(:,i(n),j(n),k(n));
    mc1d_on_1(:,n)  = c1d_novc_on(:,i(n),j(n),k(n)) ;
end

%% flag=2
[i,j]=find(flag==2);
k=ceil(j/max(Nch));
j=mod(j,max(Nch));
j(j==0) = max(Nch);

figure
for n=1:length(i)
    subplot(2,1,1), title('STN <-> STN coherence (no correction)')
    semilogx(f,c1d_off(:,i(n),j(n),k(n)) ,'-', 'color', col(n,:))
    hold all, ylim([0 1])    
    semilogx(f,c1d_on(:,i(n),j(n),k(n)) ,'--','color', col(n,:))
    
    subplot(2,1,2), title('STN <-> STN coherence (w/o volume conduction)')
    semilogx(f,c1d_novc_off(:,i(n),j(n),k(n)),'-', 'color', col(n,:))
    hold all, ylim([0 1])
    semilogx(f,c1d_novc_on(:,i(n),j(n),k(n)),'--', 'color', col(n,:))
    
    mc1d_off_2(:,n) = c1d_novc_off(:,i(n),j(n),k(n));
    mc1d_on_2(:,n)  = c1d_novc_on(:,i(n),j(n),k(n)) ;
end