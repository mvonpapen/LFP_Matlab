%% Calculate all needed variables (coherence, PSD, phases) to check phase locking OFF/ON
clear all

%% Load LFP data
[t, dat, CH, art, rigor, ndat, pat, LFPactivity] = ...
    load_dat_from_DATA ( 'experiment', 'Ruhe', 'activity', 0, 'commontime', 1, 'noEMG', 1);

load /afs/.geo.uni-koeln.de/usr/neuro/DATA/DATA_last.mat

%% Parameters
Nf = 30;
f = logspace(0, log10(80), Nf);

%% Determine sites and analysis
% % Sites with rest ON OFF and rigor>=30
% i0   = [2  6 11 18 26 32 36 41 47 65];
% iend = [5 10 17 25 31 68 40 46 55 67];
% % Sites with hold ON OFF and nLFP>1
% i0   = [1 3 5  9 11 21];
% iend = [2 4 6 10 23 22];
% % Sites with fist ON OFF and rigor>=30
% i0   = [1 3 6 8 10 11 14 18];
% iend = [2 4 7 9 20 12 15 19];
i0=6; iend=10;
Nex = length(i0);



for i=1:Nex
    DataOFF{i}=dat{i0(i)};
    ArtOFF{i}=art{i0(i)};
    
    DataON{i}=dat{iend(i)};
    ArtON{i}=art{iend(i)};
    
    Channel{i}=CH{i0(i)};
    Patient{i}=DATA(ndat(iend(i))).patient;
    Nch(i)=size(dat{iend(i)},2);
    LFPact{i}=LFPactivity{iend(i)};
end

clear i i0 iend DATA dat art CH LFPactivity t ndat pat

   

%% Calculate wavelet transformation
Poff            = NaN(Nf,max(Nch),Nex);
Pon             = NaN(Nf,max(Nch),Nex);
c1d_off         = NaN(Nf,max(Nch),max(Nch),Nex);
c1d_on          = NaN(Nf,max(Nch),max(Nch),Nex);
c1d_novc_off    = NaN(Nf,max(Nch),max(Nch),Nex);
c1d_novc_on     = NaN(Nf,max(Nch),max(Nch),Nex);

for i=1:Nex
    
    [x,Woff,coi_off,Poff(:,1:Nch(i),i)] = procdata(DataOFF{i}, 'freq', f, ...
        'filter', [], 'art', ArtOFF{i});
    [x,Won, coi_on, Pon(:,1:Nch(i),i)]  = procdata(DataON{i},  'freq', f, ...
        'filter', [], 'art', ArtON{i} );
    
    if Nch(i)>1
        % Calculate coherent spectra of all channels w.r.t. non-STN
        [C_off,Wxy_off] = wave_acoh ( Woff, f );
        [C_on, Wxy_on]  = wave_acoh ( Won, f  );
        [c1d_off(:,1:Nch(i),1:Nch(i),i), c1d_novc_off(:,1:Nch(i),1:Nch(i),i)] ...
            = acoh1d( f, C_off, coi_off, Wxy_off );
        [c1d_on(:,1:Nch(i) ,1:Nch(i),i), c1d_novc_on(:,1:Nch(i) ,1:Nch(i),i)] ...
            = acoh1d( f, C_on , coi_on,  Wxy_on  );
    end
                  
end


% %% Put coherences together
% C1D_OFF     =[]; 
% C1D_ON      =[];
% C1D_NOVC_OFF=[]; 
% C1D_NOVC_ON =[];
% for i=1:4
%     for j=i+1:4
%         for k=1:10
%             C1D_OFF     = [C1D_OFF      c1d_off(:,i,j,k)     ];
%             C1D_ON      = [C1D_ON       c1d_on(:,i,j,k)      ];
%             C1D_NOVC_OFF= [C1D_NOVC_OFF c1d_novc_off(:,i,j,k)];
%             C1D_NOVC_ON = [C1D_NOVC_ON  c1d_novc_on(:,i,j,k) ];
%         end
%     end
% end
% C1D_OFF(C1D_OFF==0)             =NaN;
% C1D_ON(C1D_ON==0)               =NaN;
% C1D_NOVC_OFF(C1D_NOVC_OFF==0)   =NaN;
% C1D_NOVC_ON(C1D_NOVC_ON==0)     =NaN;
% 
% C1D_OFF     =C1D_OFF(:,~isnan(C1D_OFF(20,:)));
% C1D_ON      =C1D_ON(:,~isnan(C1D_ON(20,:)));
% C1D_NOVC_OFF=C1D_NOVC_OFF(:,~isnan(C1D_NOVC_OFF(20,:)));
% C1D_NOVC_ON =C1D_NOVC_ON(:,~isnan(C1D_NOVC_ON(20,:)));
% 
% [ind_off, Cmean_off]            = ssca(C1D_OFF,     3, [14:25], 0);
% [ind_on, Cmean_on]              = ssca(C1D_ON,     3, [14:25], 0);
% [ind_novc_off, Cmean_novc_off]  = ssca(C1D_NOVC_OFF,3, [14:25], 0);
% [ind_novc_on, Cmean_novc_on]    = ssca(C1D_NOVC_ON, 3, [14:25], 0);
% 
% Cdiff     =[nanmean(C1D_OFF     ./C1D_ON,     2) nanstd(C1D_OFF'     ./C1D_ON')'];
% Cdiff_novc=[nanmean(C1D_NOVC_OFF./C1D_NOVC_ON,2) nanstd(C1D_NOVC_OFF'./C1D_NOVC_ON')'];