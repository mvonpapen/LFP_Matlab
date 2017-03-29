%% Check changes in power OFF/ON
clear all

%% Load LFP data
[t, dat, CH, art, rigor, ndat, pat, LFPactivity] = ...
    load_dat_from_DATA ( 'experiment', 'Ruhe', 'activity', 0, 'commontime', 1, 'noEMG', 1);
for i=1:length(dat)
    nLFP(i)=size(dat{i},2);
end
load /afs/.geo.uni-koeln.de/usr/neuro/DATA/DATA_last.mat

for typ=2:3

%% Parameters
Nf = 50;
f = logspace(0, 2.6, Nf);
delta=find(f>=1 & f<4);
theta=find(f>=4 & f<7);
alpha=find(f>=7 & f<13);
beta=find(f>=13 & f<30);
beta1=find(f>=13 & f<20);
beta2=find(f>=20 & f<30);
gamma=find(f>=30 & f<70);
gamma2=find(f>=60 & f<80);
HF=find(f>=250 & f<350);

%% Determine sites and analysis
LFPtype={'raw'; 'Ploc'; 'dLFP'};
LFPtype=LFPtype{typ};
% Sites with rest ON OFF and nLFP>1
i0   = [2  6 11 26 32 41 65];
iend = [5 10 17 31 68 46 67];
ref = [2 3 3 3 3 1 1];
% % Sites with hold ON OFF and nLFP>1
% i0   = [1 3 5  9 11 21];
% iend = [2 4 6 10 23 22];
% ref = [2 3 3 3 3 1];
% % Sites with fist ON OFF and nLFP>1
% i0   = [1 3 8 10 18];
% iend = [2 4 9 20 19];
% ref = [2 3 3 3 1];
Nex = length(i0);



for i=1:Nex
    DataOFF{i}=dat{i0(i)};
    RigorOFF(i)=rigor(i0(i));
    ChannelOFF{i}=CH{i0(i)};
    NdatOFF(i)=ndat(i0(i));
    ArtOFF{i}=art{i0(i)};
%     Patient{i}=DATA(ndat(i0(i))).patient,
    DataON{i}=dat{iend(i)};
    RigorON(i)=rigor(iend(i));
    ChannelON{i}=CH{iend(i)};
    NdatON(i)=ndat(iend(i));
    ArtON{i}=art{iend(i)};
    Patient{i}=DATA(ndat(iend(i))).patient;
    NLFP(i)=size(dat{iend(i)},2);
    LFPact{i}=LFPactivity{iend(i)};
end

clear i i0 iend



switch LFPtype
    case 'dLFP'
        %% Calculate dLFP=x_1-x_2
        i=0;

        % BEMI
        i=i+1;
        D{i,1}(:,1)=(DataOFF{i}(:,1)-DataOFF{i}(:,2)); %CM*
        D{i,1}(:,2)=(DataOFF{i}(:,1)-DataOFF{i}(:,3)); %CA
        D{i,1}(:,3)=(DataOFF{i}(:,2)-DataOFF{i}(:,3)); %MA*
        D{i,2}(:,1)=(DataON{i}(:,1)-DataON{i}(:,2)); %CM*
        D{i,2}(:,2)=(DataON{i}(:,1)-DataON{i}(:,3)); %CA
        D{i,2}(:,3)=(DataON{i}(:,2)-DataON{i}(:,3)); %MA*
        flag([1 3],i)=1; %flag denotes non-STN data
        Nch(i)=3;

        % BIMA
        i=i+1;
        D{i,1}(:,1)=(DataOFF{i}(:,1)-DataOFF{i}(:,2)); %CA
        D{i,1}(:,2)=(DataOFF{i}(:,1)-DataOFF{i}(:,3)); %CP*
        D{i,1}(:,3)=(DataOFF{i}(:,1)-DataOFF{i}(:,4)); %CL
        D{i,1}(:,4)=(DataOFF{i}(:,2)-DataOFF{i}(:,3)); %AP*
        D{i,1}(:,5)=(DataOFF{i}(:,2)-DataOFF{i}(:,4)); %AL
        D{i,1}(:,6)=(DataOFF{i}(:,3)-DataOFF{i}(:,4)); %PL*
        D{i,2}(:,1)=(DataON{i}(:,1)-DataON{i}(:,2)); %CA
        D{i,2}(:,2)=(DataON{i}(:,1)-DataON{i}(:,3)); %CP*
        D{i,2}(:,3)=(DataON{i}(:,1)-DataON{i}(:,4)); %CL
        D{i,2}(:,4)=(DataON{i}(:,2)-DataON{i}(:,3)); %AP*
        D{i,2}(:,5)=(DataON{i}(:,2)-DataON{i}(:,4)); %AL
        D{i,2}(:,6)=(DataON{i}(:,3)-DataON{i}(:,4)); %PL*
        flag([2 4 6],i)=1;
        Nch(i)=6;

        % GRFR, no fist available
        i=i+1;
        D{i,1}(:,1)=(DataOFF{i}(:,1)-DataOFF{i}(:,2)); %CA
        D{i,1}(:,2)=(DataOFF{i}(:,1)-DataOFF{i}(:,3)); %CM
        D{i,1}(:,3)=(DataOFF{i}(:,2)-DataOFF{i}(:,3)); %AM
        D{i,2}(:,1)=(DataON{i}(:,1)-DataON{i}(:,2)); %CA
        D{i,2}(:,2)=(DataON{i}(:,1)-DataON{i}(:,3)); %CM
        D{i,2}(:,3)=(DataON{i}(:,2)-DataON{i}(:,3)); %AM
        Nch(i)=3;

        % REGE
        i=i+1;
        D{i,1}(:,1)=(DataOFF{i}(:,1)-DataOFF{i}(:,2)); %CL
        D{i,1}(:,2)=(DataOFF{i}(:,1)-DataOFF{i}(:,3)); %CA*
        D{i,1}(:,3)=(DataOFF{i}(:,1)-DataOFF{i}(:,4)); %CP
        D{i,1}(:,4)=(DataOFF{i}(:,2)-DataOFF{i}(:,3)); %LA*
        D{i,1}(:,5)=(DataOFF{i}(:,2)-DataOFF{i}(:,4)); %LP
        D{i,1}(:,6)=(DataOFF{i}(:,3)-DataOFF{i}(:,4)); %AP*
        D{i,2}(:,1)=(DataON{i}(:,1)-DataON{i}(:,2)); %CL
        D{i,2}(:,2)=(DataON{i}(:,1)-DataON{i}(:,3)); %CA*
        D{i,2}(:,3)=(DataON{i}(:,1)-DataON{i}(:,4)); %CP
        D{i,2}(:,4)=(DataON{i}(:,2)-DataON{i}(:,3)); %LA*
        D{i,2}(:,5)=(DataON{i}(:,2)-DataON{i}(:,4)); %LP
        D{i,2}(:,6)=(DataON{i}(:,3)-DataON{i}(:,4)); %AP*
        flag([2 4 6],i)=1;
        Nch(i)=6;

        % RIRO
        i=i+1;
        D{i,1}(:,1)=(DataOFF{i}(:,1)-DataOFF{i}(:,2)); %CP
        D{i,1}(:,2)=(DataOFF{i}(:,1)-DataOFF{i}(:,3)); %CM*
        D{i,1}(:,3)=(DataOFF{i}(:,2)-DataOFF{i}(:,3)); %PM*
        D{i,2}(:,1)=(DataON{i}(:,1)-DataON{i}(:,2)); %CP
        D{i,2}(:,2)=(DataON{i}(:,1)-DataON{i}(:,3)); %CM*
        D{i,2}(:,3)=(DataON{i}(:,2)-DataON{i}(:,3)); %PM*
        flag([2 3],i)=1;
        Nch(i)=3;

        % SUGI, no fist/hold available
        i=i+1;
        D{i,1}(:,1)=(DataOFF{i}(:,1)-DataOFF{i}(:,2)); %CA*
        D{i,1}(:,2)=(DataOFF{i}(:,1)-DataOFF{i}(:,3)); %CM*
        D{i,1}(:,3)=(DataOFF{i}(:,1)-DataOFF{i}(:,4)); %CP*
        D{i,1}(:,4)=(DataOFF{i}(:,2)-DataOFF{i}(:,3)); %AM*
        D{i,1}(:,5)=(DataOFF{i}(:,2)-DataOFF{i}(:,4)); %AP*
        D{i,1}(:,6)=(DataOFF{i}(:,3)-DataOFF{i}(:,4)); %MP
        D{i,2}(:,1)=(DataON{i}(:,1)-DataON{i}(:,2)); %CA*
        D{i,2}(:,2)=(DataON{i}(:,1)-DataON{i}(:,3)); %CM*
        D{i,2}(:,3)=(DataON{i}(:,1)-DataON{i}(:,4)); %CP*
        D{i,2}(:,4)=(DataON{i}(:,2)-DataON{i}(:,3)); %AM*
        D{i,2}(:,5)=(DataON{i}(:,2)-DataON{i}(:,4)); %AP*
        D{i,2}(:,6)=(DataON{i}(:,3)-DataON{i}(:,4)); %MP
        flag(1)=2;
        flag(2:5,i)=1;
        Nch(i)=6;

        % MEGE
        i=i+1;
        D{i,1}(:,1)=(DataOFF{i}(:,1)-DataOFF{i}(:,2)); %CM*
        D{i,1}(:,2)=(DataOFF{i}(:,1)-DataOFF{i}(:,3)); %CP*
        D{i,1}(:,3)=(DataOFF{i}(:,2)-DataOFF{i}(:,3)); %MP
        D{i,2}(:,1)=(DataON{i}(:,1)-DataON{i}(:,2)); %CM*
        D{i,2}(:,2)=(DataON{i}(:,1)-DataON{i}(:,3)); %CP*
        D{i,2}(:,3)=(DataON{i}(:,2)-DataON{i}(:,3)); %MP
        flag(1:2,i)=1;
        Nch(i)=3;

        %% Set Artefacts for these channels
        i=0;
        % BEMI
        i=i+1;
        A{i,1}{1}=[ArtOFF{i}{1}; ArtOFF{i}{2}];
        A{i,1}{2}=[ArtOFF{i}{1}; ArtOFF{i}{3}];
        A{i,1}{3}=[ArtOFF{i}{2}; ArtOFF{i}{3}];
        A{i,2}{1}=[ArtON{i}{1}; ArtON{i}{2}];
        A{i,2}{2}=[ArtON{i}{1}; ArtON{i}{3}];
        A{i,2}{3}=[ArtON{i}{2}; ArtON{i}{3}];

        % BIMA
        i=i+1;
        A{i,1}{1}=[ArtOFF{i}{1}; ArtOFF{i}{2}];
        A{i,1}{2}=[ArtOFF{i}{1}; ArtOFF{i}{3}];
        A{i,1}{3}=[ArtOFF{i}{1}; ArtOFF{i}{4}];
        A{i,1}{4}=[ArtOFF{i}{2}; ArtOFF{i}{3}];
        A{i,1}{5}=[ArtOFF{i}{2}; ArtOFF{i}{4}];
        A{i,1}{6}=[ArtOFF{i}{3}; ArtOFF{i}{4}];
        A{i,2}{1}=[ArtON{i}{1}; ArtON{i}{2}];
        A{i,2}{2}=[ArtON{i}{1}; ArtON{i}{3}];
        A{i,2}{3}=[ArtON{i}{1}; ArtON{i}{4}];
        A{i,2}{4}=[ArtON{i}{2}; ArtON{i}{3}];
        A{i,2}{5}=[ArtON{i}{2}; ArtON{i}{4}];
        A{i,2}{6}=[ArtON{i}{3}; ArtON{i}{4}];

        % GRFR, no fist available
        i=i+1;
        A{i,1}{1}=[ArtOFF{i}{1}; ArtOFF{i}{2}];
        A{i,1}{2}=[ArtOFF{i}{1}; ArtOFF{i}{3}];
        A{i,1}{3}=[ArtOFF{i}{2}; ArtOFF{i}{3}];
        A{i,2}{1}=[ArtON{i}{1}; ArtON{i}{2}];
        A{i,2}{2}=[ArtON{i}{1}; ArtON{i}{3}];
        A{i,2}{3}=[ArtON{i}{2}; ArtON{i}{3}];

        % REGE
        i=i+1;
        A{i,1}{1}=[ArtOFF{i}{1}; ArtOFF{i}{2}];
        A{i,1}{2}=[ArtOFF{i}{1}; ArtOFF{i}{3}];
        A{i,1}{3}=[ArtOFF{i}{1}; ArtOFF{i}{4}];
        A{i,1}{4}=[ArtOFF{i}{2}; ArtOFF{i}{3}];
        A{i,1}{5}=[ArtOFF{i}{2}; ArtOFF{i}{4}];
        A{i,1}{6}=[ArtOFF{i}{3}; ArtOFF{i}{4}];
        A{i,2}{1}=[ArtON{i}{1}; ArtON{i}{2}];
        A{i,2}{2}=[ArtON{i}{1}; ArtON{i}{3}];
        A{i,2}{3}=[ArtON{i}{1}; ArtON{i}{4}];
        A{i,2}{4}=[ArtON{i}{2}; ArtON{i}{3}];
        A{i,2}{5}=[ArtON{i}{2}; ArtON{i}{4}];
        A{i,2}{6}=[ArtON{i}{3}; ArtON{i}{4}];

        % RIRO
        i=i+1;
        A{i,1}{1}=[ArtOFF{i}{1}; ArtOFF{i}{2}];
        A{i,1}{2}=[ArtOFF{i}{1}; ArtOFF{i}{3}];
        A{i,1}{3}=[ArtOFF{i}{2}; ArtOFF{i}{3}];
        A{i,2}{1}=[ArtON{i}{1}; ArtON{i}{2}];
        A{i,2}{2}=[ArtON{i}{1}; ArtON{i}{3}];
        A{i,2}{3}=[ArtON{i}{2}; ArtON{i}{3}];

        % SUGI, no hold/fist available
        i=i+1;
        A{i,1}{1}=[ArtOFF{i}{1}; ArtOFF{i}{2}];
        A{i,1}{2}=[ArtOFF{i}{1}; ArtOFF{i}{3}];
        A{i,1}{3}=[ArtOFF{i}{1}; ArtOFF{i}{4}];
        A{i,1}{4}=[ArtOFF{i}{2}; ArtOFF{i}{3}];
        A{i,1}{5}=[ArtOFF{i}{2}; ArtOFF{i}{4}];
        A{i,1}{6}=[ArtOFF{i}{3}; ArtOFF{i}{4}];
        A{i,2}{1}=[ArtON{i}{1}; ArtON{i}{2}];
        A{i,2}{2}=[ArtON{i}{1}; ArtON{i}{3}];
        A{i,2}{3}=[ArtON{i}{1}; ArtON{i}{4}];
        A{i,2}{4}=[ArtON{i}{2}; ArtON{i}{3}];
        A{i,2}{5}=[ArtON{i}{2}; ArtON{i}{4}];
        A{i,2}{6}=[ArtON{i}{3}; ArtON{i}{4}];

        % MEGE
        i=i+1;
        A{i,1}{1}=[ArtOFF{i}{1}; ArtOFF{i}{2}];
        A{i,1}{2}=[ArtOFF{i}{1}; ArtOFF{i}{3}];
        A{i,1}{3}=[ArtOFF{i}{2}; ArtOFF{i}{3}];
        A{i,2}{1}=[ArtON{i}{1}; ArtON{i}{2}];
        A{i,2}{2}=[ArtON{i}{1}; ArtON{i}{3}];
        A{i,2}{3}=[ArtON{i}{2}; ArtON{i}{3}];
        
    case {'Ploc', 'raw'}
        %% Calculate phases relative to non-STN channel
        Nch = NLFP;
        for i=1:Nex
            D{i,1} = DataOFF{i};
            D{i,2} = DataON {i};
            A{i,1} = ArtOFF{i};
            A{i,2} = ArtON {i};
        end
        
end


%% Calculate wavelet transformation
Poff = NaN(Nf,max(Nch),Nex);
Pon  = NaN(Nf,max(Nch),Nex);

switch LFPtype
    case 'dLFP'
        Pcoh_off = NaN(Nf,max(Nch),max(Nch),Nex);
        Pcoh_on  = NaN(Nf,max(Nch),max(Nch),Nex);
        Pinc_off = NaN(Nf,max(Nch),max(Nch),Nex);
        Pinc_on  = NaN(Nf,max(Nch),max(Nch),Nex);
        Pvc_off  = NaN(Nf,max(Nch),max(Nch),Nex);
        Pvc_on   = NaN(Nf,max(Nch),max(Nch),Nex);
    case {'Ploc', 'raw'}
        Pcoh_off = NaN(Nf,max(Nch),Nex);
        Pcoh_on  = NaN(Nf,max(Nch),Nex);
        Pinc_off = NaN(Nf,max(Nch),Nex);
        Pinc_on  = NaN(Nf,max(Nch),Nex);
        Pvc_off  = NaN(Nf,max(Nch),Nex);
        Pvc_on   = NaN(Nf,max(Nch),Nex);
end

for i=1:Nex
    [x,Woff,coi_off,Poff(:,1:Nch(i),i)] = procdata(D{i,1}, 'freq', f, ...
        'filter', [45 55; 90 110], 'art', A{i,1});
    [x,Won, coi_on, Pon(:,1:Nch(i),i)]  = procdata(D{i,2}, 'freq', f, ...
        'filter', [45 55; 90 110], 'art', A{i,2});
    clear x
    
    switch LFPtype
        
        case 'dLFP'
            % Calculate coherent spectra of all channels w.r.t. each other
            [C_off,Wxy_off] = wave_acoh ( Woff,f );
            [C_on ,Wxy_on ] = wave_acoh ( Won ,f );
            [Pcoh_off(:,1:Nch(i),1:Nch(i),i), Pinc_off(:,1:Nch(i),1:Nch(i),i), ...
                Pvc_off(:,1:Nch(i),1:Nch(i),i)] = psd_coh( f, Woff, C_off, ...
                coi_off, 0.5, Wxy_off);
            [Pcoh_on(:,1:Nch(i),1:Nch(i),i) , Pinc_on(:,1:Nch(i),1:Nch(i),i) , ...
                Pvc_on(:,1:Nch(i),1:Nch(i),i) ] = psd_coh( f, Won , C_on , ...
                coi_on , 0.5, Wxy_on );
            
        case 'Ploc'
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
           
end


%% Determine ratio
Prel=NaN(Nf,max(Nch),Nex);
for i=1:Nex
    
    switch LFPtype
        
        case 'dLFP'
            ind=flag(:,i)==1;
            %% Ptotal
            Prel(:,ind,i) = Poff(:,ind,i)./Pon(:,ind,i);
%             %% coh-vc
%             Prel(:,ind,i)=squeeze( nanmean( (Pcoh_off(:,ind,ind,i)-Pvc_off(:,ind,ind,i)) ...
%                 ./(Pcoh_on(:,ind,ind,i)-Pvc_on(:,ind,ind,i)), 3 ) );
%             %% inc
%             Prel(:,ind,i)=squeeze( nanmean(
%             Pinc_off(:,ind,ind,i)./Pinc_on(:,ind,ind,i), 3 ) );

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
subplot(3,1,typ)
% figure
boxplot(log10(Prelband),'labels', {'1-4Hz', '4-7Hz', '7-13Hz', '13-30Hz', ...
    '13-20Hz', '20-30Hz', '30-70Hz', '60-80Hz', '250-350Hz'})
ylim([-1 1]), grid on
ylabel('log_{10}(P_{OFF}/P_{ON})')
hold all
for i=1:9
    p(i)=signrank(log(Prelband(:,i)));
    if p(i)<0.01;        
        plot([i-0.05 i+0.05],[0.9 0.9],'*k');
    elseif p(i)<0.05;
        plot(i,0.9,'*k');
    end
end
prob(:,typ)=p;
bandpow{typ}=Prelband;
relP{typ}=Prel;

clear p Prel D Prelband i x DataOFF DataON NdatOFF NdatON NLFP Patient ...
    ArtOFF ArtON LFPact ChannelOFF ChannelON RigorOFF RigorON Nex Nch flag A
end

save 'bandpower.mat' prob bandpow relP