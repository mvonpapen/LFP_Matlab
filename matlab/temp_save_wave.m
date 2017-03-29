%% Load data from variable DATA and put it in the pipe!
clear all
load '/afs/geo.uni-koeln.de/usr/neuro/DATA/DATA_last.mat'
[t, dat, ndat, rigor, CHlab, art]=...
    load_dat_from_DATA ( DATA, 'value1', 'Ruhe', 'activity', 1);

% FF=[1 4; 4 7; 7 13; 13 30; 30 70; 60 80; 250 350];

% for ii=1:length(FF)

%% Set variables
nEMG=3;
nf=50;
t_smo=3; % number of scales used to smooth
sc_smo=0; % number of neighboring scales to smooth in freq-domain
ns=length(ndat);
LFP=cell(ns,1);
EMG=cell(ns,1);
WLFP=cell(ns,1);
coiLFP=cell(ns,1);
PLFP=NaN(nf,5,ns);
PEMG=NaN(nf,nEMG,ns);
WEMG=cell(ns,1);
coiEMG=cell(ns,1);
T=cell(ns,1);
% f=logspace(log10(FF(ii,1)),log10(FF(ii,2)),nf);
f=logspace(0,2.7,nf);
dt=1/2500;

%% Define LFP and EMG channels
LFPch={'C', 'L', 'P', 'A', 'M'};
EMGch={'EDCre', 'EDCli', 'FDIre', 'FDIli', 'FDLre', 'FDLli'};

%% Brain frequency bands
delta=find(f>=1 & f<4);
theta=find(f>=4 & f<7);
alpha=find(f>=7 & f<13); 
beta=find(f>=13 & f<30);
gamma=find(f>=30 & f<70);
HF=find(f>=250 & f<350);

for n=1:ns
    
    site(n)=DATA(ndat(n)).site;
    
    %% Determine common time vector for LFP and EMG
    for i=1:length(t{n})
        ti(i,:)=minmax(t{n}{i});
    end
    ti=[max(ti(:,1)) min(ti(:,2))];
    T{n} = t{n}{i}( t{n}{i}>=ti(1) & t{n}{i}<=ti(2) );
    
    %% LFP data
    iLFP=find(ismember(CHlab{n},LFPch));
    nLFP(n)=length(iLFP);
    k=1;
    for i=iLFP
        ni =  t{n}{i}>=ti(1) & t{n}{i}<=ti(2);
        tmp(:,k)=dat{n}{i}(ni);
        art{n}{i}=art{n}{i}-ni(1)+1;
        k=k+1;
    end
    [LFP{n}, WLFP{n}, coiLFP{n}, PLFP(:,1:nLFP(n),n)]...
            = procdata (tmp, 'freq', f, 'filter', [45 55; 90 110], ...
                        'art', art{n});
    clear tmp ni i k
    
    %% dLFP data
%     nLFP(n)=length(t{n})-nEMG;
%     if nLFP(n)<=1
%         continue
%     end
%     ni2 =  t{n}{nLFP(n)}>=ti(1) & t{n}{nLFP(n)}<=ti(2);
%     for i=1:nLFP(n)-1
%         ni1 =  t{n}{i}>=ti(1) & t{n}{i}<=ti(2);
%         [datdf{n}(:,i), WdLFP{n}(:,:,i), coidLFP{n}(:,i), PdLFP(:,i,n)]...
%             = procdata (dat{n}{i}(ni1)-dat{n}{nLFP(n)}(ni2),...
%                 'freq', f, 'filter', []); 
%     end
    
    %% EMG data
    iEMG=find(ismember(CHlab{n},EMGch));
    nEMG(n)=length(iEMG);
    if nEMG(n)>0
        k=1;
        for i=iEMG;
            ni =  t{n}{i}>=ti(1) & t{n}{i}<=ti(2);
            tmp(:,k)=dat{n}{i}(ni);
            k=k+1;
        end
        [EMG{n}, WEMG{n}, coiEMG{n}, PEMG(:,1:nEMG(n),n)]...
                =procdata(tmp,'filter', [0 60; 90 110], 'freq', f, ...
                'art', art{n}, 'rect', 1);
        clear tmp ni i k
    end
    
    %% Coherency
    if nLFP(n)>0 && nEMG(n)>0
        [CLE{n}, WLE{n}] = wave_coh(WLFP{n}, WEMG{n}, f, t_smo, sc_smo);
%     [CLL{n}, WLL{n}] = wave_coh(WLFP{n}, WLFP{n},...
%         f, t_smo, sc_smo);
%     [dCLE{n}, dWLE{n}] = wave_coh(WdLFP{n}, WEMG{n},f);
    end
   

    %% Phases
%     clear i
%     for j=1:nLFP(n)
%         for k=1:nLFP(n) %nEMG
%             % LFP - EMG
%             WLE{n}(:,:,j,k)=squeeze(WLFP{n}(:,:,j)...
%                 .*conj(WEMG{n}(:,:,k)));
%             phase=phases(f,WLE{n}(:,:,j,k),minmax(f),...
%                 'coi',min([coiLFP{n}(:,j)';coiEMG{n}(:,k)']));
% %             % LFP - LFP
% %             if j==k
% %                 continue
% %             end
% %             WLL{n}(:,:,j,k)=squeeze(WLFP{n}(:,:,j)...
% %                 .*conj(WLFP{n}(:,:,k)));
% %             phase=phases(f,WLL{n}(:,:,j,k),minmax(f),...
% %                 'coi',min([coiLFP{n}(:,j)';coiLFP{n}(:,k)']));
%             [Hph(:,j,k,ii,n), bin]=hist(phase/pi*180,[-170:20:170]);
%             Hph(:,j,k,ii,n)=Hph(:,j,k,ii,n)/sum(Hph(:,j,k,ii,n));
%             jj=~isnan(phase);
%             PLV(j,k,ii,n)=abs(sum(exp(i*phase(jj))))/length(phase(jj));
%         end
end
% end
clear ni nia nib i j n ti ni2



%% Plot Routines

% %% Coherency-Phase Scalogram
% for n=1:ns
%     fig1=figure('Papertype', 'A4', ...
%         'Paperunits', 'normalized', 'PaperPosition', [0 0 1 1], ...
%         'PaperPositionmode', 'manual', 'Visible', 'off', ...
%         'PaperOrientation', 'landscape');
%     k=1;
%     for i=1:nLFP(n);
%         for j=1:3; 
%             subplot(nLFP(n),3,k)
%             plot_coh_phase( T{n},f,squeeze(CLE{n}(:,:,i,j)),WLE{n}(:,:,i,j),...
%               'ahead',2, 'titel', [CHlab{n}{i} '<->' CHlab{n}{j+nLFP(n)}],...
%               'coi', min([coiLFP{n}(:,i)';coiEMG{n}(:,j)']), 'sig', 0.5 );
%             k=k+1;
%         end
%     end
%     set(gcf, 'renderer', 'painters')
%     saveas(fig1,['coh_' DATA(ndat(n)).patient '_S' num2str(DATA(ndat(n)).site) ...
%         '.eps']);
%     close(fig1)
% end

% %% Coherency-Phase Scalogram for dLFP
% for n=1:ns
%     figure
%     k=1;
%     for i=1:2;
%         for j=1:3; 
%             subplot(2,3,k)
%             plot_coh_phase(T{n},f,squeeze(dCLE{n}(:,:,i,j)),dWLE{n}(:,:,i,j),...
%                 'ahead',2, 'sig', 0.5);
%             k=k+1;
%         end
%     end
% end


% %% Plot resulting PSD
% figure
% for i=1:max(nLFP)
%     subplot(1,max(nLFP),i)
%     loglog(f,squeeze(PLFP(:,i,11:15)))
% %     semilogx(f,PLFP(:,i,2)./PLFP(:,i,1),f,PLFP(:,i,3)./PLFP(:,i,1),...
% %         f,PLFP(:,i,4)./PLFP(:,i,1),f,PLFP(:,i,5)./PLFP(:,i,1));
% %     semilogx( f,squeeze(PLFP(:,i,2:5)-PLFP(:,i,1:4)),'-',...
% %         f,squeeze(PLFP(:,i,7:10)-PLFP(:,i,6:9)),'--' );
%     xlabel('f [Hz]')
%     ylabel('PSD V^2/Hz]')
%     title(CHlab{1}{i})
% end

% figure
% for i=1:2;
%     subplot(2,1,i)
%     semilogx(f,dPLFP(:,i,2)./dPLFP(:,i,1),f,dPLFP(:,i,3)./dPLFP(:,i,1),...
%         f,dPLFP(:,i,4)./dPLFP(:,i,1),f,dPLFP(:,i,5)./dPLFP(:,i,1)); 
%     xlabel('f [Hz]')
%     ylabel('dP_i/dP_1')
% end
% 
% figure
% for i=1:3;
%     subplot(3,1,i)
%     loglog(f,squeeze(PEMG(:,i,:)))
% %     semilogx(f,PEMG(:,i,2)./PEMG(:,i,1),f,PEMG(:,i,3)./PEMG(:,i,1),...
% %         f,PEMG(:,i,4)./PEMG(:,i,1),f,PEMG(:,i,5)./PEMG(:,i,1)); 
%     xlabel('f [Hz]')
%     ylabel('PSD V^2/Hz]')
% %     ylabel('P_i/P_1')
% end

% %% 1D Coherency spectra
% figure
% C1d=NaN(nf,ns,max(nLFP),3);
% for s=1:ns;
%     for i=1:nLFP(s);
%         for j=1:3;
%             C1d(:,s,i,j)=nanmean(coi2nan(f,squeeze(CLE{s}(:,:,i,j)),...
%                 min([coiLFP{s}(:,i)';coiEMG{s}(:,j)'])),2);
%         end
%     end
% end
% k=1;
% for i=1:nLFP(s)
%     for j=1:3
%         subplot(max(nLFP),3,k)
%         plot(f,C1d(:,:,i,j), '-', [0 1e3], [0.21 0.21], '--k')%,'-',f,C1d(:,6:10,i,j),'--')
%         k=k+1;
%     end
% end

% %% Plot all PSD for almost no (a) and stronger (b) rigor amelioration
% a=find(rigor<=10)
% b=find(rigor>10)
% subplot(2,1,1)
% for i=1:5; loglog(f,squeeze(PLFP(:,i,a))), hold all, end
% subplot(2,1,2)
% for i=1:5; loglog(f,squeeze(PLFP(:,i,b))), hold all, end


% %% Plot ratio of power ON/OFF in specified freq. band
% % Begin and end of experiments
% i0   = [1  6 11 17 25 31 36 41 50];
% iend = [5 10 16 24 30 34 40 49 56];
% 
% % For LFP channels
% k=1;
% for i=1:length(i0)
%     for j=1:nLFP(i0(i))
%         P0(k,ii)=nanmean(PLFP(:,j,i0(i)));
%         Pend(k,ii)=nanmean(PLFP(:,j,iend(i)));
%         k=k+1;
%     end
% end
% % % For EMG channels
% % k=1;
% % for i=1:length(i0)
% %     for j=1:nEMG(i0(i))
% %         P0_EMG(k,ii)=mean(PEMG(:,j,i0(i)));
% %         Pend_EMG(k,ii)=mean(PEMG(:,j,iend(i)));
% %         k=k+1;
% %     end
% % end
% % clear PLFP datf
% end
% figure, boxplot(Pend./P0,'labels',...
%     {'1-4Hz', '4-7Hz', '7-13Hz', '13-30Hz', '30-70Hz', '60-80Hz', '250-350Hz'})
% % figure, boxplot(Pend_EMG./P0_EMG,'labels',...
% %     {'1-4Hz', '4-7Hz', '7-13Hz', '13-30Hz', '30-70Hz', '60-80Hz', '250-350Hz'})

% %% Boxplot PLV changes
% i0=[1 6 11 16 24 30];
% iend=[5 10 15 23 29 38];
% for ii=1:6
%     HH=NaN(1,max(iend-i0)+1);
%     m=1;
%     for i=1:length(i0)-1;
%         for j=1:3; 
%             for k=1:3;
%                 HH(m,1:(iend(i)-i0(i)+1))=...
%                     squeeze(PLV(j,k,ii,i0(i):iend(i)));
%                 m=m+1;
%             end
%         end
%     end
%     i=find(HH==0);
%     HH(i)=NaN;
%     figure
%     boxplot(HH)
% end