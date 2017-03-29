%% Test mono- and bipolar PSD changes and coherence from OFF to ON


%% Load LFP data
[t, dat, CH, art, rigor, ndat, x, LFPactivity] = ...
    load_dat_from_DATA ( 'experiment', 'Ruhe', 'activity', 0, 'commontime', 1, 'noEMG', 0);

load DATA_last.mat

%% Parameters
Nf      = 23;
f       = linspace(1, 46, Nf);
sig     = 0.495;
phthres = 10;
t_smo   = 3;




%% Determine sites and analysis
% Sites with rest ON OFF, Nelec>1 and rigor>=30
i0   = [1  6 11 26 32 41 56 65]; % last entry: no EMG channels
iend = [5 10 17 31 68 46 62 67];
% % Sites with hold ON OFF and nLFP>1
% i0   = [1 3 5  9 11 21];
% iend = [2 4 6 10 23 22];
% % Sites with fist ON OFF and rigor>=30
% i0   = [1 3 8 10 18];
% iend = [2 4 9 20 19];

Nex = length(i0);

LFPch = {'C', 'L', 'A', 'M', 'P'};
EMGch = {'EDCre', 'FDLre', 'FDIre', 'EDCli', 'FDLli', 'FDIli'};
for i=1:Nex
    Channel{i}=CH{i0(i)};
    iLFP = find(ismember(Channel{i}, LFPch));
    iEMG = find(ismember(Channel{i}, EMGch));
    
    LFP_OFF{i}=dat{i0(i)}(:,iLFP);
    EMG_OFF{i}=dat{i0(i)}(:,iEMG);
    LFP_ON{i}=dat{iend(i)}(:,iLFP);
    EMG_ON{i}=dat{iend(i)}(:,iEMG);
    
    for j=1:length(iLFP)
        ArtLFP_OFF{i}{j}=art{i0(i)}{iLFP(j)};
        ArtLFP_ON{i}{j}=art{iend(i)}{iLFP(j)};
    end
    for j=1:length(iEMG)
        ArtEMG_OFF{i}{j}=art{i0(i)}{iEMG(j)};
        ArtEMG_ON{i}{j}=art{iend(i)}{iEMG(j)};
    end
        
    Patient{i}=DATA(ndat(i0(i))).patient;
    Patient2{i}=DATA(ndat(iend(i))).patient;
%     Nch(i)=size(dat{iend(i)},2);
    Nch(i) = length(iLFP);
    LFPact{i}=LFPactivity{iend(i)};
    T_off(i) = max(t{i0(i)}) - min(t{i0(i)});
    T_on(i)  = max(t{iend(i)}) - min(t{iend(i)});
    
    
    if ~strcmp(Patient{i},Patient2{i})
        error('ERROR: i0 and iend refer to different patients for i=%i!', i)
    end
end
rigor=rigor(iend);

clear i i0 iend DATA dat art CH LFPactivity t ndat pat Patient2 x



%% Calculate spectral powers: total, coherent, incoherent and volume
%% conducted
nch = max(Nch);

P_EMG_off    = NaN(Nf, 3 ,Nex);
P_tot_off    = NaN(Nf,nch,Nex);
P_nvc_off    = NaN(Nf,nch,Nex);
P_bip_off    = NaN(Nf, 6 ,Nex);
cLL_tot_off  = NaN(Nf,nch,nch,Nex);
cLL_nvc_off  = NaN(Nf,nch,nch,Nex);
cLL_bip_off  = NaN(Nf, 6 , 6 ,Nex);
cLE_tot_off  = NaN(Nf,nch,3,Nex);
% cLE_nvc_off  = NaN(Nf,nch,3,Nex);
cLE_bip_off  = NaN(Nf, 6 ,3,Nex);

P_EMG_on     = NaN(Nf, 3 ,Nex);
P_tot_on     = NaN(Nf,nch,Nex);
P_nvc_on     = NaN(Nf,nch,Nex);
P_bip_on     = NaN(Nf, 6 ,Nex);
cLL_tot_on   = NaN(Nf,nch,nch,Nex);
cLL_nvc_on   = NaN(Nf,nch,nch,Nex);
cLL_bip_on   = NaN(Nf, 6 , 6 ,Nex);
cLE_tot_on   = NaN(Nf,nch,3,Nex);
% cLE_nvc_on   = NaN(Nf,nch,3,Nex);
cLE_bip_on   = NaN(Nf, 6 ,3,Nex);

STNelec      = NaN(6,Nex);
Nbc          = NaN(Nex,1);
Nemg         = NaN(Nex,1);
BiPolCh      = cell(Nex,1);


for w0 = 10 %[10 14 20];


for i=1 %:Nex
    
    Nemg(i)=size(EMG_OFF{i},2);
    
    % Determine bipolar combinations w/o pure nonSTN channels
    combo = nchoosek(1:Nch(i),2);
    STNelec(1:size(combo,1),i) = sum(LFPact{i}(combo)>0,2);
    combo = combo( STNelec(:,i)>0, : ); %no bipolar PSD of nonSTN electrodes
    Nbc(i)= size(combo,1);
    BiPolCh{i} = Channel{i}(combo);
    
    
    
    %% OFF
    
    % LFP monopolar
    [x,W1,coi1,P_tot_off(:,1:Nch(i),i)] = procdata(LFP_OFF{i}, 'freq', f, 'w0', w0, ...
        'filter', [], 'art', ArtLFP_OFF{i});
    [C, Wxy] = wave_acoh ( W1, f );
    [cLL_tot_off(:,1:Nch(i),1:Nch(i),i), cLL_nvc_off(:,1:Nch(i),1:Nch(i),i)] ...
            = acoh1d( f, C, coi1, Wxy, phthres );
    [Pcoh, Pinc, Pvc] = psd_coh ( f, W1, C, coi1, sig, Wxy );
    P_nvc_off(:,1:Nch(i),i) = P_tot_off(:,1:Nch(i),i) - squeeze(mean(Pvc,3));
%     [cLL_tot_off(:,1:Nch(i),1:Nch(i),i), cLL_nvc_off(:,1:Nch(i),1:Nch(i),i)] ...
%         = acohere1d ( f, W1, coi1, Wxy, C, phthres );
        
%      % LFP bipolar
%      x = LFP_OFF{i}(:,combo(:,1))-LFP_OFF{i}(:,combo(:,2));
%      for j=1:Nbc(i);
%          a{j} = [ArtLFP_OFF{i}{combo(j,1)}; ArtLFP_OFF{i}{combo(j,2)}];
%      end
%      [x,W2,coi2,P_bip_off(:,1:Nbc(i),i)] = procdata(x, 'freq', f, 'w0', w0, 'filter', [], 'art', a);
%      [C, Wxy] = wave_acoh ( W2, f );
%      cLL_bip_off(:,1:Nbc(i),1:Nbc(i),i) = acoh1d( f, C, coi2, Wxy, phthres );
%  %     cLL_bip_off(:,1:Nbc(i),1:Nbc(i),i) = acohere1d ( f, W2, coi2 );
%      
%      % EMG
%      if Nemg(i)>0
%          [x,WE,coiE,P_EMG_off(:,1:Nemg(i),i)] = procdata(EMG_OFF{i}, 'freq', f, 'w0', w0, ...
%              'filter', [0 60; 90 110], 'art', ArtEMG_OFF{i}, 'rect', 1);
%          C = wave_coh ( W1, WE, f );
%          cLE_tot_off(:,1:Nch(i),1:Nemg(i),i) = coh1d( f, C, coi1, coiE );
%          C = wave_coh ( W2, WE, f );
%          cLE_bip_off(:,1:Nbc(i),1:Nemg(i),i) = coh1d( f, C, coi2, coiE );
%  %         cLE_tot_off(:,1:Nch(i),1:Nemg(i),i) = cohere1d( f, W1, WE, coi1, coiE );
%  %         cLE_bip_off(:,1:Nbc(i),1:Nemg(i),i) = cohere1d( f, W2, WE, coi2, coiE );
%      end
%      
%      
%      %% ON
%      
%      % LFP monopolar
%      [x,W1,coi1,P_tot_on(:,1:Nch(i),i)] = procdata(LFP_ON{i}, 'freq', f, 'w0', w0, ...
%          'filter', [], 'art', ArtLFP_ON{i});
%      [C, Wxy] = wave_acoh ( W1, f );
%      [cLL_tot_on(:,1:Nch(i),1:Nch(i),i), cLL_nvc_on(:,1:Nch(i),1:Nch(i),i)] ...
%              = acoh1d( f, C, coi1, Wxy, phthres );
%      [Pcoh, Pinc, Pvc] = psd_coh ( f, W1, C, coi1, sig, Wxy );
%      P_nvc_on(:,1:Nch(i),i) = P_tot_on(:,1:Nch(i),i) - squeeze(mean(Pvc,3));
%  %     [cLL_tot_on(:,1:Nch(i),1:Nch(i),i), cLL_nvc_on(:,1:Nch(i),1:Nch(i),i)] ...
%  %         = acohere1d ( f, W1, coi1, Wxy, C, phthres );
%          
%      % LFP bipolar
%      x = LFP_ON{i}(:,combo(:,1))-LFP_ON{i}(:,combo(:,2));
%      for j=1:Nbc(i);
%          a{j} = [ArtLFP_ON{i}{combo(j,1)}; ArtLFP_ON{i}{combo(j,2)}];
%      end
%      [x,W2,coi2,P_bip_on(:,1:Nbc(i),i)] = procdata(x, 'freq', f, 'w0', w0, 'filter', [], 'art', a);
%      [C, Wxy] = wave_acoh ( W2, f );
%      cLL_bip_on(:,1:Nbc(i),1:Nbc(i),i) = acoh1d( f, C, coi2, Wxy, phthres );
%  %     cLL_bip_on(:,1:Nbc(i),1:Nbc(i),i) = acohere1d ( f, W2, coi2 );
%      
%      % EMG
%      if Nemg(i)>0
%          [x,WE,coiE,P_EMG_on(:,1:Nemg(i),i)] = procdata(EMG_ON{i}, 'freq', f, 'w0', w0, ...
%              'filter', [0 60; 90 110], 'art', ArtEMG_ON{i}, 'rect', 1);
%          C = wave_coh ( W1, WE, f );
%          cLE_tot_on(:,1:Nch(i),1:Nemg(i),i) = coh1d( f, C, coi1, coiE );
%          C = wave_coh ( W2, WE, f );
%          cLE_bip_on(:,1:Nbc(i),1:Nemg(i),i) = coh1d( f, C, coi2, coiE );
%  %         cLE_tot_on(:,1:Nch(i),1:Nemg(i),i) = cohere1d( f, W1, WE, coi1, coiE );
%  %         cLE_bip_on(:,1:Nbc(i),1:Nemg(i),i) = cohere1d( f, W2, WE, coi2, coiE );
%      end
    
    
end

%  clear LFP_OFF EMG_OFF LFP_ON EMG_ON ArtOFF ArtON x a W1 W2 WE Wxy C coi1 coi2 coiE combo
%  
%  
%  
%  %% Total power in frequency bands
%  delta = find(f>=1  & f<=4 );
%  theta = find(f>=4  & f<=7 );
%  alpha = find(f>=7  & f<=13);
%  beta  = find(f>=13 & f<=30);
%  beta1 = find(f>=13 & f<=20);
%  beta2 = find(f>=20 & f<=30);
%  gamma = find(f>=30 & f<=70);
%  % HFO   = find(f>=250 & f<=350);
%  
%  x=P_nvc_off(:,:);
%  PL_nvc_off(:,1) = squeeze(sum(x(delta,:),1));
%  PL_nvc_off(:,2) = squeeze(sum(x(theta,:),1));
%  PL_nvc_off(:,3) = squeeze(sum(x(alpha,:),1));
%  PL_nvc_off(:,4) = squeeze(sum(x(beta,:),1));
%  PL_nvc_off(:,5) = squeeze(sum(x(beta1,:),1));
%  PL_nvc_off(:,6) = squeeze(sum(x(beta2,:),1));
%  PL_nvc_off(:,7) = squeeze(sum(x(gamma,:),1));
%  
%  x=P_nvc_on(:,:);
%  PL_nvc_on(:,1) = squeeze(sum(x(delta,:),1));
%  PL_nvc_on(:,2) = squeeze(sum(x(theta,:),1));
%  PL_nvc_on(:,3) = squeeze(sum(x(alpha,:),1));
%  PL_nvc_on(:,4) = squeeze(sum(x(beta,:),1));
%  PL_nvc_on(:,5) = squeeze(sum(x(beta1,:),1));
%  PL_nvc_on(:,6) = squeeze(sum(x(beta2,:),1));
%  PL_nvc_on(:,7) = squeeze(sum(x(gamma,:),1));
%  
%  x=P_bip_off(:,:);
%  PL_bip_off(:,1) = squeeze(sum(x(delta,:),1));
%  PL_bip_off(:,2) = squeeze(sum(x(theta,:),1));
%  PL_bip_off(:,3) = squeeze(sum(x(alpha,:),1));
%  PL_bip_off(:,4) = squeeze(sum(x(beta,:),1));
%  PL_bip_off(:,5) = squeeze(sum(x(beta1,:),1));
%  PL_bip_off(:,6) = squeeze(sum(x(beta2,:),1));
%  PL_bip_off(:,7) = squeeze(sum(x(gamma,:),1));
%  
%  x=P_bip_on(:,:);
%  PL_bip_on(:,1) = squeeze(sum(x(delta,:),1));
%  PL_bip_on(:,2) = squeeze(sum(x(theta,:),1));
%  PL_bip_on(:,3) = squeeze(sum(x(alpha,:),1));
%  PL_bip_on(:,4) = squeeze(sum(x(beta,:),1));
%  PL_bip_on(:,5) = squeeze(sum(x(beta1,:),1));
%  PL_bip_on(:,6) = squeeze(sum(x(beta2,:),1));
%  PL_bip_on(:,7) = squeeze(sum(x(gamma,:),1));
%  
%  for i=1:7
%      sig_nvc(i) = signrank (PL_nvc_off(:,i)./PL_nvc_on(:,i),1);
%      sig_bip(i) = signrank (PL_bip_off(:,i)./PL_bip_on(:,i),1);
%  end
%  
%  
%  %% Reshape matrices and one-sided coherence
%  % n1=1; n2=1;
%  % for k=1:Nex
%  %     for j=1:Nbc(k);
%  %         for i=j+1:Nbc(k)
%  %             coh_LFP_on(:,n1)      = c1d_LFP_on(:,i,j,k);
%  %             coh_LFP_off(:,n1)     = c1d_LFP_off(:,i,j,k);
%  %             coh_LFP_novc_on(:,n1) = c1d_LFP_novc_on(:,i,j,k);
%  %             coh_LFP_novc_off(:,n1)= c1d_LFP_novc_off(:,i,j,k);
%  %             n1=n1+1;
%  %         end
%  %         for i=1:Nemg(k)
%  %             coh_LE_on(:,n2)      = c1d_LE_on(:,j,i,k);
%  %             coh_LE_off(:,n2)     = c1d_LE_off(:,j,i,k);
%  %             n2=n2+1;
%  %         end
%  %     end
%  % end
%  % clear n1 n2 i j k
%  % P_LFP_off = P_LFP_off(:,:);
%  % P_LFP_on  = P_LFP_on(:,:);
%  % P_EMG_off = P_EMG_off(:,:);
%  % P_EMG_on  = P_EMG_on(:,:);
%  % 
%  % 
%  % % % STN power relative to nonSTN-reference-electrode
%  % % j = STNelec==1;
%  % % pow1_off = P_off(:,j);
%  % % pow1_on  = P_on(:,j);
%  % % % STN power relative to STN-reference-electrode
%  % % j = STNelec==2;
%  % % pow2_off = P_off(:,j);
%  % % pow2_on  = P_on(:,j);
%  % % 
%  % % P_off = [pow1_off pow2_off];
%  % % P_on  = [pow1_on  pow2_on ];
%  % 
%  %     
%  % 
%  % %% Significant power of coherence changes?
%  % for i=1:30;
%  % %     sig_pow(i,1)=signrank( pow1_off(i,:) ./ pow1_on(i,:), 1 );
%  % %     sig_pow(i,2)=signrank( pow2_off(i,:) ./ pow2_on(i,:), 1 );
%  % %     sig_pow(i,3)=signrank( P_off(i,:) ./ P_on(i,:), 1 );
%  % %     sig_coh(i,1)=signrank( coh_off(i,:) - coh_on(i,:), 0 );
%  % %     sig_coh(i,2)=signrank( coh_novc_off(i,:) - coh_novc_on(i,:), 0 );
%  %     sig_pow(i,1)=signrank( P_LFP_off(i,:) ./ P_LFP_on(i,:), 1 );
%  %     sig_pow(i,2)=signrank( P_EMG_off(i,:) ./ P_EMG_on(i,:), 1 );
%  %     sig_coh(i,1)=signrank( coh_LFP_novc_off(i,:) - coh_LFP_novc_on(i,:), 0 );
%  %     sig_coh(i,2)=signrank( coh_LE_off(i,:) - coh_LE_on(i,:), 0 );
%  % 
%  % end
%  % 
%  % 
%  % % %% Plot results
%  % % figure
%  % % subplot(3,3,1)
%  % % loglog( f, pow3_off ./ pow3_on )
%  % % title('Total power change during wrist'), ylim([0.05 20])
%  % % ylabel('P_{OFF}/P_{ON}'), xlabel('f [Hz]')
%  % % hold all
%  % % plot( f, nanmean( mpsd_tot_off ./ mpsd_tot_on, 2 ), '--k', 'linew', 2 )
%  % % i= sig_tot<0.05; if any(i); plot(f(i),10,'*k','markers',5); end
%  % % i= sig_tot<0.01; if any(i); plot(f(i),13,'+k','markers',5); end
%  % % 
%  % % subplot(2,1,2)
%  % % loglog(f,(mpsd_tot_off-mpsd_vc_off) ./ ...
%  % %     (mpsd_tot_on-mpsd_vc_on))
%  % % title('Change of power w/o vol. cond. during wrist'), ylim([0.05 20])
%  % % ylabel('P_{OFF}/P_{ON}'), xlabel('f [Hz]')
%  % % hold all
%  % % plot(f, nanmean((mpsd_tot_off-mpsd_vc_off) ./ ...
%  % %     (mpsd_tot_on-mpsd_vc_on),2), '--k', 'linew', 2)
%  % % i= sig_novc<0.05; if any(i); plot(f(i),10,'*k','markers',5); end
%  % % i= sig_novc<0.01; if any(i); plot(f(i),13,'+k','markers',5); end
%  
%  
%  
%  %% Plot results for each patient
%  for i=1:Nex
%      fig = figure('PaperUnits','centimeters', 'visi', 'off', ...
%          'PaperPosition', [1 1 19 27], ...
%          'PaperSize', [21 29.7], ...
%          'Position',[0 0 21 29.7]/2);
%      subplot(4,2,1)
%      semilogy(f,P_nvc_off(:,:,i),'-',f,P_nvc_on(:,:,i),'--')
%      title('Mono - vol.cond.'), xlim([min(f) max(f)])%, ylim([1e-8 1e-6])
%      xlabel('f [Hz'), ylabel('PSD [V^2/Hz]')
%      subplot(4,2,2)
%      semilogy(f,P_bip_off(:,:,i),'-',f,P_bip_on(:,:,i),'--')
%      title('Bipolar'), xlim([min(f) max(f)])%, ylim([1e-8 1e-6])
%      xlabel('f [Hz'), ylabel('PSD [V^2/Hz]')
%      subplot(4,2,3)
%      semilogy(f,P_tot_off(:,:,i),'-',f,P_tot_on(:,:,i),'--')
%      title('mono'), xlim([min(f) max(f)])%, ylim([1e-8 1e-6])
%      xlabel('f [Hz'), ylabel('PSD [V^2/Hz]')
%      subplot(4,2,4)
%      semilogy(f,P_EMG_off(:,:,i),'-',f,P_EMG_on(:,:,i),'--')
%      title('EMG'), xlim([min(f) max(f)])%, ylim([1e-11 1e-8])
%      xlabel('f [Hz'), ylabel('PSD [V^2/Hz]')
%      for j=1:Nch(i)-1
%          subplot(4,2,5)
%          plot(f,cLL_nvc_off(:,j+1:end,j,i),'-',f,cLL_nvc_on(:,j+1:end,j,i),'--')
%          title('Coh Mono - vol.cond.'), xlim([min(f) max(f)]), ylim([0 1])
%          xlabel('f [Hz'), ylabel('Coherence'), hold all
%          subplot(4,2,6)
%          plot(f,cLL_bip_off(:,j+1:end,j,i),'-',f,cLL_bip_on(:,j+1:end,j,i),'--')
%          title('Coh Bipolar'), xlim([min(f) max(f)]), ylim([0 1])
%          xlabel('f [Hz'), ylabel('Coherence'), hold all
%      end
%      for j=1:Nemg(i)
%          subplot(4,2,7)
%          plot(f,cLE_tot_off(:,:,j,i),'-',f,cLE_tot_on(:,j+1:end,j,i),'--')
%          title('Coh LE Mono'), xlim([min(f) max(f)]), ylim([0 1])
%          xlabel('f [Hz'), ylabel('Coherence'), hold all
%          subplot(4,2,8)
%          plot(f,cLE_bip_off(:,:,j,i),'-',f,cLE_bip_on(:,j+1:end,j,i),'--')
%          title('Coh LE Bipolar'), xlim([min(f) max(f)]), ylim([0 1])
%          xlabel('f [Hz'), ylabel('Coherence'), hold all
%      end
%      saveas(fig, ['pow_coh_P' num2str(i,'%02i') 'w0_' num2str(w0,'%02i') '.pdf']);
%      close(fig)
%  end

end