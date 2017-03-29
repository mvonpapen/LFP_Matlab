%%  Print a booklet with all time series and scalograms
%%
%%
clear all
load /afs/.geo.uni-koeln.de/usr/neuro/DATA/DATA_last.mat
list=dir('/afs/.geo.uni-koeln.de/usr/neuro/DATA/*S*.mat');
whereto = '/afs/.geo.uni-koeln.de/usr/neuro/Databooklet/';


LFPch={'C', 'L', 'A', 'M', 'P'};
EMGch={'EDCre', 'EDCli', 'FDIre', 'FDIli', 'FDLre', 'FDLli'};
limL=[-9 -4]; %in log10(V^2/Hz)
limE=[-11 -6];

f=logspace(0,2.6,40);
dt=1/2500;

ind=1;
for n=1:length(DATA)
    
    if isempty(DATA(n).site)
        continue
    end
    
    %% Load data
    patient=DATA(n).patient;
    site=DATA(n).site;
    apo=DATA(n).Apo;
    LFPact=DATA(n).LFPactivity;
    exp=DATA(n).experiment;
    if strcmp(exp, 'Apo') || strcmp(exp, 'Placebo')
        continue
    end
    [t, dat, ndat, rigor, CH, art] = load_dat_from_DATA ( DATA, ...
        'field1', 'site', 'value1', site, 'value2', patient, 'apomorphine', apo );
    if isempty(t{1}{1})
        continue
    end
    art=art{1};
    t=t{1};
    dat=dat{1};
    CH=CH{1};
    nch=length(CH);
    
    
    %% Convert artefact times into data points
    for i=1:length(art{1})
        iart{i}=[];
        iart2{i}=[];
        [a1 a2]=size(art{1}{i});
        for j=1:a1
            i0=max([1 find(t{i}<=art{1}{i}(j,1),1,'last')]);
            iend=min([length(t{i}) find(t{i}>=art{1}{i}(j,2),1,'first')]);
            iart{i}=[iart{i}; i0 iend];
            iart2{i}=[iart2{i} i0:iend]; % used later to set data to NaN
        end
    end
    
    
    %% LFP data
    iLFP=find(ismember(CH,LFPch));
    nLFP=length(iLFP);
    dLFP=cell(nLFP,1);
    WLFP=cell(nLFP,1);
    coiLFP=cell(nLFP,1);
    PLFP=zeros(length(f),nLFP);
    for i=1:nLFP
        [dLFP{i}, WLFP{i}, coiLFP{i}, PLFP(:,i)]...
                = procdata (dat{iLFP(i)}, 'freq', f, 'filter', [], ...
                  'art', art{iLFP(i)});
    end
    iLFP2=find(ismember(DATA(n).channel,LFPch));
    LFPact=LFPact(iLFP2);
    clear iLFP2
    
%     %% EMG data
%     iEMG=find(ismember(CH,EMGch));
%     nEMG=length(iEMG);
%     if nEMG>3
%         iEMG=iEMG(1:3);
%         nEMG=length(iEMG);
%     end
%     dEMG=cell(nEMG,1);
%     WEMG=cell(nEMG,1);
%     coiEMG=cell(nEMG,1);
%     PEMG=NaN(length(f),nEMG);
%     for i=1:nEMG;
%         [dEMG{i},WEMG{i},coiEMG{i},PEMG(:,i)]...
%             =procdata(dat{iEMG(i)},'filter', [],...
%             'rect', 0, 'freq', f, 'art', art{iEMG(i)});
%     end
    
    %% Calculate power in spectral bands
    delta=round(100*nanmean(PLFP(f>=1 & f<=4,:))*1e6)/100;
    theta=round(100*nanmean(PLFP(f>=4 & f<=7,:))*1e6)/100;
    alpha=round(100*nanmean(PLFP(f>=7 & f<=13,:))*1e6)/100;
    beta1=round(1000*nanmean(PLFP(f>=13 & f<=20,:))*1e6)/1000;
    beta2=round(1000*nanmean(PLFP(f>=20 & f<=30,:))*1e6)/1000;
    gamma=round(1e6*nanmean(PLFP((f>=30 & f<=45) | (f>=55 & f<=70),:))*1e6);
    HF=round(1e6*nanmean(PLFP(f>=250 & f<=350,:))*1e6);
    
    
    %% Page 1
    %% Plot time series
    maxLFP=0; maxEMG=0; minmaxT=[];
    for i=1:nLFP
%         j=setxor([1:length(dLFP{i})], iart2{iLFP(i)});
%         maxLFP=max([maxLFP abs(dLFP{i}(j))']);
        maxLFP=max([maxLFP abs(dLFP{i})']);
        minmaxT = minmax( [minmaxT t{iLFP(i)}] );
    end
    for i=1:nEMG
%         j=setxor([1:length(dEMG{i})], iart2{iEMG(i)});
%         maxEMG=max([maxEMG abs(dEMG{i}(j))']);
        maxEMG=max([maxEMG abs(dEMG{i})']);
        minmaxT = minmax( [minmaxT t{iEMG(i)}] );
    end
    fig1=figure('Papertype', 'A4', ...
        'Paperunits', 'normalized', 'PaperPosition', [0 0 1 1], ...
        'PaperPositionmode', 'manual', 'Visible', 'off', ...
        'PaperOrientation', 'portrait');
    for i=1:nLFP
        subplot(4,2,2*i-1)
        plot(t{iLFP(i)},dLFP{i}*1e3, '-k')
        xlim( minmaxT )
        ylim([-maxLFP maxLFP]*1e3)
        ylabel('LFP [mV]')
        title(CH{iLFP(i)})
    end
    for i=1:nEMG
        subplot(4,2,2*i)
        plot(t{iEMG(i)},dEMG{i}*1e3, '-k')
        xlim( minmaxT )
        ylim([-maxEMG maxEMG]*1e3)
        ylabel('EMG [mV]')
        title(CH{iEMG(i)})
    end
    
    %% Write mean band power
    annotation(fig1,'textbox', [0.55 0.1 0.5 0.17], 'String', ...
        {['Mean spectral power in LFP''s'], ...
         ['\delta = ' num2str(delta) ' [mV^2/Hz]'], ...
         ['\theta = ' num2str(theta) ' [mV^2/Hz]'], ...
         ['\alpha = ' num2str(alpha) ' [mV^2/Hz]'], ...
         ['\beta_1 = ' num2str(beta1) ' [mV^2/Hz]'], ...
         ['\beta_2 = ' num2str(beta2) ' [mV^2/Hz]'], ...
         ['\gamma = ' num2str(gamma) ' [\muV^2/Hz]'], ...
         ['HF = ' num2str(HF) ' [\muV^2/Hz]']}, ...
         'LineStyle', 'none');
    
    
    %% Page 2
    %% Plot scalograms
    fig2=figure('Papertype', 'A4', ...
        'Paperunits', 'normalized', 'PaperPosition', [0 0 1 1], ...
        'PaperPositionmode', 'manual', 'Visible', 'off', ...
        'PaperOrientation', 'portrait');
    
    for i=1:nLFP;
        subplot(4,2,2*i-1)
        plot_sca(t{iLFP(i)}, f, 2*dt*abs(WLFP{i}).^2, 'coi', coiLFP{i}, 'lim', limL, ...
            'shading', 'flat');
        xlim( minmax( t{iLFP(i)} ) )
        ylim( minmax( f(:)' ) )
        ylabel('f [Hz]', 'Fontsize', 8)
        title(['LFP ' CH{iLFP(i)}], 'Fontsize', 8)
    end
    
    %% Plot EMG
    for i=1:nEMG;
        subplot(4,2,2*i)
        plot_sca(t{iEMG(i)}, f, 2*dt*abs(WEMG{i}).^2, 'coi', coiEMG{i}, 'lim', limE, ...
            'shading', 'flat');
        xlim( minmax( t{iEMG(i)} ) )
        ylim( minmax( f(:)' ) )
        ylabel('f [Hz]', 'Fontsize', 8)
        title(CH{iEMG(i)}, 'Fontsize', 8)
    end
    
    %% Plot power spectra
    subplot(4,2,8)
    if ~isempty(PEMG)
        loglog(f,PLFP,f,PEMG)
    else
        loglog(f,PLFP)
    end
    legend(CH, 'Location', [0.9 0.25 0.01 0.01])
    legend boxoff
    xlim(minmax(f))
    ylim(minmax([PLFP(:); PEMG(:)]'))
    xlabel('f [Hz]')
    ylabel('PSD [mV^2/Hz]')
        
    
    %% Plot annotations and save
    annotation(fig1,'textbox', [0.1 0.95 0.8 0.02],...
        'String', ['\bf Name: ' patient num2str(site)...
            ', Experiment: ' exp ', Rigor: ', num2str(rigor)],...
        'HorizontalAlign', 'center', 'LineStyle', 'none', 'Fontsize', 14);
    annotation(fig2,'textbox', [0.25 0.95 0.5 0.02],...
        'String', ['\bf Name: ' patient num2str(site)...
            ', LFPact: ', num2str(LFPact(iLFP)), ', Apo: ', num2str(apo)],...
        'HorizontalAlign', 'center', 'LineStyle', 'none', 'Fontsize', 14);
    set(fig2,'renderer','painters');
    saveas(fig1, [whereto 'test-' num2str(2*ind-1,'%03u') '.pdf']);
    saveas(fig2, [whereto 'test-' num2str(2*ind,'%03u') '.pdf']);
    close(fig1);
    close(fig2);
    ind=ind+1
end