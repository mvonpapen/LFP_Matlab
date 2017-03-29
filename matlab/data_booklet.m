%%  Print a booklet with all time series and scalograms
%%
%%
clear all
load /afs/.geo.uni-koeln.de/usr/neuro/DATA/DATA_last.mat
list=dir('/afs/.geo.uni-koeln.de/usr/neuro/DATA/*S*.mat');
whereto = '/afs/.geo.uni-koeln.de/usr/neuro/Databooklet/';
ind=1;
for n=1:length(list)
    
%     if isempty(DATA(n).interval)
%         continue
%     end
    
    
    %% Load data
    load(['/afs/.geo.uni-koeln.de/usr/neuro/DATA/' list(n).name])
    f=logspace(0,2.6,40);
    dt=1/2500;
    t=[1:length(data)-1]/2500;
    data=data(2:end,1:end-1);
    nch=length(CH)-1;
    patient=list(n).name(1:6);
    ndat1 = structfind(DATA,'patient',patient);
    ndat2 = structfind(DATA,'site',site);
    ndat = intersect( ndat1, ndat2 );
    if isempty(ndat)
        continue
    end
    clear ndat1 ndat2
    rigor=DATA(ndat).rigor;
    apo=DATA(ndat).Apo;
    LFPact=DATA(ndat).LFPactivity;
    exp=DATA(ndat).experiment;
    
	%% Sort into LFP and EMG data    
    k=1;
    for i=1:nch
        if any(strcmp(CH{i},{'A','C','M','P','L'}))
            CHn(k)=i;
            CHlab{k}=CH{i};
            if isempty(DATA(ndat).interval)
                int(:,k)=minmax(t);
            else
                if isempty(DATA(ndat).interval{i})
                    int(:,k)=minmax(t);
                else
                    int(:,k)=DATA(ndat).interval{i};
                end
            end
            k=k+1;
        end
    end
    nLFP=k-1;
    j=1;
    for i=1:nch
        if any(strncmp(CH{i},{'EDC','FDL','FDI'},3)) && j<=3
            CHn(k)=i;
            CHlab{k}=CH{i};
            if isempty(DATA(ndat).interval)
                int(:,k)=minmax(t);
            else
                if isempty(DATA(ndat).interval{i})
                    int(:,k)=minmax(t);
                else
                    int(:,k)=DATA(ndat).interval{i};
                end
            end
            k=k+1;
            j=j+1; %count until 3 EMG channels are found
        end
    end
    nEMG=k-1-nLFP;
    data=data(:,CHn);
    clear k i nch CHn
    
    
    %% Plot time series
    fig1=figure('Papertype', 'A4', ...
        'Paperunits', 'normalized', 'PaperPosition', [0 0 1 1], ...
        'PaperPositionmode', 'manual', 'Visible', 'off', ...
        'PaperOrientation', 'portrait');
    tslim=repmat([-30 30]',1,nLFP);
    tslim=[tslim repmat([-20 20]',1,nEMG)];
    plot_LFP_EMG(t,data(:,1:nLFP),data(:,nLFP+1:nLFP+nEMG), 'CH', CHlab,...
        'nosca', 1, 'interval', int, 'tslim', tslim);

    %% Plot scalograms
    fig2=figure('Papertype', 'A4', ...
        'Paperunits', 'normalized', 'PaperPosition', [0 0 1 1], ...
        'PaperPositionmode', 'manual', 'Visible', 'off', ...
        'PaperOrientation', 'portrait');
    [LFPf, WLFP, coi1]=procdata(data(:,1:nLFP), 'freq', f);
    [EMGf, WEMG, coi2]=procdata(data(:,nLFP+1:nLFP+nEMG), 'filter', [0 60; 90 110],...
        'rect', 1, 'freq', f);
    WLFP=2*dt*abs(WLFP).^2;
    WEMG=2*dt*abs(WEMG).^2;
    limL=[-9 -4]; %in log10(V^2/Hz)
    limE=[-11 -6];
    
    for i=1:nLFP;
        subplot(max([nLFP,nEMG]),2,2*i-1)
        plot_sca(t, f, WLFP(:,:,i), 'coi', coi, 'lim', limL, ...
            'interval', int(:,i), 'shading', 'flat');
        xlim( minmax( t(:)' ) )
        ylim( minmax( f(:)' ) )
        ylabel('f [Hz]', 'Fontsize', 8)
        title(['LFP ' CHlab{i}], 'Fontsize', 8)
    end
    
    %% Plot EMG
    for i=1:nEMG;
        subplot(max([nLFP,nEMG]),2,2*i)
        plot_sca(t, f, WEMG(:,:,i), 'coi', coi, 'lim', limE, ...
            'interval', int(:,nLFP+i), 'shading', 'flat');
        xlim( minmax( t(:)' ) )
        ylim( minmax( f(:)' ) )
        ylabel('f [Hz]', 'Fontsize', 8)
        title(CHlab{nLFP+i}, 'Fontsize', 8)
    end
    
    %% Plot annotations and save
    annotation(fig1,'textbox', [0.25 0.95 0.5 0.02],...
        'String', ['Name: ' patient num2str(site)...
            ', Experiment: ' exp ', rigor: ', num2str(rigor)],...
        'HorizontalAlign', 'center');
    annotation(fig2,'textbox', [0.25 0.95 0.5 0.02],...
        'String', ['Name: ' patient num2str(site)...
            ', LFPact: ', num2str(LFPact), ', Apo: ', num2str(apo)],...
        'HorizontalAlign', 'center');
    set(fig2,'renderer','painters');
    saveas(fig1, [whereto 'test-' num2str(2*ind-1,'%03u') '.pdf']);
    saveas(fig2, [whereto 'test-' num2str(2*ind,'%03u') '.pdf']);
    close(fig1);
    close(fig2);
    ind=ind+1;
    clear CHlab int tslim
end