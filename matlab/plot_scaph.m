%% Plot wavelet scalogram

function h=plot_scaph(t,period,Wp,varargin)

args=struct('lim',[0 360],...
            'coi',-1,...
            'linear',0,...
            'nintp', 1000,...
            'phasecor',0);
args=parseArgs(varargin,args);

linear=args.linear;
coi=args.coi;
lim=args.lim;
nintp=args.nintp;

%% Interpolate data to smaller set
if nintp~=0
    ti = linspace(min(t),max(t),nintp);
    for i=1:length(period)
        Wpi(i,:) = interp1(t,Wp(i,:),ti);
    end
    if coi~=1
        coii=interp1(t,coi,ti);
        coi=coii;
    end
    Wp=Wpi; t=ti;
    clear Wi ti coii
end

%% Phase correction: leads to nonchanging 'phase' for periodic signal
if args.phasecor==1
    PhC=Wp;
    for i=1:length(period)
        PhC(i,:)=mod(t./period(i)*360,360);
    end
    Wp=mod(Wp-PhC,360);
end

%% Plot phase scalogram
h=pcolor(t,1./period,Wp);
set(gca,'YDir','Reverse');
if linear == 0
    set(gca,'YScale','log');
end
shading interp
caxis(lim)
if coi ~= -1
    hold all
    plot(t,1./coi,'k','linew',2);
end
hold off

%% Set axis
set(gca,'TickDir','out')
xlim(minmax(t))
ylim([min(1./coi) max(1./period)])
xlabel('Time [s]')
ylabel('Frequency [Hz]')