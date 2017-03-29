%% Plot wavelet coherence scalogram

function plot_coh(t,f,C,coi,nintp)

if nargin<5; nintp=100; end
if nargin<4; coi=gencoi(length(t)); end

if nintp~=0
    ti = linspace(min(t),max(t),nintp);
    Ci=NaN(length(f),length(ti));
    for i=1:length(f)
        Ci(i,:) = interp1(t,C(i,:),ti);
    end
    coii=interp1(t,coi,ti);
    coi=coii;
    C=Ci; t=ti;
    clear Ci ti coii
end

pcolor(t,f,C);
set(gca,'YDir','Reverse');
shading interp
caxis([0 1])
hold all
plot(t,1./coi,'k','linew',1);
hold off

%% Set axis
set(gca,'TickDir','out')
xlim([minmax(t)])
ylim([minmax(f)])
xlabel('Time [s]')
ylabel('Frequency [Hz]')