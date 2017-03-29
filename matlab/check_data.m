%%  Plot and analyze a single channel
%%
%%
%%  Input:
%%          data        data with n elements
%%          T           time window to analyze
%%
%% Author:  M. von Papen
%% Date:    05.12.2014

%% Set which data to analyze
type='L';
ch=1;
% T=[0 60];


%% Define data vector and set parameters
dt=1/2500;
t=[0:length(data)-1]*dt;
i=find(t>=T(1) & t<=T(2));
t=t(i);
dat=data(i,ch)*1e3; % in mV umgewandelt
f=logspace(-1,2.5,50);


%% Filtered data set and wavelet transform
switch type
    case 'L'
        argstr=['''filter'',[], ''order'', 5,'...
            '''freq'', f, ''inc'', 1'];
        [datf, W, coi, P]=procdata(dat,'filter',[],...
            'order', 5, 'freq', f, 'inc', 1);
        flim=[0.1 300];
        lim=[-4 0];
    case 'E'
        [datf, W, coi, P]=procdata(dat,'filter',[0 60; 90 110], 'rect', 1, 'freq', f);
        argstr='''filter'',[0 60; 90 110], ''freq'', f, ''rect'', 1';
        flim=[0.1 300];
        lim=[-4 0];
end
    


%% Plot time series and scalogram

%% Labels
ylab ='mV';
yplab ='PSD [mV^2/Hz]';
tlab ='Time [s]';
flab ='f [Hz]';


%% Plot
ha(1)=subplot(2,2,1);
h(1)=plot(t, dat,'black');
xlabel(tlab);
ylabel(ylab);
t_str=mat2str(minmax(t));
title(['minmax(t)=' t_str]);
subplot(2,2,3)
h(3)=plot(t, datf,'black');
xlabel(tlab);
ylabel(ylab);
title('Processed data')


%% Allow interactive editing of data
%% Brushed data is set to NaN after deletion
for i=[1 3]
    set(h(i),'XDataSOURCE','t');
end
set(h(1),'YDataSOURCE','dat');
set(h(3),'YDataSOURCE','datf');
brush on;
linkdata on;

%% Plot data that shall not be linked
subplot(2,2,2);
loglog(f, P,'black');
xlabel(flab);
ylabel(yplab);
xlim(flim);
hold all
ha(4)=subplot(2,2,4);
plot_sca(t, 1./f, 2*dt*abs(W).^2, 'lim',lim,'coi',coi);


%% Insert user interface controls
he=0.02; wi=0.04;
uicontrol('Style', 'pushbutton',...
    'String', 'Refresh',...
    'Units', 'normalized', 'Position', [0.05 0.01 wi he],...
    'Callback', ['i=~isnan(t); t=t(i); dat=dat(i); clear i;'...
        '[datf, W, coi, P]=procdata(dat,' argstr ');'...
        'refreshdata; subplot(2,2,2); plot(f,P); subplot(ha(4));'...
        'plot_sca(t, 1./f, 2*dt*abs(W).^2, ''lim'',lim, ''coi'',coi);'...
        't_str=mat2str(minmax(t)); title(ha(1), [''minmax(t)='' t_str])']);
uicontrol('Style', 'pushbutton',...
    'String', 'Original',...
    'Units', 'normalized', 'Position', [0.1 0.01 wi he],...
    'Callback', ['t=[0:length(data)-1]*dt; i=find(t>=T(1) & t<=T(2));'...
        't=t(i); dat=data(i,ch);']);
uicontrol('Style', 'pushbutton',...
    'String', 'C-',...
    'Units', 'normalized', 'Position', [0.92 0.13 0.02 0.05],...
    'Callback', 'lim=lim-0.5; caxis(ha(4),lim);');
uicontrol('Style', 'pushbutton',...
    'String', 'C+',...
    'Units', 'normalized', 'Position', [0.92 0.18 0.02 0.05],...
    'Callback', 'lim=lim+0.5; caxis(ha(4),lim);');

%% Force figure toolbar
set( gcf, 'toolbar', 'figure' )