%%  Plot EMG and LFP data after load
%%
%%
%%  Input:
%%          data        data matrix n x m with number of channels m=nLFP+nEMG
%%                      and number of data points n
%%          nLFP        number of LFP channels
%%          nEMG        number of EMG channels
%%          T           time window to analyze
%%
%% Author:  M. von Papen
%% Date:    28.10.2014

% function check_LFP_EMG(data,nLFP,nEMG,T)

% %% Ckeck input
% if nargin<4; T=[0 9e9]; end
% if nargin<3; nEMG=3; end
% if nargin<2; nLFP=5; end




%% Set data vector
dt=1/2456;
lim=[-7 -3];
t=[1:length(data)-1]*dt;
data=data(2:end,:);
i=find(t>=T(1) & t<=T(2));
t=t(i);
LFP=data(i,1:nLFP); %[1 4 5]); %
EMG=data(i,nLFP+1:nLFP+nEMG); %6:8); %

% %% Filtered data set
LFPf=procdata(LFP,'filter',[45 55; 90 110]);
EMGf=procdata(EMG,'filter',[0 60; 90 110], 'rect', 1);

%% Wavelet transform
f=logspace(-1,2,50);
% W=abs(WLFP(:,:,1)).^2;

%% Plot time series and scalogram
l=LFP;
e=EMG;
[h, hs]=plot_LFP_EMG(t,l,e);


%% Allow interactive editing of data
%% Brushed data is set to NaN after deletion
for i=1:nLFP
    set(h(i),'XDataSOURCE','t');
    set(h(i),'YDataSOURCE',['l(:,' mat2str(i) ')']);
end
for i=1:nEMG
    set(h(nLFP+i),'XDataSOURCE','t');
    set(h(nLFP+i),'YDataSOURCE',['e(:,' mat2str(i) ')']);
end
linkdata on;
brush on;


%% Insert user interface controls
he=0.02; wi=0.04;
uicontrol('Style', 'pushbutton',...
    'String', 'Refresh',...
    'Units', 'normalized', 'Position', [0.05 0.01 wi he],...
    'Callback', ['i=~isnan(t); t=t(i); LFP=LFP(i,:); EMG=EMG(i,:); clear i;'...
        '[LFPf, WLFP, coi]=procdata(LFP, ''freq'', f);'...
        '[EMGf, WEMG]=procdata(EMG,''filter'',[0 60;90 110], ''rect'', 1, ''freq'', f);'...
        'W=2*dt*abs(WLFP(:,:,1)).^2; l=LFP; e=EMG; refreshdata']);
%     '[h, hs]=plot_LFP_EMG(t,LFP,EMG);'
uicontrol('Style', 'pushbutton',...
    'String', 'Original',...
    'Units', 'normalized', 'Position', [0.15 0.01 wi he],...
    'Callback', 'l=LFP; e=EMG; refreshdata');
uicontrol('Style', 'pushbutton',...
    'String', 'Filtered',...
    'Units', 'normalized', 'Position', [0.25 0.01 wi he],...
    'Callback', 'l=LFPf; e=EMGf; refreshdata');

%% Buttons for scalogram plot
for j=1:nLFP
    uicontrol('Style', 'pushbutton',...
        'String', ['LFP' mat2str(j)],...
        'Units', 'normalized', 'Position', [0.45+j*0.05 0.02+he wi he],...
        'Callback', ['W=2*dt*abs(WLFP(:,:,' mat2str(j) ')).^2; subplot(hs);'...
            'plot_sca(t, f, W, ''lim'',lim, ''coi'',coi,'...
            '''linear'',0);']);
end
for j=1:nEMG
    uicontrol('Style', 'pushbutton',...
        'String', ['EMG' mat2str(j)],...
        'Units', 'normalized', 'Position', [0.45+j*0.05 0.01 wi he],...
        'Callback', ['W=abs(WEMG(:,:,' mat2str(j) ')).^2; subplot(hs);'...
            'plot_sca(t, f, W, ''lim'',lim, ''coi'',coi,'...
            '''linear'',0);']);
end
uicontrol('Style', 'pushbutton',...
    'String', 'C-',...
    'Units', 'normalized', 'Position', [0.92 0.13 0.02 0.05],...
    'Callback', 'lim=lim-0.5; caxis(hs,lim);');
uicontrol('Style', 'pushbutton',...
    'String', 'C+',...
    'Units', 'normalized', 'Position', [0.92 0.18 0.02 0.05],...
    'Callback', 'lim=lim+0.5; caxis(hs,lim);');

%% Force figure toolbar
set( gcf, 'toolbar', 'figure' )

% clear EMG EMGf LFP LFPf e l