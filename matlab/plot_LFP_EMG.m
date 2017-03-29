%%  Plot EMG and LFP data after load
%%
%%
%%  Input:
%%          t           time vector, n elements
%%          LFP         LFP data, n x m1 elements
%%          EMG         EMG data, n x m2 elements
%%
%%  Optional:
%%          LFP2        second set of LFP data, n x m1 elements
%%          EMG2        second set of EMG data, n x m2 elements
%%          SCG         Wave power of one time series to plot as scalogram, s x n elements
%%          period      needed to plot scalogram, s elemnts
%%          lim         color axis limits for scalogram
%%          coi         cone of influence corresponding to W, n elements
%%          linear      linear=1 plots y-axis of scalogram linear
%%          nosca       nosca=1 if no scalogram is desired
%%
%% Author:  M. von Papen
%% Date:    28.10.2014

function [h, hs]=plot_LFP_EMG(t,LFP,EMG,varargin)

args=struct('LFP2',[],...
            'EMG2',[],...
            'SCG',[],...
            'period',[],...
            'lim',[-5 2],...
            'tslim',[],...
            'coi',-1,...
            'linear',0,...
            'nosca',0,...
            'interval',[],...
            'CH',[]);
args=parseArgs(varargin,args);


if args.nosca==1;
    scg=0;
else
    scg=1; %set counter one up if necessary to make room for scalogram
end

%% Parameters
nLFP=size(LFP,2);
nEMG=size(EMG,2);

%% Labels
ylab1='LFP [mV]';
ylab2='EMG [mV]';
xlab ='Time [s]';


%% Plot LFP
lim=[min(min(LFP*1e3)) max(max(LFP*1e3))];
for i=1:nLFP;
    subplot(max([nLFP,nEMG+scg]),2,2*i-1)
    h(i)=plot(t, LFP(:,i)*1e3,'black');
    hold all
    if ~isempty(args.LFP2)
        plot(t,args.LFP2(:,i)*1e3,'--')
    end
    if ~isempty(args.interval)
        plot(args.interval(1,i)*[1 1],10*lim,'--k')
        plot(args.interval(2,i)*[1 1],10*lim,'--k')
    end
    if i==nLFP
        xlabel(xlab);
    end
    ylabel(ylab1);
    xlim( minmax( t ) );
    if ~isempty(args.tslim)
        ylim( args.tslim(:,i) );
    end
    if isempty(args.CH)
        if i==1
            title('LFP Channels');
        end
    else
        title(args.CH{i})
    end
end

%% Plot EMG
lim=[min(min(EMG*1e3)) max(max(EMG*1e3))];
for i=1:nEMG;
    subplot(max([nLFP,nEMG+scg]),2,2*i);
    h(nLFP+i)=plot(t, EMG(:,i)*1e3,'black');
    hold all
    if ~isempty(args.EMG2)
        plot(t,args.EMG2(:,i)*1e3,'--')
    end
    if ~isempty(args.interval)
        plot(args.interval(1,nLFP+i)*[1 1],10*lim,'--k')
        plot(args.interval(2,nLFP+i)*[1 1],10*lim,'--k')
    end
    if i==nEMG
        xlabel(xlab);
    end
    ylabel(ylab2);
    xlim( minmax( t ) );
    if ~isempty(args.tslim)
        ylim( args.tslim(:,nLFP+i) );
    end
    if isempty(args.CH)
        if i==1
            title('EMG Channels');
        end
    else
        title(args.CH{i+nLFP})
    end
end


if scg==1
    hs=subplot(max([nLFP,nEMG+1]),2,2*(nEMG+1));
end
if ~isempty(args.SCG)
    plot_sca(t, args.period, args.SCG, 'lim',args.lim,...
        'coi',args.coi, 'linear',args.linear);
end