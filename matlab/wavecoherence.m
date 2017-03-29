function varargout=wavecoherence(t,x,y,varargin)
% Possible variable input parameters are: 
% 'pad' (zero padding: 0 or 1, default 0),'dj' (scale step, default 0.2),
% 's0' (smallest scale, default 2*dt),
% 'j1'(maximum scale, or leave blank to use j1=floor(log2(n*dt/s0)/dj)),
% 'mother' ('morlet','paul','dog'),
% 'smooth' ('gauss' or 'box', default gauss),
% 'plot' (0 off, 1 on, default 0),
% 'scaleSmooth' (0 off, 1 on, default 0)
% 'param' (wavelet parameter, defaults are: morlet=6,paul=4,dog=2)
% 
% Possible variable output parameters are:
% [coherence, scales, periods, CrossWaveletCoefficients]

%% input parameters
fSize=14;
n = length(x);
dt = t(2)-t(1) ;

Args=struct('pad',0,'dj',0.2, 's0',2*dt,'j1',[],'mother','Morlet',...
    'smooth','gauss','plot',0,'scaleSmooth',0,'param',-1);
Args=parseArgs(varargin,Args);
if isempty(Args.j1)
        Args.j1=floor(log2(n*dt/Args.s0)/Args.dj); %auto maxscale
end


%% Wavelet transform:
[wavex,period,scale,coix] = wavelet(x,dt,Args.pad,Args.dj,Args.s0,Args.j1,Args.mother);
[wavey,period,scale,coiy] = wavelet(y,dt,Args.pad,Args.dj,Args.s0,Args.j1,Args.mother);
coi=min(coix,coiy);

wavexy=conj(wavex).*wavey;
powx=abs(wavex).^2;
powy=abs(wavey).^2;

%% Wavelet Smoothing
sX=smoothed(powx,dt,scale,Args.dj,'smooth',Args.smooth,'scaleSmooth',Args.scaleSmooth);
sY=smoothed(powy,dt,scale,Args.dj,'smooth',Args.smooth,'scaleSmooth',Args.scaleSmooth);
sXY=smoothed(wavexy,dt,scale,Args.dj,'smooth',Args.smooth,'scaleSmooth',Args.scaleSmooth);

%% calculation of the coherence

coherencesmooth=abs(sXY).^2./(sX.*sY);

%% fourier frequencies

if Args.param==-1
    switch upper(Args.mother)
        case 'MORLET'
            f=1.03./period;        
            coi=1.03./coi;
        case 'PAUL'
            f=1.3963./period;
            coi=1.3963./coi;
        case 'DOG'
            f=3.9738./period;
            coi=3.9738./coi;
    end
    

else
    switch upper(Args.mother)
        case 'MORLET'
            four=4*pi/(Args.param+sqrt(2+Args.param^2));
            f=four./period;        
            coi=four./coi;
        case 'PAUL'
            four=4*pi/(2*Args.param+1);
            f=four./period;        
            coi=four./coi;
        case 'DOG'
            four=2*pi/sqrt(Args.param+0.5);
            f=four./period;        
            coi=four./coi;
    end
end

%% plotting

if Args.plot==1
    titlee=strcat(Args.mother,' wavelet,');
    titlee=strcat(titlee,strcat('scalesmoothing=',int2str(Args.scaleSmooth)));
    titlee=strcat(titlee,')');
     pcolor(t,log(f),coherencesmooth), shading flat;
%     imagesc(t,log(f),coherencesmooth);
    if strcmpi(Args.smooth,'gauss')==1
        titlee=strcat('Coherence (gaussian smoothing, ',titlee);
%         title('coherence with smoothing gaussian ('mother' wavelet','FontSize',fSize);
        title(titlee,'FontSize',fSize);  
    elseif strcmpi(Args.smooth,'box')==1
        titlee=strcat('Coherence (box smoothing, ',titlee);
        title(titlee,'FontSize',fSize);
    end
    hold all;
%     coi=log(coi);
    plot(t,log(coi),'k','LineWidth',3);
%     fill(ones(length(t),1).*min(coi),ones(length(t),1).*(2*min(coi)))
    [mincoi,posmincoi]=min(coi);
    fill([t(1) t(1:posmincoi)],log([mincoi coi(1:posmincoi)]),'w','FaceAlpha',0.5)
    fill([t(posmincoi:length(t)) t(length(t))],log([coi(posmincoi:length(coi)) mincoi]),'w','FaceAlpha',0.5)
    ylabel('f','FontSize',fSize);
    xlabel('t','FontSize',fSize);
    colorbar;
    caxis([0 1]);
    set(gca,'FontSize',fSize);

     ylim([log(min(coi)),log(max(f))]);
%      set(gca,'ytick',fliplr(log(f)))
%      labels=num2cell(f);
%      set(gca,'yticklabel',labels)
    
%          ytick=get(gca,'YTick')';
%          h=gca;
        yticks=[1; 2; 5; 10; 20; 40; 60; 100];
%         h.xtick=[3];
%     ytickmode{manual};
%           set(gca,'yScale','log');
     set(gca,'YTick',log(yticks));
%       set(gca,'YDir','reverse');
%       axis ij;
%     disp(num2str(yticks(2)))
    ticklabel=num2str(yticks);
%     disp(ticklabel)
% 	for i=1:length(yticks)
%         ticklabel(i)=num2str(yticks(i));
%     end
%     disp(ticklabel)
       set(gca,'YTicklabel',ticklabel);


end

%% output
coher=coherencesmooth;
varargout={coher,scale,period,wavexy};
varargout=varargout(1:length(varargout));
end
