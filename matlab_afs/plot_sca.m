function h = plot_sca(t,f,W,varargin)
%% PLOT_SCA Plot wavelet scalogram
% Using: W = 2*dt*abs(WT).^2
%   h = PLOT_SCA(t,f,W,varargin)
% 
%   INPUT:
%           t:          Time vector (nt)
%           f:          Frequency vector (nf)
%           W:          Wavelet Matrix (nf x nt)
% 
%   OPTIONAL INPUT:
%           lim:        limits for colorbar in log10
%           coi:        Cone of influence vector (nt)
%           linear:     Linear y-axis (frequency) on/off (default off)
%           interval:   vector of time stamps, vertical black line will be drawn
%           nintp:      Interpolate coefficients to smaller grid to improve 
%                       visualization speed (nintp = 500)
%           shading:    define shading for pcolor (default=flat)
%           nolog:      if nolog=1, do not use logarithmic power
%           alphamask:  alphamask matrix (nf x nt) to set opaque levels
% 
% Author: Michael von Papen
% 
% Date: 14.10.15


args = struct('lim',[-9 -4],...
            'coi',-1,...
            'linear',0,...
            'interval',[],...
            'nintp', 200,...
            'shading', 'flat',...
            'nolog',0, ...
            'alphamask', []);
args = parseArgs(varargin,args);

if ~isreal(W)
    W = 2/2456*abs(W).^2;
end

linear = args.linear;
coi = args.coi;
lim = args.lim;
nintp = args.nintp;
C = args.alphamask;

%% Interpolation of coefficients
if nintp ~= 0
    ti = linspace(min(t),max(t),nintp);
    Wi = NaN(length(f),length(ti));
    Ci = NaN(length(f),length(ti));
    for i = 1:length(f)
        Wi(i,:) = interp1(t,W(i,:),ti);
    end
    if coi ~=- 1
        coii = interp1(t,coi,ti);
        coi = coii;
    else
        coi = gencoi(length(ti),(max(ti)-min(ti))/nintp);
    end
    if ~isempty(C)
        for i = 1:length(f)
            Ci(i,:) = interp1(t,C(i,:),ti);
        end
    end
    W = Wi; t = ti; C = Ci;
    clear Wi ti coii
end


%% Plotting
if args.nolog == 0
    h = pcolor(t,f,log10(W));
else
    h = pcolor(t,f,W);
end
if ~isempty(args.alphamask)
    alpha(C)
    alim([0 1])
end
shading(gca, args.shading)
caxis(lim)
colormap('jet')
hold all
plot(t,1./coi,'k','linew',1);
if ~isempty(args.interval)
    plot(args.interval(1)*[1 1],[min(f) max(f)],'--k')
    plot(args.interval(2)*[1 1],[min(f) max(f)],'--k')
end

%% Set axis
set(gca,'YDir','Reverse');
if linear == 0
    set(gca,'YScale','log');
end
set(gca,'TickDir','out')
xlim([min(t) max(t)])
ylim([min(f) max(f)])
xlabel('Time [s]')
ylabel('f [Hz]')