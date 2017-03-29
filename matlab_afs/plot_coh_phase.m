function h=plot_coh_phase( t, f, C, Wxy, varargin )
%% PLOT_COH_PHASE Plot wavelet coherence scalogram with phase arrows
% 
%   PLOT_COH_PHASE( t, f, C, Wxy, varargin )
% 
%   INPUT:
%           t:          Time vector (nt)
%           f:          Frequency vector (nf)
%           C:          Coherence Matrix (nf x nt)
%           Wxy:        Cross wavelet coefficient matrix (nf x nt)
% 
%   OPTIONAL INPUT:
%           coi:        Cone of influence vector (nt)
%           nintp:      Interpolate coefficients to smaller grid to improve 
%                       visualization speed (nintp=500)
%           sig:        ?????????TODO?????
%           dt:         Time increment
%           aheadsize:  ?????????TODO?????
%           asize:      ?????????TODO?????
%           adensity:   ?????????TODO?????
%           shading:    ?????????TODO?????
%           title:      Define plot title
%           clim:       Colormap limit
% 
% Author: Michael von Papen
% 
% Date: 14.10.15

args = struct('coi',-1,...
            'nintp', 500,...
            'sig',0.495,...  %p05 = [0.86 0.645 0.495] for n_sigma = [1 2 3]
            'dt',1/2500,...
            'aheadsize', 2,...
            'asize', 1/30,...
            'adensity',[40 40],...
            'shading', 'flat',...
            'titel', [],...
            'clim', [0 1]);
        
args = parseArgs(varargin,args);

dt = 1/2500;
coi = args.coi;
nintp = args.nintp;

if nintp ~= 0
    ti = linspace(min(t),max(t),nintp);
    Wi = NaN(length(f),length(ti));
    for i = 1:length(f)
        Wi(i,:) = interp1(t,Wxy(i,:),ti);
        Ci(i,:) = interp1(t,C(i,:),ti);
    end
    if coi ~= -1
        coii = interp1(t,coi,ti);
        coi = coii;
    end
    Wxy = Wi;
    t = ti;
    C = Ci;
    dt = (t(end)-t(1))/nintp;
    clear Wi Ci ti coii
end

if coi == -1
    coi = gencoi(length(t),dt);
end


ArrowDensity = args.adensity;
ArrowHeadSize = args.aheadsize;
ArrowSize = args.asize;
sig = args.sig;

h=pcolor(t,log10(f),C);
caxis(args.clim)
shading(gca, args.shading)

set(gca,'YLim',log10(minmax(f)), ...
    'YDir','rev');%, 'layer','top', ...
%     'YTick',log2(Yticks(:)), ...
%     'YTickLabel',num2str(Yticks'), ...
%     'layer','top')
ylabel('log_1_0(f [Hz]')
hold all
plot(t,log10(1./coi),'k','linew',2);
title(args.titel)

%phase plot
aWxy = angle(Wxy);
aWxy(abs(C) < sig) = NaN;

phs_dt = round(length(t)/ArrowDensity(1));
tidx = max(floor(phs_dt/2),1):phs_dt:length(t);
phs_dp = round(length(f)/ArrowDensity(2));
pidx = max(floor(phs_dp/2),1):phs_dp:length(f);
phaseplot(t(tidx),log10(f(pidx)),aWxy(pidx,tidx),...
    ArrowSize,ArrowHeadSize);
%% arrows point in following directions:
% 90° up, 180° left, 270° down, 0° right
% 90° means signal 1 is 90° before signal 2, with Wxy(:,:,sig1,sig2)


% tt = [t([1 1])'-dt*.5;t';t([end end])'+dt*.5];
% hcoi = fill(tt,log10([f([1 end]) 1./coi f([end 1])]),'w');
% set(hcoi,'alphadatamapping','direct','facealpha',.5)
% hold off

%% Set axis
set(gca,'TickDir','out')
xlim(minmax(t))
ylim(minmax(log10(f)))
xlabel('Time [s]')
ylabel('f [Hz]')