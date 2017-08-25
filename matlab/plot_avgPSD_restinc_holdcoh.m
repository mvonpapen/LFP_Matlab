% Plots two subplots. 1st: incoherent signal during rest OFF/ON, 
%                     2nd: coherent signal during hold OFF/ON

fignam  = ['avg_psd_restinc_holdcoh_v2'];
fnamPCC = 'LFP_pow_PCC_v1_filt_';


fig = figure('Papersize', [20 9], 'PaperPosition', [0.75 0.5 18.5 8], ...
   'PaperPositionmode', 'manual', 'Visible', 'off');

%% Load data and average
load([fnamPCC 'rest.mat']);
Mtot_rest = [nanmean(P_tot_off(:,:),2) nanmean(P_tot_on(:,:),2)  ]';
N         = [sum(~isnan(P_tot_off(:,:)),2) sum(~isnan(P_tot_on(:,:)),2)];
Etot_rest = [ nanstd(P_tot_off(:,:),0,2) nanstd(P_tot_on(:,:),0,2) ]'./sqrt(N');
Minc_rest  = [nanmean(P_inc_off(:,:),2) nanmean(P_inc_on(:,:),2)  ]';
N         = [sum(~isnan(P_inc_off(:,:)),2) sum(~isnan(P_inc_on(:,:)),2)];
Einc_rest = [ nanstd(P_inc_off(:,:),0,2) nanstd(P_inc_on(:,:),0,2) ]'./sqrt(N');
load([fnamPCC 'hold.mat'])
Mtot_hold = [nanmean(P_tot_off(:,:),2) nanmean(P_tot_on(:,:),2)  ]';
N         = [sum(~isnan(P_tot_off(:,:)),2) sum(~isnan(P_tot_on(:,:)),2)];
Etot_hold = [ nanstd(P_tot_off(:,:),0,2) nanstd(P_tot_on(:,:),0,2) ]'./sqrt(N');
Mcoh_hold   = [nanmean(P_coh_off(:,:),2) nanmean(P_coh_on(:,:),2)  ]';
N         = [sum(~isnan(P_coh_off(:,:)),2) sum(~isnan(P_coh_on(:,:)),2)];
Ecoh_hold = [ nanstd(P_coh_off(:,:),0,2) nanstd(P_coh_on(:,:),0,2) ]'./sqrt(N');
     
        
%% MSEB PLOTS
clear lineprops
subplot(1,2,1)
lineprops.col={'k'; [0.7 0.7 0.7]};
% lineprops.edgestyle='-';
mseb(f, Mtot_rest, Etot_rest, lineprops);
hold on

lineprops.col={'r'; [1 0.7 0.7]};
Einc_rest = min(cat(3, Einc_rest, Minc_rest-(1e-10)),[],3);
Einc_rest(Minc_rest-(1e-10)-Einc_rest<=0) = NaN;
% lineprops.edgestyle='-';
mseb(f, Minc_rest, Einc_rest, lineprops, 0);
set(gca, 'YSca', 'log')
xlabel('f [Hz]', 'Interp', 'Latex')
ylabel('PSD [V$^2$/Hz]', 'Interp', 'Latex')
h = legend('$P_\mathrm{tot}$(OFF)', '$P_\mathrm{tot}$(ON)', ...
           '$P_\mathrm{inc}$(OFF)', '$P_\mathrm{inc}$(ON)');
set(h, 'Interp', 'Latex')
xlim([1 45])
ylim([1e-9 4e-6])
title('Average PSD during rest', 'Interp', 'Latex')


clear lineprops
subplot(1,2,2)
lineprops.col={'k'; [0.7 0.7 0.7]};
% lineprops.edgestyle='-';
mseb(f, Mtot_hold, Etot_hold, lineprops, 0);
hold on

lineprops.col={'b'; [0.7 0.7 1]};
Einc_rest = min(cat(3, Ecoh_hold, Mcoh_hold-(1e-10)),[],3);
Einc_rest(Minc_rest-(1e-10)-Einc_rest<=0) = NaN;
% lineprops.edgestyle='-';
mseb(f, Mcoh_hold, Ecoh_hold, lineprops, 0);
set(gca, 'YSca', 'log')
xlabel('f [Hz]', 'Interp', 'Latex')
ylabel('PSD [V$^2$/Hz]', 'Interp', 'Latex')
h = legend('$P_\mathrm{tot}$(OFF)', '$P_\mathrm{tot}$(ON)', ...
           '$P_\mathrm{coh}$(OFF)', '$P_\mathrm{coh}$(ON)');
set(h, 'Interp', 'Latex')
xlim([1 45])
ylim([1e-9 4e-6])
title('Average PSD during hold', 'Interp', 'Latex')



%% Save figure
print(fig, [fignam '.eps'], '-depsc')
saveas(fig, [fignam, '.fig'])
close(fig)