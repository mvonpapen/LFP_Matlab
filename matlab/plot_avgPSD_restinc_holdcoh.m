% Plots two subplots. 1st: incoherent signal during rest OFF/ON, 
%                     2nd: coherent signal during hold OFF/ON


task = 'hold';
fnam   = ['avg_psd_' task '_coh_v1.eps'];
fnamPCC = 'LFP_pow_PCC_v1_filt_';
fnamIC  = 'LFP_pow_IC_p01_v1_filt_';
fnamPLI = 'LFP_pow_wPLI_p01_v1_filt_';
load([fnamPCC task])

%% Set variables and parameters
Mtot = zeros(length(f), 2);
Mpcc = zeros(length(f), 2);
Mic  = zeros(length(f), 2);
Mpli = zeros(length(f), 2);


fig = figure('Papersize', [7 6], 'PaperPosition', [0.75 0.5 12.5 5], ...
   'PaperPositionmode', 'manual', 'Visible', 'on');

%% Load data and average
load([fnamPCC task '.mat']);
Mtot(:,:) = [nanmean(P_tot_off(:,:),2) nanmean(P_tot_on(:,:),2)  ];
Mpcc(:,:)  = [nanmean(P_coh_off(:,:),2) nanmean(P_coh_on(:,:),2)  ];
load([fnamIC task '.mat'])
Mic(:,:)   = [nanmean(P_coh_off(:,:),2) nanmean(P_coh_on(:,:),2)  ];
load([fnamPLI task '.mat'])
Mpli(:,:) = [nanmean(P_coh_off(:,:),2) nanmean(P_coh_on(:,:),2)  ];
     
        
%% PLOT
h = semilogy(f, Mtot);
hold on
set(h(1), 'color', 'k')
set(h(2), 'color', [0.5 0.5 0.5])
h = semilogy(f, Mpcc);
set(h(1), 'color', 'b')
set(h(2), 'color', [0.5 0.5 1])
h = semilogy(f, Mic);
set(h(1), 'color', 'b', 'linestyle', '--')
set(h(2), 'color', [0.5 0.5 1], 'linestyle', '--')
h = semilogy(f, Mpli);
set(h(1), 'color', 'b', 'linestyle', '-.')
set(h(2), 'color', [0.5 0.5 1], 'linestyle', '-.')
xlabel('f [Hz]', 'Interp', 'Latex')
ylabel('PSD [V$^2$/Hz]', 'Interp', 'Latex')

h = legend('$P_\mathrm{tot}$(OFF)', '$P_\mathrm{tot}$(ON)', ...
           '$P_\mathrm{coh,pcc}$(OFF)', '$P_\mathrm{coh,pcc}$(ON)', ...
           '$P_\mathrm{coh,ic}$(OFF)', '$P_\mathrm{coh,ic}$(ON)', ...
           '$P_\mathrm{coh,pli}$(OFF)', '$P_\mathrm{coh,pli}$(ON)');
set(h, 'Interp', 'Latex')
xlim([1 45])
ylim([1e-9 4e-6])
title(['Average PSD during ' task])

% % Estimate significance with signrank
% PcohON = P_coh_on(:,:); PcohOFF = P_coh_off(:,:);
% for i=1:length(f); p(i) = signrank(PcohOFF(i,:)-PcohON(i,:)); end
% f(p<0.05)
% % PincON = P_inc_on(:,:); PincOFF = P_inc_off(:,:);
% % for i=1:length(f); p(i) = signrank(PincOFF(i,:)-PincON(i,:)); end
% % f(p<0.05)
% hold all
% plot([13 21], [1 1]*5e-7, 'r-', 'linewidth', 2)
% plot(17, 6e-7, 'r*', 'markersize', 5)