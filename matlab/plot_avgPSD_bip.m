% Plots the average PSD of bipolar recordings in an mseb plot

function h = plot_avgPSD_bip( task, fnam, offon )

if nargin<2
    fnam = 'LFP_pow_PCC_v1_filt_';
end
     
legstr = {'$P_\mathrm{bip}$(OFF)', '$P_\mathrm{bip}$(ON)'};


%% MSEB PLOT of bipolar recordings
load([fnam task '.mat'])

N = [sum(~isnan(P_bip_off(:,:)),2) sum(~isnan(P_bip_on(:,:)),2)];
M = [nanmean(P_bip_off(:,:),2) nanmean(P_bip_on(:,:),2)]'+1e-10;
E = [ nanstd(P_bip_off(:,:),0,2) nanstd(P_bip_on(:,:),0,2)]'./sqrt(N');
E = min(cat(3, E, M),[],3);
lineprops.col={'k'; [0.5 0.5 0.5]};

h = mseb(f, M(offon,:), E(offon,:), lineprops, 1);
set(gca, 'YSca', 'log')
xlabel('f [Hz]', 'Interp', 'Latex')
ylabel('PSD [V$^2$/Hz]', 'Interp', 'Latex')

h1 = legend({legstr{offon}});
set(h1, 'Interp', 'Latex')
xlim([1 45])
ylim([1e-9 4e-6])
title(['bipolar rec. during ' task])