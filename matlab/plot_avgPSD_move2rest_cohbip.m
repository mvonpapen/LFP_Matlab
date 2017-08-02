% Plots movement induced changes of coherent signal for PCC

fignam  = ['avg_psd_move2rest_cohbip_v1.eps'];
fnam = 'LFP_pow_PCC_v1_filt_';

fig = figure('Papersize', [30 10], 'PaperPosition', [0.75 0.5 28.5 9], ...
   'PaperPositionmode', 'manual', 'Visible', 'on');


%% Plot
subplot(2,2,1)
plot_avgPSD_PCC('rest', fnam, 3, 1);
plot_avgPSD_PCC('fist', fnam, 3, 1);
subplot(2,2,3)
plot_avgPSD_PCC('rest', fnam, 3, 2);
plot_avgPSD_PCC('fist', fnam, 3, 2);
subplot(2,2,2)
plot_avgPSD_bip('rest', fnam, 1);
plot_avgPSD_bip('fist', fnam, 1);
subplot(2,2,4)
plot_avgPSD_bip('rest', fnam, 2);
plot_avgPSD_bip('fist', fnam, 2);

% print(fig, fignam, '-depsc')