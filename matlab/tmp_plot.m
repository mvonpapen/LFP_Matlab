i=1;
load pow_coh_ds10_pc15_w6n3_w0_06_rest
subplot(2,2,1)
semilogy(f,[nanmean(P_tot_off(:,:,i),2) nanmean(P_inc_off(:,:,i),2) nanmean(P_coh_off(:,:,i),2) nanmean(P_vc_off(:,:,i),2)])
title('w6n3 OFF'), ylim([1e-9 1e-5]), xlim([1 45])
subplot(2,2,2)
semilogy(f,[nanmean(P_tot_on(:,:,i),2) nanmean(P_inc_on(:,:,i),2) nanmean(P_coh_on(:,:,i),2) nanmean(P_vc_on(:,:,i),2)])
title('w6n3 ON'), ylim([1e-9 1e-5]), xlim([1 45])
load pow_coh_ds10_pc15_w0_12_rest
subplot(2,2,3)
semilogy(f,[nanmean(P_tot_off(:,:,i),2) nanmean(P_inc_off(:,:,i),2) nanmean(P_coh_off(:,:,i),2) nanmean(P_vc_off(:,:,i),2)])
title('w12 OFF'), ylim([1e-9 1e-5]), xlim([1 45])
subplot(2,2,4)
semilogy(f,[nanmean(P_tot_on(:,:,i),2) nanmean(P_inc_on(:,:,i),2) nanmean(P_coh_on(:,:,i),2) nanmean(P_vc_on(:,:,i),2)]), legend('tot', 'inc', 'coh', 'vc')
title('w12 ON'), ylim([1e-9 1e-5]), xlim([1 45])