%% Plots power composition for bipolar

load pow_coh_bip-cohinc_rest
tit = 'Rest, bipolar rec.,';

Rcoh_off = P_coh_off(:,:) ./ P_bip_off(:,:);
Rinc_off = P_inc_off(:,:) ./ P_bip_off(:,:);
Rcoh_on  = P_coh_on(:,:)  ./ P_bip_on(:,:);
Rinc_on  = P_inc_on(:,:)  ./ P_bip_on(:,:);

subplot(1,2,1)
h1=shadedErrorBar(f, nanmean(Rinc_off,2)*100, nanstd(Rinc_off,1,2)*100, '-r', 1);
hold all
h2=shadedErrorBar(f, nanmean(Rcoh_off,2)*100, nanstd(Rcoh_off,1,2)*100, '-b', 1);
ylim([0 100]), xlim([1 85])
ylabel('P_i/P_{tot} [%]')
xlabel('f [Hz]')
title([tit ' OFF'])

subplot(1,2,2)
h1=shadedErrorBar(f, nanmean(Rinc_on,2)*100, nanstd(Rinc_on,1,2)*100, '-r', 1);
hold all
h2=shadedErrorBar(f, nanmean(Rcoh_on,2)*100, nanstd(Rcoh_on,1,2)*100, '-b', 1);
ylim([0 100]), xlim([1 85])
ylabel('P_i/P_{tot} [%]')
xlabel('f [Hz]')
title([tit ' ON'])
legend([h1.patch h2.patch], 'Incoherent', 'Coherent')