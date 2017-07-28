%% Plot comparison of true/false pos. rate between PCC, IC and PLI


% Parameters
bin1 = 0:0.01:1;
bin2 = 0:0.002:0.2;

% Imaginary part of Coherency
load PSD_synth_IC_p01_ns6_w12_nl3.mat
Htpr_IC = hist(TPR, bin1);
Hfpr_IC = hist(FPR, bin2);

% Imaginary part of Coherency
load PSD_synth_wPLI_p01_ns6_w12_nl3.mat
Htpr_PLI = hist(TPR, bin1);
Hfpr_PLI = hist(FPR, bin2);

% PCC
load PSD_synth_PCC_ns6_w12_nl3_v2.mat
Htpr_PCC = hist(TPR, bin1);
Hfpr_PCC = hist(FPR, bin2);

bin1 = bin1*100;
bin2 = bin2*100;

% Plot result
subplot(1,2,1)
plot( bin1, Htpr_IC,  '-m', 'DisplayName', 'TPR(IC)' )
hold on
plot( bin1, Htpr_PLI, '-c', 'DisplayName', 'TPR(PLI)' )
plot( bin1, Htpr_PCC, '-k', 'DisplayName', 'TPR(PCC)' )
legend('show', 'Location', 'NW')
xlabel('TPR=TP/P [%]')
ylabel('# trials')
title('True positive rate (sensitivity)')
xlim([0 80])

subplot(1,2,2)
plot( bin2, Hfpr_IC,  '-m', 'DisplayName', 'FPR(IC)' )
hold on
plot( bin2, Hfpr_PLI, '-c', 'DisplayName', 'FPR(PLI)' )
plot( bin2, Hfpr_PCC, '-k', 'DisplayName', 'FPR(PCC)' )
legend('show')
xlabel('FPR=FP/N [%]')
ylabel('# trials')
title('False positive rate')
xlim([0 20])  