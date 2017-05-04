
%% Set 1% significance threshold for coherence in time-freqeuncy domain
%% according to sig_test_coh_table.m

function sig = sig_coh_thresh ( w0, nsig )

load('threshold_1perc_w0_ns_table.mat', 'thres')

sig = thres(ceil(w0/2), nsig);