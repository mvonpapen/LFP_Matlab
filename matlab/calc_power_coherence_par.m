%% Test mono- and bipolar PSD changes and coherence from OFF to ON


mypool = parpool(4);

TASK = {'fist', 'hold', 'rest'};

parfor ind = 1:3
    task = TASK{ind};
    calc_pow_v2(task);
%     calc_power_coherence(task);
%     calc_pow_coh_4_bip(task);

end %type for loop end

delete(mypool)