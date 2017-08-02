% Plots coh of PCC, IC and wPLI

task = 'hold';
fignam  = ['avg_psd_' task '_coh_v1'];
fnam{1} = 'LFP_pow_PCC_v1_filt_';
fnam{2} = 'LFP_pow_IC_p01_v1_filt_';
fnam{3} = 'LFP_pow_wPLI_p01_v1_filt_';
methlabel = {'PCC', 'IC', 'wPLI'};

fig = figure('Papersize', [30 10], 'PaperPosition', [0.75 0.5 28.5 9], ...
   'PaperPositionmode', 'manual', 'Visible', 'on');

for i=1:3
    load([fnam{i} task])
    
    % Significance as a function of f and in high beta band
    Poff = P_coh_off(:,:);
    Pon  = P_coh_on(:,:);
    for k=1:length(f)
        sig_f(k,i) = signrank ( Pon(k,:)-Poff(k,:) );
    end
    Poff = nanmean(Poff(20:30,:));
    Pon  = nanmean(Pon(20:30,:));
    [sig,~,tmp] = signrank ( Pon-Poff );
    fprintf('%s: p=%f, W=%u\n', methlabel{i}, sig, tmp.signedrank)
    
    
    %% MSEB PLOTS
    subplot(1,3,i)
    N = [sum(~isnan(P_tot_off(:,:)),2) sum(~isnan(P_coh_off(:,:)),2)];
    M = [nanmean(P_tot_off(:,:),2) nanmean(P_coh_off(:,:),2)]'+1e-10;
    E = [ nanstd(P_tot_off(:,:),0,2) nanstd(P_coh_off(:,:),0,2) ]'./sqrt(N');
    E = min(cat(3, E, M-(1e-10)),[],3);
    lineprops.style='-';
    lineprops.col={'k'; 'b'};
    mseb(f, M, E, lineprops, 1);
    hold all

    N = [sum(~isnan(P_tot_on(:,:)),2) sum(~isnan(P_coh_on(:,:)),2)];
    M = [nanmean(P_tot_on(:,:),2) nanmean(P_coh_on(:,:),2)]'+1e-10;
    E = [ nanstd(P_tot_on(:,:),0,2) nanstd(P_coh_on(:,:),0,2) ]'./sqrt(N');
    E = min(cat(3, E, M),[],3);
    lineprops.style='--';
    lineprops.col={[0.5 0.5 0.5]; [0.5 0.5 1]};

    mseb(f, M, E, lineprops, 1);
    set(gca, 'YSca', 'log')
    xlabel('f [Hz]', 'Interp', 'Latex')
    ylabel('PSD [V$^2$/Hz]', 'Interp', 'Latex')

    if i==3
        h = legend('$P_\mathrm{tot}$(OFF)', '$P_\mathrm{coh}$(OFF)', ...
                   '$P_\mathrm{tot}$(ON)',  '$P_\mathrm{coh}$(ON)');
        set(h, 'Interp', 'Latex')
    end
    xlim([1 45])
    ylim([1e-9 2e-6])
    title(['Average PSD for \textbf{' methlabel{i} '} during ' task], ...
        'Interp', 'Latex')
    
    % Annotate significant freq.range
    plot(f(sig_f(:,i)<0.05), 2e-9, '*k')
    
end

print(fig, fignam, '-depsc')
saveas(fig, [fignam '.fig'])