%% Plot phase thres tests

figname = 'synth_pc15_ns6_w12_nl3';
name = 'synth_pc15_ns6_w12_nl3.mat';
% name = 'test_pc15_ns3_w12_P07_1-2';
realdata = 0;

% Load
alpha = 10;
beta  = 20;
gamma = 50;
t     = (0:1e4)/2500;
A     = 500;
load(name)

% Open Figure
fig1 = figure('Papersize', 0.8*[8 5], 'PaperPosition', 0.8*[0.75 0.5 6.5 4], ...
        'PaperPositionmode', 'manual', 'Visible', 'off'); 
    
%% Plot data
% % h1 = shadedErrorBar(f,mean(Ptot,2),min([mean(Ptot,2)'-1e-20; std(Ptot,0,2)']),'k',1);
% % hold all
% % h2 = shadedErrorBar(f,mean(Pinc,2),min([mean(Pinc,2)'-1e-20; std(Pinc,0,2)']),'r',1);
% % h3 = shadedErrorBar(f,mean(Pcoh+Pvc,2),min([mean(Pcoh+Pvc,2)'-1e-20; std(Pcoh,0,2)']),'b',1);
% % h4 = shadedErrorBar(f,mean(Pvc ,2)+1e-19,min([mean(Pvc,2)'-1e-20; std(Pvc,0,2)']),'g',1);
% % x  = sin(2*pi*alpha*t'-dp) + sin(2*pi*beta*t') + sin(2*pi*gamma*t');
% % x  = x / A;
% % [~,~,~,tmp] = procdata(x, 'freq', f, 'w0', w0, 'filt', 0.5);
% % h5 = plot(f,tmp,'k--');




M = [mean(Ptot,2) mean(Pinc,2) mean(Pcoh,2) mean(Pvc,2)+1e-19]';

%% Synthetic data
if realdata~=1
    lineprops.col={'k'; 'r'; 'b'; 'g'};
    E = [ min([mean(Ptot,2)'-1e-20; std(Ptot,0,2)']); ...
          min([mean(Pinc,2)'-1e-20; std(Pinc,0,2)']);...
          min([mean(Pcoh,2)'-1e-20; std(Pcoh,0,2)']);...
          min([mean(Pvc, 2)'-1e-20; std(Pvc, 0,2)']) ]; %% CHECKEN!
    mseb(f, M, E, lineprops);
    hold all
% %     x  = sin(2*pi*alpha*t'-dp) + sin(2*pi*beta*t') + sin(2*pi*gamma*t');
%     xa = [zeros(nt/5,1); sin(2*pi*alpha*t(nt/5+1:4*nt/5)'); zeros(nt/5,1)];
%     xb = sin(2*pi*beta*(1+t'/(nt*dt)*0.25).*t');
%     xg = sin(2*pi*gamma*t')/sqrt(2);
%     Pn = 3 * powlawnoise(length(t),1)';
%     x  = xa + xb + xg + Pn;
%     x  = x / A;
    x = synth_ts_PM_paper( nt, nl, A, dp );
    [~,~,~,tmp] = procdata(x, 'freq', f, 'w0', w0, 'filt', []);
    loglog(f,tmp,'k--');
    h = legend('total (sine+noise)', 'incoherent', 'coherent', ...
        'vol.cond.', 'sine w/o noise');
    set(h, 'Interpreter', 'Latex');
    ylim([1e-8 1e-5])
end

%% Real data
if realdata==1
    set(gca, 'ColorOrder', [0 0 0; 1 0 0; 0 0 1; 0 1 0], 'NextPlot', 'replacechildren');
    plot(f, M, 'Linew', 2)
    h = legend('total (sine+noise)', 'incoherent', 'coherent', 'vol.cond.');
    set(h, 'Interpreter', 'Latex');
    ylim([1e-9 5e-6])
end


% Set Axes
set(gca, 'Ysca', 'log')
xlim([1 60])
xlabel('f [Hz]', 'Interpreter', 'Latex')
ylabel('PSD [V$^2$/Hz]', 'Interpreter', 'Latex')
tit_str = ['$\omega_0=' num2str(w0) ', n_\sigma=' num2str(nsig) '$'];
title(tit_str, 'Interpreter', 'Latex');
% legend([h1.patch h2.patch h3.patch h4.patch h5], ...
%     'Sine+Noise', 'Inc', 'Coh', 'VC', 'Sine')

% Save and close
print(fig1, figname, '-depsc')
close(fig1)