%% Plot the outcome of several trials trying to replicate PSD of synthetic time series

%% PCC
fname = 'PSD_synth_PCC_ns6_w12_nl3_v2';
load(fname)
tit_str = ['PCC, $\omega_0=' num2str(w0) ', n_\sigma=' num2str(nsig) '$'];

%% IC
% fname = 'PSD_synth_IC_p05_ns6_w12_nl3';
% load(fname)
% tit_str = ['IC, $p=0.05, \omega_0=' num2str(w0) ', n_\sigma=' num2str(nsig) '$'];

%% PLI
% fname = 'PSD_synth_PLI_p05_ns6_w12_nl3';
% load(fname)
% tit_str = ['PLI, $p=0.05, \omega_0=' num2str(w0) ', n_\sigma=' num2str(nsig) '$'];

%%


fig1 = figure('Papersize', [16 10], 'PaperPosition', [0.75 0.5 14.5 9], ...
        'PaperPositionmode', 'manual', 'Visible', 'on'); 

    
%% Parameters
ds    = 1;
ds2   = 150; %downsampling for plot
dt    = 1/2500;
nt    = 20/dt;
t     = (0:nt-1)*dt;
t     = t(1:ds:end);
limy  = [1e-8 1e-6];



%% PCC
M = [mean(Ptot,2) mean(Pinc,2) mean(Pcoh,2) mean(Pvc,2)+1e-19]';
lineprops.col={'k'; 'r'; 'b'; 'g'};
E = [ min([mean(Ptot,2)'-1e-20; std(Ptot,0,2)']); ...
    min([mean(Pinc,2)'-1e-20; std(Pinc,0,2)']);...
    min([mean(Pcoh,2)'-1e-20; std(Pcoh,0,2)']);...
    min([mean(Pvc, 2)'-1e-20; std(Pvc, 0,2)']) ]; %% CHECKEN!

%% IC or PLI
% M = [mean(Ptot,2) mean(Pinc,2) mean(Pcoh,2)]';
% lineprops.col={'k'; 'r'; 'b'};
% E = [ std(Ptot,0,2) std(Pinc,0,2) std(Pcoh,0,2) ]'; %% CHECKEN!
% E = min(cat(3, E, M-(1e-10)),[],3);

%%

mseb(f, M, E, lineprops);
xlim([1 60])
hold all
xa    = [zeros(nt/5,1); ...
         sin(2*pi*alpha*(1:3*nt/5)'*dt); zeros(nt/5,1)];
xb1   = [sin(2*pi*beta*(1+(1:nt/2)'/nt*0.25).*(1:nt/2)'*dt); ...
         sin(2*pi*beta*(1+(nt/2+1:nt)'/nt*0.25).*(nt/2+1:nt)'*dt-dp)*0];
xb2   = [sin(2*pi*beta*(1+(1:nt/2)'/nt*0.25).*(1:nt/2)'*dt)*0; ...
         sin(2*pi*beta*(1+(nt/2+1:nt)'/nt*0.25).*(nt/2+1:nt)'*dt-dp)];
xg    = [sin(2*pi*gamma*t(1:nt/2)'+dp); ...
         sin(2*pi*gamma*t(nt/2+1:end)')*0];
xa    = xa/A;
xb1   = xb1/A;
xb2   = xb2/A;
xg    = xg/A;

clear Psig
[~,~,~,Psig(:,1)] = procdata(xa , 'freq', f, 'w0', w0, 'filt', [], 'dt', dt);
[~,~,~,Psig(:,2)] = procdata(xb1, 'freq', f, 'w0', w0, 'filt', [], 'dt', dt);
[~,~,~,Psig(:,3)] = procdata(xb2, 'freq', f, 'w0', w0, 'filt', [], 'dt', dt);
[~,~,~,Psig(:,4)] = procdata(xg , 'freq', f, 'w0', w0, 'filt', [], 'dt', dt);
loglog(f,Psig,'k--');

%% PCC
h = legend('total (sine+noise)', 'incoherent', 'coherent', ...
    'vol.cond.', 'sine w/o noise');

%% IC or PLI
% h = legend('total (sine+noise)', 'incoherent', 'coherent', ...
%     'sine w/o noise');

%%
set(h, 'Interpreter', 'Latex');
ylim([1e-8 1e-6])

% Set Axes
set(gca, 'Ysca', 'log')
xlim([1 70])
xlabel('f [Hz]', 'Interpreter', 'Latex')
ylabel('PSD [V$^2$/Hz]', 'Interpreter', 'Latex')
title(tit_str, 'Interpreter', 'Latex');

% Save and close
print(fig1, fname, '-depsc')
% close(fig1)
set(fig1, 'Visi', 'On')