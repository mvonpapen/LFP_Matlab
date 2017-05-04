%% Test explicitly how much energy goes into vc and coh for varying phase differences

f       = 15;  %sets coherent frequency of phase variation
nl      = 5;
dt      = 1/2500;
nt      = 20/dt;
t       = (1:nt)*dt;
x       = sin(2*pi*f*t') + nl * randn(nt,1);
ynoise  = nl*randn(nt,1);
sig     = 0.43;
nsig    = 6;
w0      = 12;
scale   = (w0+sqrt(2+w0^2))/4/pi ./ f;
phthres = 10;
dp      = [0:40];

j=1;
for i = dp
    y  = sin( 2*pi*f*t' + i/180*pi ) + ynoise;
    [~,W,coi,tmp] = procdata([x y],'freq',f,'filter',[], 'w0', w0);
    Ptot(j) = tmp(1);
    [C,Wxy]=wave_cohere(W,scale,nsig);
    [tmp1, ~, tmp3] = psd_acoh ( f, W, C, coi/nsig, sig, Wxy, 0, phthres );
    Pcoh(j) = tmp1(:,1,2);
    Pvc(j)  = tmp3(:,1,2);
    j=j+1;
end
clear tmp1 tmp3 W coi tmp Wxy


%% Plot results
plot(dp, Pcoh./Ptot, dp, Pvc./Ptot)
hold all, plot(phthres*[1 1], [0 1], '--k')
grid on
xlim([0 40]), ylim([0 1])
ylabel('Relative power (P_i/P_{tot})')
xlabel('\Delta\Phi between signals x and y')
title('Identified components for varying phase difference \Delta\Phi')
legend('coherent power (p<0.01)', 'volume conduction')
% Create textbox
annotation(gcf,'textbox',...
    [0.686575364667747 0.217261904761905 0.19997244732577 0.428571428571429],...
    'String',['$\Delta\Phi_c = 10^\circ$', sprintf('\n'), ...
              '$n_\sigma = ' num2str(nsig) '$', sprintf('\n'), ...
              '$\omega_0 = ' num2str(w0) '$', sprintf('\n'), ...
              'SNR$ = 1/50$',sprintf('\n'), ...
              '$f = ' num2str(f) '$ Hz', sprintf('\n'), ...
              '$T = ' num2str(nt*dt) '$ s',sprintf('\n'), ...
              '$\Delta t = 1/2500$ s'],...
    'FitBoxToText','on', 'Interpreter', 'Latex');