%%  Plot relative and absolute (percentual) ON-OFF difference of PSDs

type = 'hold';
load(['pow_coh_ds10_pc15_w6n3_w0_06_' type '.mat'])
fnam = ['psd_reldiff_peak_w6n3_' type '.eps'];
label = {'\theta', '\alpha', '\beta', '\beta1', '\beta2', '\gamma'};


nch = 4;
cutoff = 1e-12;
P_inc_off(P_inc_off<=cutoff) = 0;
P_inc_on(P_inc_on<=cutoff)   = 0;
P_coh_off(P_coh_off<=cutoff) = 0;
P_coh_on(P_coh_on<=cutoff)   = 0;
P_vc_off(P_vc_off<=cutoff)   = 0;
P_vc_on(P_vc_on<=cutoff)     = 0;


fband{1} = find(f>=4  & f<=7 );
fband{2} = find(f>=7  & f<=13);
fband{3} = find(f>=13 & f<=30);
fband{4} = find(f>=13 & f<=20);
fband{5} = find(f>=20 & f<=30);
fband{6} = find(f>=60 & f<=85);
Nb = length(fband);

fig1 = figure('PaperSize', [20 14], ...
        'PaperPositionmode', 'manual', 'PaperPosition', [1 1 18 12], ...
        'Visible', 'off');



%%%%%%%% vvv RELATIVE DIFFERENCES vvv %%%%%%%%%%%%%%%
% Specify patients
clear tmp
for i=1:length(Patient)
    tmp(i)=any(strcmp(Patient{i}, {'P01_L', 'P02_R', 'P04_L'}));
end
P_tot_off(:,:,tmp) = NaN;
P_tot_on(:,:,tmp) = NaN;


% Total monopolar
totx=P_tot_off(:,:);
toty=P_tot_on(:,:);


% Coherent
x=P_coh_off(:,:);
y=P_coh_on(:,:);
for i=1:Nb
    dPrel_coh(:,i) = squeeze(nanmean(y(fband{i},:)./toty(fband{i},:) ...
                                   -x(fband{i},:)./totx(fband{i},:),1))*100;
end

% Volume Conduction
x=P_vc_off(:,:);
y=P_vc_on(:,:);
for i=1:Nb
    dPrel_vc(:,i) = squeeze(nanmean(y(fband{i},:)./toty(fband{i},:) ...
                                   -x(fband{i},:)./totx(fband{i},:),1))*100;
end

% Incoherent
x=P_inc_off(:,:);
y=P_inc_on(:,:);
for i=1:Nb
    dPrel_inc(:,i) = squeeze(nanmean(y(fband{i},:)./toty(fband{i},:) ...
                                    -x(fband{i},:)./totx(fband{i},:),1))*100;
end


% Wilcoxon sign-rank test
for i=1:Nb
    sig_coh(i) = signrank ( dPrel_coh(:,i) );
    sig_inc(i) = signrank ( dPrel_inc(:,i) );
    sig_vc(i)  = signrank ( dPrel_vc(:,i)  );
end


% Monopolar inc
subplot(3,2,1)
boxplot( dPrel_inc, ...
    'labels', label, 'whisker', 50)
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex');
set(gca,'TickLabelInterpreter','tex');
ylim([-1 1]*50), grid on
ylabel('\DeltaP_{rel} [%]')
title('Monopolar incoherent')
hold all
for i=1:Nb
    if sig_inc(i)<0.01;        
        plot([i-0.05 i+0.05],[0.9*50 0.9*50],'*k');
    elseif sig_inc(i)<0.05;
        plot(i,0.9*50,'*k');
    end
end
% Monopolar coherent
subplot(3,2,3)
boxplot( dPrel_coh, ...
    'labels', label, 'whisker', 50)
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex');
set(gca,'TickLabelInterpreter','tex');
ylim([-1 1]*50), grid on
ylabel('\DeltaP_{rel} [%]')
title('Monopolar coherent')
hold all
for i=1:Nb
    if sig_coh(i)<0.01;        
        plot([i-0.05 i+0.05],[0.9*50 0.9*50],'*k');
    elseif sig_coh(i)<0.05;
        plot(i,0.9*50,'*k');
    end
end
% Monopolar volume conduction
subplot(3,2,5)
boxplot( dPrel_vc, ...
    'labels', label, 'whisker', 50)
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex');
set(gca,'TickLabelInterpreter','tex');
ylim([-1 1]*50), grid on
ylabel('\DeltaP_{rel} [%]')
title('Monopolar, volume cond.')
hold all
for i=1:Nb
    if sig_vc(i)<0.01;        
        plot([i-0.05 i+0.05],[0.9*50 0.9*50],'*k');
    elseif sig_vc(i)<0.05;
        plot(i,0.9*50,'*k');
    end
end



%%%%%%%% vvv ABSOLUTE DIFFERENCES vvv %%%%%%%%%%%%%%%


% Incoherent
x=P_inc_off(:,:);
y=P_inc_on(:,:);
for i=1:Nb
    dPabs_inc(:,i) = squeeze( nanmean( y(fband{i},:)-x(fband{i},:), 1 ) ...
                            ./ nanmean(totx(fband{i},:), 1 ) )*100;
end

% Coherent
x=P_coh_off(:,:);
y=P_coh_on(:,:);
for i=1:Nb
    dPabs_coh(:,i) = squeeze( nanmean( y(fband{i},:)-x(fband{i},:), 1 ) ...
                            ./ nanmean( totx(fband{i},:), 1 ) );
end

% Volume Conduction
x=P_vc_off(:,:);
y=P_vc_on(:,:);
for i=1:Nb
    dPabs_vc(:,i)  = squeeze( nanmean( y(fband{i},:)-x(fband{i},:), 1 ) ...
                            ./ nanmean( totx(fband{i},:), 1 ) );
end


% Wilcoxon sign-rank test
for i=1:Nb
    sig_inc(i) = signrank ( dPabs_inc(:,i) );
    sig_coh(i) = signrank ( dPabs_coh(:,i) );
    sig_vc(i)  = signrank ( dPabs_vc(:,i)  );
end


% Monopolar inc
subplot(3,2,2)
boxplot( dPabs_inc, ...
    'labels', label, 'whisker', 50)
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex');
set(gca,'TickLabelInterpreter','tex');
ylim([-1 1]*50), grid on
ylabel('\DeltaP_{abs} [%]')
title('Monopolar incoherent')
hold all
for i=1:Nb
    if sig_inc(i)<0.01;        
        plot([i-0.05 i+0.05],[0.9*50 0.9*50],'*k');
    elseif sig_inc(i)<0.05;
        plot(i,0.9*50,'*k');
    end
end
% Monopolar coherent
subplot(3,2,4)
boxplot( dPabs_coh*100, ...
    'labels', label, 'whisker', 50)
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex');
set(gca,'TickLabelInterpreter','tex');
ylim([-1 1]*50), grid on
ylabel('\DeltaP_{abs} [%]')
title('Monopolar coherent')
hold all
for i=1:Nb
    if sig_coh(i)<0.01;        
        plot([i-0.05 i+0.05],[0.9*50 0.9*50],'*k');
    elseif sig_coh(i)<0.05;
        plot(i,0.9*50,'*k');
    end
end
% Monopolar volume conduction
subplot(3,2,6)
boxplot( dPabs_vc*100, ...
    'labels', label, 'whisker', 50)
h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex');
set(gca,'TickLabelInterpreter','tex');
ylim([-1 1]*50), grid on
ylabel('\DeltaP_{abs} [%]')
title('Monopolar, volume cond.')
hold all
for i=1:Nb
    if sig_vc(i)<0.01;        
        plot([i-0.05 i+0.05],[0.9*50 0.9*50],'*k');
    elseif sig_vc(i)<0.05;
        plot(i,0.9*50,'*k');
    end
end

% annotation(fig1,'textbox', [0.32 0.96 0.4 0.03], 'String', type, ...
%         'Fontsize', 16, 'HorizontalAlign', 'center', 'Fontw', 'bold', ...
%         'EdgeColor', 'none');
% annotation(fig1,'textbox', [0.05 0.94 0.5 0.03], 'String', 'Relative', ...
%     'Fontsize', 12, 'HorizontalAlign', 'center', 'Fontw', 'bold', ...
%     'EdgeColor', 'none');
% annotation(fig1,'textbox', [0.5 0.94 0.5 0.03], 'String', 'Absolute', ...
%     'Fontsize', 12, 'HorizontalAlign', 'center', 'Fontw', 'bold', ...
%     'EdgeColor', 'none');
print(fig1, fnam, '-depsc')
close(fig1)