% Plots the average PSD of PCC in an mseb plot

function [h1, h2] = plot_avgPSD_PCC( task, fnam, comp, offon )

if nargin<4
    offon = 1:2;
end
if nargin<3
    comp = 1:3;
end
if nargin<2
    fnam = 'LFP_pow_PCC_v1_filt_';
end

colors1 = {'k'; 'r'; 'b'};
colors2 = {[0.5 0.5 0.5]; [1 0.5 0.5]; [0.5 0.5 1]};
legstr  = {'$P_\mathrm{tot}$(OFF)', '$P_\mathrm{inc}$(OFF)', '$P_\mathrm{coh}$(OFF)', ...
           '$P_\mathrm{tot}$(ON)',  '$P_\mathrm{inc}$(ON)',  '$P_\mathrm{coh}$(ON)'};
     
%% MSEB PLOT of PCC
load([fnam task '.mat'])
if any(ismember(offon,1))
    N = [sum(~isnan(P_tot_off(:,:)),2) sum(~isnan(P_inc_off(:,:)),2) ...
        sum(~isnan(P_coh_off(:,:)),2)];
    M = [nanmean(P_tot_off(:,:),2) nanmean(P_inc_off(:,:),2) nanmean(P_coh_off(:,:),2)]'+1e-10;
    E = [ nanstd(P_tot_off(:,:),0,2) nanstd(P_inc_off(:,:),0,2) ...
        nanstd(P_coh_off(:,:),0,2) ]'./sqrt(N'); %% CHECKEN!
    E = min(cat(3, E, M-(1e-10)),[],3);
    lineprops.style='-';
    lineprops.col={colors1{comp}};
    h1 = mseb(f, M(comp,:), E(comp,:), lineprops, 1);
    hold all
end

if any(ismember(offon,2))
    N = [sum(~isnan(P_tot_on(:,:)),2) sum(~isnan(P_inc_on(:,:)),2) ...
        sum(~isnan(P_coh_on(:,:)),2)];
    M = [nanmean(P_tot_on(:,:),2) nanmean(P_inc_on(:,:),2) nanmean(P_coh_on(:,:),2)]'+1e-10;
    E = [ nanstd(P_tot_on(:,:),0,2) nanstd(P_inc_on(:,:),0,2) ...
        nanstd(P_coh_on(:,:),0,2) ]'./sqrt(N'); %% CHECKEN!
    E = min(cat(3, E, M),[],3);
    lineprops.style='--';
    lineprops.col={colors2{comp}};
    h2 = mseb(f, M(comp,:), E(comp,:), lineprops, 1);
end

set(gca, 'YSca', 'log')
xlabel('f [Hz]', 'Interp', 'Latex')
ylabel('PSD [V$^2$/Hz]', 'Interp', 'Latex')

legind = comp' + (offon-1)*3;
h = legend({legstr{legind(:)}});
set(h, 'Interp', 'Latex')
xlim([1 45])
ylim([1e-9 4e-6])
title(['Avg. PSD for PCC during ' task])