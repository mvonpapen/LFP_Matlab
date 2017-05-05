%% Plots power composition

task = 'rest';

load(['LFP_pow_v3_' task])
tit = ['Percentage of total signal during ' task];
col = 1;
STNact = 1;


% % Specify patients
% for i=1:length(Patient)
%     tmp(i)=any(strcmp(Patient{i}, {'AUIN_L', 'BEMI_R', 'GRFR_L'}));
% end
% P_tot_off(:,:,~tmp) = NaN;
% P_tot_on(:,:,~tmp) = NaN;
% clear tmp

% P_tot_off = P_coh_off + P_inc_off + P_vc_off;
% P_tot_on  = P_coh_on + P_inc_on + P_vc_on;
    
% Only use channels with STN activity
for i=1:length(LFPact)
    P_tot_off(:,LFPact{i}<STNact,i) = NaN;
    P_inc_off(:,LFPact{i}<STNact,i) = NaN;
    P_coh_off(:,LFPact{i}<STNact,i) = NaN;
    P_vc_off(:,LFPact{i}<STNact,i)  = NaN;
    P_bip_off(:,LFPact{i}<STNact,i) = NaN;
end

Rcoh_off = P_coh_off(:,:)./P_tot_off(:,:)*100;
Rinc_off = P_inc_off(:,:)./P_tot_off(:,:)*100;
Rvc_off  = P_vc_off(:,:) ./P_tot_off(:,:)*100;
    
Rcoh_on = P_coh_on(:,:)./P_tot_on(:,:)*100;
Rinc_on = P_inc_on(:,:)./P_tot_on(:,:)*100;
Rvc_on  = P_vc_on(:,:) ./P_tot_on(:,:)*100;
  
%% MSEB Plot
subplot(1,2,1)
N = sqrt(sum(~isnan(Rvc_off),2)'); %ones(1,length(f));
M = [nanmean(Rinc_off,2) nanmean(Rcoh_off,2) nanmean(Rvc_off,2)]';
% M = [nanmean(Rinc_off,2) nanmean(Rcoh_off+Rvc_off,2)]';
lineprops.col={'r'; 'b'; 'g'};
E = [ nanstd(Rinc_off,0,2)'; nanstd(Rcoh_off,0,2)'; nanstd(Rvc_off, 0,2)']./N([1 1 1],:);
% E = [ nanstd(Rinc_off,0,2)'; nanstd(Rcoh_off+Rvc_off, 0,2)']./N([1 1],:);
mseb(f, M, E, lineprops);
ylim([0 100]), xlim([1 85])
ylabel('P_i/P_{tot} [%]')
xlabel('f [Hz]')
title([tit ' OFF'])

subplot(1,2,2)
M = [nanmean(Rinc_on,2) nanmean(Rcoh_on,2) nanmean(Rvc_on,2)]';
% M = [nanmean(Rinc_on,2) nanmean(Rcoh_on+Rvc_on,2)]';
lineprops.col={'r'; 'b'; 'g'};
E = [ nanstd(Rinc_on,0,2)'; nanstd(Rcoh_on,0,2)'; nanstd(Rvc_on, 0,2)']./N([1 1 1],:);
% E = [ nanstd(Rinc_on,0,2)'; nanstd(Rcoh_on+Rvc_on, 0,2)']./N([1 1],:);
mseb(f, M, E, lineprops);
ylim([0 100]), xlim([1 85])
ylabel('P_i/P_{tot} [%]')
xlabel('f [Hz]')
title([tit ' ON'])
legend('incoherent', 'coherent','volume-conducted')
 
%% ApoPaper
% subplot(1,2,1)
% N = sqrt(sum(~isnan(Rinc_off),2));
% h1=shadedErrorBar(f,nannanmean(Rinc_off,2), nannanstd(Rinc_off,0,2)./N,'-r',1);
% hold all
% h2=shadedErrorBar(f,nannanmean(Rcoh_off,2), nannanstd(Rcoh_off,0,2)./N,'-b',1);
% h3=shadedErrorBar(f,nannanmean(Rvc_off,2), nannanstd(Rvc_off,0,2)./N,'-g',1);
% ylim([0 100]), xlim([1 85])
% ylabel('P_i/P_{tot} [%]')
% xlabel('f [Hz]')
% title([tit ' OFF'])
% 
% subplot(1,2,2)
% N = sqrt(sum(~isnan(Rinc_on),2));
% shadedErrorBar(f,nannanmean(Rinc_on,2), nannanstd(Rinc_on,0,2)./N,'-r',1)
% hold all
% shadedErrorBar(f,nannanmean(Rcoh_on,2), nannanstd(Rcoh_on,0,2)./N,'-b',1)
% shadedErrorBar(f,nannanmean(Rvc_on,2), nannanstd(Rvc_on,0,2)./N,'-g',1)
% ylim([0 100]), xlim([1 85])
% xlabel('f [Hz]')
% ylabel('P_i/P_{tot} [%]')
% title([tit ' ON'])
% legend([h1.patch h2.patch h3.patch], ...
%     'Incoherent', 'Coherent','Volume-conducted')