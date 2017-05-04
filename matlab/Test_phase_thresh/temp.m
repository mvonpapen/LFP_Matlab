% figure
% i=3; r=2;
% coh = repmat(Cs(5,:,i,r)>sig,50,1);
% 
% % coherence
% for j=1:3
%     tmp=CsLE(:,:,i,j);
%     tmp(~coh)=NaN;
%     c(:,1,j)=nanmean(tmp,2);
%     tmp=CsLE(:,:,i,j);
%     tmp(coh)=NaN; 
%     c(:,2,j)=nanmean(tmp,2);
% end
% subplot(1,2,1)
% % mseb(f, nanmean(c,3)', nanstd(c,[],3)'/sqrt(3));
% plot(f, c(:,:,2));
% xlim([1 40]), ylim([0 0.6])
% 
% % power
% tmp=WEs(:,:,2);
% tmp(~coh)=NaN;
% p(:,1)=nanmean(coi2nan(f,tmp,coiEs(:,2)),2);
% tmp=WEs(:,:,i);
% tmp(coh)=NaN; 
% p(:,2)=nanmean(coi2nan(f,tmp,coiEs(:,2)),2);
% subplot(1,2,2)
% semilogy(f,p,'+-'), xlim([1 40])

load data_rest_v3
f = [1:40];
%% Parameters
f       = logspace(0,1.85,100);
Nf      = length(f);
phthres = 15;
w0      = 12;
nsig    = 3;
sig     = sig_coh_thresh(w0, nsig);
usenan  = 0; % 1 = take mean only over times where signal is present
ds      = 5;
tag     = 'trem';
Nex     = length(LFP_OFF);
zero    = true; % define volume-conduction by 0°--phthres only (not 180°-phthres--180°)
scale   = (w0+sqrt(2+w0^2))/4/pi ./ f;
i=2;
[~,W1,coi1,P_tot_off(:,1:Nch(i),i)] = procdata(LFP_OFF{i}, 'freq', f, 'w0', w0, ...
'filter', [], 'art', ArtLFP_OFF{i}); %'filter', [44 56; 90 110],
[Cs, Wxys, W1s] = wave_cohere ( W1, scale, nsig, ds );
coi1s = coi1(1:ds:end,:)/nsig;
[Pcoh, Pinc, Pvc, Pnvc] = psd_acoh ( f, W1s, Cs, coi1s, sig, Wxys, usenan, phthres, zero );
Ph = angle(Wxys)/pi*180;

% figure
% for i=1:4; for j=1:4; [H(:,i,j),bin] = hist(Ph(5,:,i,j),[-180:5:180]); end, end
% for i=1:4; subplot(2,2,i), plot(bin,H(:,:,i)), end

i=2;j=3;
Pinc(Pinc<1e-10) = 1e-10;
Pcoh(Pcoh<1e-10) = 1e-10;
Pvc(Pvc<1e-10) = 1e-10;
t=[0:length(Wxys)-1]/2500*ds;

figure
X = LFP_OFF{2}(1:ds:end,:);
plot(t,[X(:,i)-0.02 X(:,j)+0.02])
figure
plot(t(500:1000),[X(500:1000,i)-0.02 X(500:1000,j)+0.02])

figure
semilogy(f,P_tot_off(:,i,2),'-k', f,Pinc(:,i,j),'-r', f,Pcoh(:,i,j),'b-', f,Pvc(:,i,j),'-g')
xlim([1 70])

figure
coh = Cs(:,:,i,j)>sig & abs(Ph(:,:,i,j))>phthres;
inc = Cs(:,:,i,j)<=sig;
vc = Cs(:,:,i,j)>sig & abs(Ph(:,:,i,j))<=phthres;
Z   = zeros(size(Wxys(:,:,i)));
Z(coh) = 1;
Z(vc ) = 2;
ds2=50;
pcolor(t(1:ds2:end),f,Z(:,1:ds2:end)), shading flat
hold all
plot(t(1:ds2:end),1./coi1s(1:ds2:end,i),'-k')
caxis([0 2])
colormap([1 0 0; 0 0 1; 0 1 0])
colorbar('YTickLabel',{'Inc'; 'Coh'; 'VC'}, 'YTick', [1/3 1 5/3])
title('PCC in time-frequency domain')