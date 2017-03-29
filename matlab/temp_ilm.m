%% Incoherent locality measure

clear all
load '/afs/geo.uni-koeln.de/usr/neuro/DATA/DATA_last.mat'
[t, dat, ndat, rigor, CHlab, art]=load_dat_from_DATA ( DATA, 'value1', 'Ruhe', ...
    'activity', 1, 'value2', 'PECH_L');


n=5;
nt=length(t{n}{1});
nf=50;
f =logspace(0,2.7,nf);
[dj,s0,j1] = scale4wavelet(f);
dt=1/2500;
sig=0.5;
nch=length(dat{n})-3;

clear LFP W P0 coi p s
for i=1:nch;
    LFP(:,i)=dat{n}{i};
    [W(:,:,i),p,s,coi(:,i)]=wavelet(LFP(:,i),dt,0,dj,s0,j1);
%     P0(:,i)=mygws(abs(W(:,:,i)).^2,dt,coi(:,i),f);
end
[C,Wxy]=wave_coh(W,W,f,dt,3,0);
Ph=angle(Wxy);

clear ts P
iLC=indvc(C,Ph);
for i=1:nch
    [ts(:,i),P(:,i)]=wave_recstr (W(:,:,i), f, iLC);
end

% F=repmat(f(:),[1 nt]);
% iall=find(C(:,:,1,2)>=sig & C(:,:,1,3)>=sig & C(:,:,2,3)>=sig ...
%     & ~(F>=35 & F<=70) & ~(F>=71 & F<=140));
% for i=1:3;
%     [tsall(:,i),Pall(:,i)]=wave_recstr (W(:,:,i), f, iall);
% end

dLFP=[ts(:,2)-ts(:,1) ts(:,3)-ts(:,1)];
eigvec=princomp(dLFP)

% 
% 
% figure
% semilogx(f,1-P(:,[1 3])./P0(:,[1 1]),f,1-P(:,[2 5])./P0(:,[2 2]), ...
%     f,1-P(:,[4 6])./P0(:,[3 3]))
% legend('1-(C_{CA}/C)','1-(A_{CA}/A)','1-(C_{CL}/C)','1-(L_{CL}/L)', ...
%     '1-(A_{AL}/A)','1-(L_{AL}/L)')
% xlabel('f [Hz]')
% ylabel('ILM: 1-(P_{coh}/P_{total})')
% xlim([0.9 510])
% ylim([0 1])