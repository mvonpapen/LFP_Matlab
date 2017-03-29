%% Structure function test
clear all
load '/afs/geo.uni-koeln.de/usr/neuro/DATA/DATA_last.mat'
[t, dat, ndat, rigor, CHlab]=load_dat_from_DATA(DATA, 'value1', 'Ruhe',...
    'activity', 2);

ns=length(ndat);
dt=1/2500;
fit=[2 500]; %in Hz
nEMG=3;
for i=1:ns
    nLFP(i)=length(CHlab{i})-nEMG;
end
M=1:7; %max order
tau=unique(round(logspace(log10(1/min(fit)/dt),log10(1/max(fit)/dt),20)));

%% Begin and end of experiments
i0   = [1  6 11 17 25 31 35 41 50];
iend = [5 10 16 24 30 34 40 49 56];

%% Calculate structure functions
j=1;
for n=1:ns;
%     if nLFP(n)==1
%         continue
%     end
%     %% Check time vector
    for ch=1:nLFP(n)
%         ti(ch,:)=minmax(t{n}{ch});
%     end
%     ti=[max(ti(:,1)) min(ti(:,2))];
%     for ch=1:nLFP(n)-1;
%         ni1 =  t{n}{ch}>=ti(1) & t{n}{ch}<=ti(2);
%         ni2 =  t{n}{ch+1}>=ti(1) & t{n}{ch+1}<=ti(2);
%         x=dat{n}{ch}(ni1)-dat{n}{ch+1}(ni2);
        x=dat{n}{ch};
        for m=M;
            i=1;
            for l=tau;
                SL(i,m,j,ch)=nanmean(abs(x(l+1:end)-x(1:end-l)).^m);
                i=i+1;
            end
        end
    end
    j=j+1;
%     for ch=1:nEMG;
%         x=dat{n}{nLFP(n)+ch};
%         for m=M;
%             i=1;
%             for l=tau;
%                 SE(i,m,n,ch)=mean(abs(x(l+1:end)-x(1:end-l)).^m);
%                 i=i+1;
%             end
%         end
%     end
end

%% Calculate Slope of structure functions
for n=1:ns;
    for ch=1:nLFP(n)
        for m=M;
            xiL(m,:,n,ch)=fitkappa(1./(tau*dt),SL(:,m,n,ch),fit);
        end
    end
%     for ch=1:nEMG;
%         for m=M;
%             xiE(m,:,n,ch)=fitkappa(1./(tau*dt),SE(:,m,n,ch),fit);
%         end
%     end
end

%% Plot results: xi(m) vs. m
for i=1:length(i0);
    for j=1:nLFP(i0(i))
        plot(M,squeeze(xiL(:,1,i0(i):iend(i),j)),'-')
        hold all
    end
end
% for i=1:nEMG;
%     subplot(2,max([max(nLFP) nEMG]),i+max([max(nLFP) nEMG]))
%     if i<=2
%         plot(1:M,squeeze(xiE(:,1,1:5,i)),'-')
%         hold all
%     end
%     plot(1:M,squeeze(xiE(:,1,6:10,i)),'--',...
%         1:M,squeeze(xiE(:,1,11:15,i)),'.-')
% end

% figure
% plot(rigor(1:5), squeeze(xiL(M,1,1:5,1)),'-+', ...
%     rigor(6:10), squeeze(xiL(M,1,6:10,1)),'-+', ...
%     rigor(11:15), squeeze(xiL(M,1,11:15,1)),'-+', ...
%     rigor(16:23), squeeze(xiL(M,1,16:23,1)),'-+', ...
%     rigor(24:29), squeeze(xiL(M,1,24:29,2)),'-+',...
%     rigor(30:38), squeeze(xiL(M,1,30:38,1)),'-+')
% plot(rigor(1:5), squeeze(xiL(M,1,1:5,1)),'+', ...
%     rigor(1:5), squeeze(xiL(M,1,1:5,2)),'+', ...
%     rigor(6:10), squeeze(xiL(M,1,6:10,1)),'+', ...
%     rigor(6:10), squeeze(xiL(M,1,6:10,2)),'+', ...
%     rigor(6:10), squeeze(xiL(M,1,6:10,3)),'+', ...
%     rigor(11:15), squeeze(xiL(M,1,11:15,1)),'+', ...
%     rigor(11:15), squeeze(xiL(M,1,11:15,2)),'+', ...
%     rigor(11:15), squeeze(xiL(M,1,11:15,3)),'+', ...
%     rigor(16:23), squeeze(xiL(M,1,16:23,1)),'+', ...
%     rigor(24:29), squeeze(xiL(M,1,24:29,2)),'+', ...
%     rigor(24:29), squeeze(xiL(M,1,24:29,3)),'+')

%% Boxplot xiL(m) against time after apo-injection

% %for dLFP
% i0=[1 6 11 16];
% iend=[5 10 15 21];

for m=4
    HH=NaN(1,max(iend-i0)+1);
    n=1;
    for i=1:length(i0);
        for j=1:max(nLFP);
            HH(n,1:(iend(i)-i0(i)+1))=...
                squeeze(xiL(m,1,i0(i):iend(i),j));
            n=n+1;
        end
    end
    HH(HH==0)=NaN;
    figure
    boxplot(HH)
    xlim([0.5 9.5])
end

%% Boxplot with grouped variables
for m=2;
    xi=squeeze(xiL(m,1,:,:));
    xi(xi==0)=NaN;
    G=NaN(size(repmat(rigor',1,3)));
    G(rigor>=60,:)=2;
    G(rigor>=40 & rigor<60,:)=3;
    G(rigor>=20 & rigor<40,:)=2;
    G(rigor<30,:)=1;
    G=G(:); xi=xi(:);
    figure
    boxplot(xi,G)
    clear xi x ch n m
end