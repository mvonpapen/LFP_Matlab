%% Plots differences of relative power from OFF to ON as mseb-plot
%% use pow_coh_Pw0_**_type.mat as input
task = 'rest';
nB = 1;


load(['pow_coh_ds10_pc15_w0_12_' task])
tit = {'Incoherent', 'Coherent', 'Vol.cond.'};

coh_off = P_coh_off(:,:)./P_tot_off(:,:)*100;
inc_off = P_inc_off(:,:)./P_tot_off(:,:)*100;
vc_off  = P_vc_off(:,:) ./P_tot_off(:,:)*100;
    
coh_on = P_coh_on(:,:)./P_tot_on(:,:)*100;
inc_on = P_inc_on(:,:)./P_tot_on(:,:)*100;
vc_on  = P_vc_on(:,:) ./P_tot_on(:,:)*100;

clear X
X(:,:,3) =  vc_on'- vc_off';
X(:,:,2) = coh_on'-coh_off';
X(:,:,1) = inc_on'-inc_off';
X(abs(X)<1) = 0;
col={{'r'}, {'b'}, {'g'}};

%% MSEB Plot
figure
for i=1:3
    subplot(1,3,i)
%     x = smoothts(X(:,:,i),'b',nB);
    x = X(:,:,i);
    N = sqrt(sum(~isnan(x)));
    lineprops.col=col{i};
    mseb(f, nanmean(x), nanstd(x)./N, lineprops);
    ylim([-30 30]), xlim([1 85])
    ylabel('\Delta(P_i/P_{tot}) [%]')
    xlabel('f [Hz]')
    title([task ' - ' tit{i}])
    grid
    % significance
    for j=1:length(f)
        p(j) = signrank(x(:,j))*length(f)/nB;
    end
    if nB==1
        p = p/length(f);
    end
    hold all
    plot(f(p<0.05),nanmean(x(:,p<0.05)),'*k')
end