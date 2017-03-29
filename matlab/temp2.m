f=logspace(-1,2,30);
dt=1/2500;
for i=1:length(ndat); 
    for j=1:3;
        ni=T{i}{j}>=t(i,1) & T{i}{j}<=t(i,2);
        [LFPf(:,j,i), WLFP(:,:,j,i), coi, PLFP(:,j,i)]=procdata(data{i}{j}(ni));
    end
    for j=1:3;
        ni=T{i}{j+3}>=t(i,1) & T{i}{j+3}<=t(i,2);
        [EMGf(:,j,i), WEMG(:,:,j,i), coi, PEMG(:,j,i)]=procdata(data{i}{j+3}(ni));
    end
    for j=1:2;
        nia=T{i}{j}>=t(i,1) & T{i}{j}<=t(i,2);
        nib=T{i}{j+1}>=t(i,1) & T{i}{j+1}<=t(i,2);
        [dLFPf(:,j,i), dWLFP(:,:,j,i), coi, dPLFP(:,j,i)]=...
            procdata(data{i}{j}(nia)-data{i}{j+1}(nib));
    end
    [CLE(:,:,:,:,i), WLE(:,:,:,:,i)] = wave_coh(WLFP(:,:,:,i), WEMG(:,:,:,i),f,dt);
    [dCLE(:,:,:,:,i), dWLE(:,:,:,:,i)] = wave_coh(dWLFP(:,:,:,i),...
        WEMG(:,:,:,i),f,dt);
end


%% Plot resulting PSD
for i=1:3;
    subplot(3,1,i)
    loglog(f,squeeze(PLFP(:,i,:)))
%     semilogx(f,PLFP(:,i,2)./PLFP(:,i,1),f,PLFP(:,i,3)./PLFP(:,i,1),...
%         f,PLFP(:,i,4)./PLFP(:,i,1),f,PLFP(:,i,5)./PLFP(:,i,1));
    xlabel('f [Hz]')
    ylabel('PSD V^2/Hz]')
%     ylabel('P_i/P_1')
end

% figure
% for i=1:2;
%     subplot(2,1,i)
%     semilogx(f,dPLFP(:,i,2)./dPLFP(:,i,1),f,dPLFP(:,i,3)./dPLFP(:,i,1),...
%         f,dPLFP(:,i,4)./dPLFP(:,i,1),f,dPLFP(:,i,5)./dPLFP(:,i,1)); 
%     xlabel('f [Hz]')
%     ylabel('dP_i/dP_1')
% end
% 
figure
for i=1:3;
    subplot(3,1,i)
    loglog(f,squeeze(PEMG(:,i,:)))
%     semilogx(f,PEMG(:,i,2)./PEMG(:,i,1),f,PEMG(:,i,3)./PEMG(:,i,1),...
%         f,PEMG(:,i,4)./PEMG(:,i,1),f,PEMG(:,i,5)./PEMG(:,i,1)); 
    xlabel('f [Hz]')
    ylabel('PSD V^2/Hz]')
%     ylabel('P_i/P_1')
end

% %% 1D Coherency spectra
% for s=1:5;
%     for i=1:3;
%         for j=1:3;
%             C1d(:,i,j)=nanmean(coi2nan(f,squeeze(CLE(:,:,i,j,s)),coi),2);
%         end
%     end
% end

% %% Coherency-Phase Scalogram
% s=1; % "site" number
% k=1;
% for i=1:3;
%     for j=1:3; 
%         subplot(3,3,k)
%         plot_coh_phase(ti,f,squeeze(CLE(:,:,i,j,s)),WLE(:,:,i,j,s),'ahead',2);
%         k=k+1;
%     end
% end

clear i  j