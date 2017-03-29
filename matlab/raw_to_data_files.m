%% Test mono- and bipolar PSD changes and coherence from OFF to ON


%% Load LFP data

for TYPE = {'Ruhe'}; %{'Faust', 'Halte', 'Ruhe'};
    type = TYPE{1};
    
    [t, dat, CH, art, rigor, ndat, ~, LFPactivity] = ...
        load_dat_from_DATA ( 'experiment', type, 'activity', 0, 'commontime', 1, 'noEMG', 0);
    
    load DATA_last.mat
    tag     = '2nd';



    %% Determine sites and analysis
    switch type
        case 'Ruhe'
            % Sites with rest ON OFF, Nelec>1 and rigor>=30
            i0   = [1  6 11 26 32 41 56 65]+1; % last entry: no EMG channels
            iend = [5 10 17 31 68 46 62 67];
        case 'Halte'
            % Sites with hold ON OFF and nLFP>1
            i0   = [1 3 5  9 11 17 21];
            iend = [2 4 6 10 23 18 22];
        case 'Faust'
            % Sites with fist ON OFF and rigor>=30
            i0   = [1 3 8 10 16 18];
            iend = [2 4 9 20 17 19];
    end
    
    Nex = length(i0);
    clear Patient
    LFPch = {'C', 'L', 'A', 'M', 'P'};
    EMGch = {'EDCre', 'FDLre', 'FDIre', 'EDCli', 'FDLli', 'FDIli'};
    for i=1:Nex
        Channel{i}=CH{iend(i)};
        iLFP = find(ismember(Channel{i}, LFPch));
        iEMG = find(ismember(Channel{i}, EMGch));
        
        LFP_OFF{i}=dat{i0(i)}(:,iLFP);
        EMG_OFF{i}=dat{i0(i)}(:,iEMG);
        LFP_ON{i}=dat{iend(i)}(:,iLFP);
        EMG_ON{i}=dat{iend(i)}(:,iEMG);
        
        for j=1:length(iLFP)
            ArtLFP_OFF{i}{j}=art{i0(i)}{iLFP(j)};
            ArtLFP_ON{i}{j}=art{iend(i)}{iLFP(j)};
        end
        for j=1:length(iEMG)
            ArtEMG_OFF{i}{j}=art{i0(i)}{iEMG(j)};
            ArtEMG_ON{i}{j}=art{iend(i)}{iEMG(j)};
        end
            
        Patient{i}=DATA(ndat(i0(i))).patient;
        Patient2{i}=DATA(ndat(iend(i))).patient;
        Nch(i) = length(iLFP);
        LFPact{i}=LFPactivity{iend(i)};
        T_off(i) = max(t{i0(i)}) - min(t{i0(i)});
        T_on(i)  = max(t{iend(i)}) - min(t{iend(i)});
        
        
        if ~strcmp(Patient{i},Patient2{i})
            error('ERROR: i0 and iend refer to different patients for i=%i!', i)
        end
    end
    rigor=rigor(iend);
    
    clear i i0 iend DATA dat art CH LFPactivity t ndat pat Patient2 x Nex j LFPch EMGch iEMG iLFP

    save(['data_' type '.mat'])
    
end %type for loop end
