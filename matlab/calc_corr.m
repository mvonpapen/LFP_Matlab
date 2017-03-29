function [C_off, C_on] = calc_corr(task)

%% Load LFP data
load(['data_' task '_v3']);



%% Parameters
Nf      = 60;
f       = logspace(0,3,Nf);
w0      = 12;
ds      = 10;
Nex     = length(LFP_OFF);



%% Preallocate
C_off = NaN( Nf, Nf, 6, Nex );
C_on  = NaN( Nf, Nf, 6, Nex );



%% Calculate correlation matrices
for i=1:Nex
        
    fprintf('Processing data set %u/%u.\n', i, Nex);

    % Determine bipolar combinations
    combo = nchoosek(1:Nch(i),2);
    Nbc(i)= size(combo,1);
    BiPolCh{i} = Channel{i}(combo);

    
    % OFF
    x = LFP_OFF{i}(:,combo(:,1))-LFP_OFF{i}(:,combo(:,2));
    a = cell(Nbc(i),1);
    for j=1:Nbc(i);
        a{j} = [ArtLFP_OFF{i}{combo(j,1)}; ArtLFP_OFF{i}{combo(j,2)}];
    end
    [~, W, coi] = procdata(x, 'freq', f, 'w0', w0, 'filter', [], 'art', a);
    for n=1:Nbc(i)
        P = log10(2/2500*abs(W(:,1:ds:end,n)).^2);
%         P = coi2nan(f,P,coi(1:ds:end,n));
        C_off(:,:,n,i) = corr(P');
    end
    
    % ON
    x = LFP_ON{i}(:,combo(:,1))-LFP_ON{i}(:,combo(:,2));
    a = cell(Nbc(i),1);
    for j=1:Nbc(i);
        a{j} = [ArtLFP_ON{i}{combo(j,1)}; ArtLFP_ON{i}{combo(j,2)}];
    end
    [~, W, coi] = procdata(x, 'freq', f, 'w0', w0, 'filter', [], 'art', a);
    for n=1:Nbc(i)
        P = log10(2/2500*abs(W(:,1:10:end,n)).^2);
%         P = coi2nan(f,P,coi(1:ds:end,n));
        C_on(:,:,n,i) = corr(P');
    end
    
    clear x a W coi P
    
end