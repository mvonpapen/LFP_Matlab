% Sift and shift cluster analysis for PSD

function [ind, Pmean, minRMS] = ssca (PSD, n, range, logsca)

if nargin<4
    logsca=1;
end

[nf, nch] = size (PSD);
if logsca==1
    PSD = log10(PSD);
end
PSD(PSD==Inf) = NaN;
PSD(PSD==-Inf) = NaN;

if nargin<3
    range=1:nf;
end

R = rand (1, nch);
R(R==1) = R(R==1)-0.01;

Pmean = zeros (nf, n);
RMS   = zeros (n, nch);
ind   = zeros (1, nch);
indnew= ind;

%% First random step
for i=1:n
    indnew (R>=(i-1)/n & R<i/n) = i;
end

ite = 0;
while norm(ind-indnew)~=0
    
    ind = indnew;
    
    %% Calculate mean PSD
    for i=1:n
        member = PSD(:,ind==i);
        if numel(member) == nf
            Pmean(:,i) = member;
        else
            Pmean(:,i) = nanmean( member,2 );
        end
    end
    
    %% Calculate RMS for each PSD and each group
    for i=1:n
        for j=1:nch
            RMS(i,j) = sqrt( nanmean( (PSD(range,j)-Pmean(range,i)).^2 ) );
        end
    end

    %% Compare index ind with index of smallest RMS
    [minRMS, indnew] = min ( RMS );
    ite = ite+1;
    
end

if logsca==1
    Pmean = 10.^Pmean;
end

fprintf('After %u iterations, number of group members:', ite)
for i=1:n
    fprintf(' %u', sum(ind==i));
end
fprintf('\n')