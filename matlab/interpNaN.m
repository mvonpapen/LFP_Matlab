%% Interpolates over NaNs and gives out interpolated data together with
%% interval positions of NaNs to feed, e.g., artcoi

function [data, arts] = interpNaN ( data )

[n, nch] = size(data);

ti=cell(nch,1);
arts=cell(nch,1);

%% Determine NaN intervals
for i=1:nch
    ti{i}=find(isnan(data(:,i)));
    if ~isempty(ti{i})
        %are there multiple NaN intervals?
        jump=find(diff(ti{i})~=1);
        if isempty(jump)
            arts{i}=[ti{i}(1) ti{i}(end)];
        else
            arts{i}=[ti{i}([1 jump(:)'+1]) ti{i}([jump(:)' end])];
        end
    end
end

for i=1:nch
    if ~isempty(arts{i})
        [m1, m2]=size(arts{i});
        for j=1:m1
            if arts{i}(j,1)==1
                data(arts{i}(j,1):arts{i}(j,end)+1,i)=...
                    data(arts{i}(j,end)+1,i);
            else
                data(arts{i}(j,1)-1:arts{i}(j,end)+1,i)=...
                    interp1([arts{i}(j,1)-1 arts{i}(j,end)+1],...
                        data([arts{i}(j,1)-1 arts{i}(j,end)+1]),...
                        arts{i}(j,1)-1:arts{i}(j,end)+1);
            end
        end
    end
end