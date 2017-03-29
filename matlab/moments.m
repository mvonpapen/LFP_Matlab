%% Calculates the first n statistical moments of time series x starting at
%% moment 2, i.e., the variance

function mom = moments ( x, n )

stdx = std (x,1);

mom = zeros (n,1);

mom(1)=mean(x);
mom(2)=var(x);
for i=3:n
    mom(i) = moment (x/stdx,i);
end