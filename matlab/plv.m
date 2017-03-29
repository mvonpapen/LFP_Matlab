%% Calculate phase-locking-value x from phases Phi

function x = plv ( Phi, coi, f, frange )

if nargin<4
    frange = [min(f) max(f)];
end
if nargin>2
    Phi = coi2nan( f, Phi, coi );
end

i   = f>=min(frange) & f<=max(frange);
Phi = Phi(i,:);

[m n] = size ( Phi );
x = NaN ( m, 1 );

for i = 1:m
    tmp = Phi ( i, ~isnan(Phi(i,:)) );
    n   = length(tmp);

    x(i) = abs( nansum( exp( sqrt(-1)*tmp ))) / n;
end