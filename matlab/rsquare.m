%% Determine R-squared for polyfit of order n

function R2 = rsquare ( x, y, n )

x = x(:);
y = y(:);

% Clear NaN's
i=~isnan(y) & ~isnan(x);
x = x(i);
y = y(i);

% Polyfit
p  = polyfit(x,y,n);
yr = polyval(p,x);
ym = mean(y);

% Sum of squares
SST = sum( ( y-ym).^2 );
SSE = sum( ( y-yr).^2 );
SSR = sum( (yr-ym).^2 );

% R-squared
R2a = SSR/SST;
R2  = 1-SSE/SST;

if round(R2a,4)~=round(R2,4)
    error('Something went wrong!')
end