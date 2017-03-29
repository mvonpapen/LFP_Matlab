%% Removes large deviation from baseline, such as a Spike, from data by
%% fitting a 3rd order polynomial and subtracting it. This works in
%% principle like a high pass filter
%%
%%  Input:
%%      x   data vector
%%
%%  Ourput:
%%      xc  corrected data vector
%%
%%  Author: Michael von Papen
%%  Date:   04.12.2014

function    xc = remspk ( x, varargin )

args=struct('order',3, 'interpol', 0);
args=parseArgs(varargin,args);


if args.interpol == 0
    
    %% Normalize data to allow for good fitting
    [z,m,v]=zscore(x(:)');
    t=[1:length(z)];

    %% Evaluate polynomial
    p = polyfit ( t, z, args.order );
    xc = ( z - polyval( p, t ) ) .* v;

    %% Refit corrected data to original one by shifting along y-axis
    d = x(1)-xc(1);
%     c=mean(d([1 end]));
    xc=xc+d;
    
else
    
    %% if polynomial doesnt work, interpolate
    xc=linspace(x(1),x(end),length(x));
    
end