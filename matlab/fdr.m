% Compute false discovery rate according to
%
% "False discovery rate control is a recommended
% alternative to Bonferroni-type adjustments in health studies"
%   by Glickman et al., 2014, Journal of Clinical Epidemiology
%
%
% Author: Mitch von Papen
% Date: 19.04.2016

function id = fdr ( p, pval )

if nargin < 2
    pval = 0.05;
end

N = numel(p);

thresh = (1:N)'/N*pval;

[p, id] = sort( p(:) );

i = find(p>thresh, 1, 'first')-1;
if i==0
    id = [];
else
    id = id(1:i);
end