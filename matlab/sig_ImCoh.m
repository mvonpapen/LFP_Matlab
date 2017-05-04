%% Significance of Imaginary Part of Coherency 
% Function determines which entries of Coherency matrix are significant at 
% a pval-level
%
% Author: M. von Papen
% Date:   May 03, 2017
%
% INPUT:
%   Coherency: Nf x Nt x Nch x Nch matrix of complex valued coherency
%   pval: p-value used to determine significance threshold, can be either p=0.01 or p=0.05 [pval=0.01]
%   nsig: averaging width over which coherency was determined, nsig=n_cy/n_co from Lachaux et al., 2002
%
% OUTPUT:
%   ICthreshd: thresholded imaginary part of coherency, ICthreshd=1 if significant, else 0
%   StdC: standard deviation for each matrix entry according to Nolte et al., 2004

function [ICthreshd, StdC] = sig_ImCoh (Coherency, pval, nsig)


if nargin<2
    pval = 0.01;
end

switch pval
    case 0.01
        std_fac = 2.575829;
    case 0.05
        std_fac = 1.959964;
    case 0.1
        std_fac = 1.644854;
end

StdC        = sqrt( (1-abs(Coherency).^2) ./ 2/(1+2*nsig) );
ICthreshd   = zeros(size(StdC));
ICthreshd( abs(imag(Coherency))>std_fac*StdC ) = 1;