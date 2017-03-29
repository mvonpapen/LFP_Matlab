function newCoi = artcoi ( coi, arts)
%% ARTCOI Updates cone of interest so that artefacts are left out.
%       Only considers Wavelet being Morlet (param=6). 
%       newCoi = ARTCOI( coi, arts )
%       OUT:
%           newCoi: new COI, same length as input coi
%       IN:
%           coi:    COI from Wavelet-transformation
%           arts:   array of intervals of artifacts in time series 
%                   (dimensions: number of artifacts x 2 (t_min, t_max) )
%
% Author: Michael von Papen, Date: 13.10.15

newCoi = coi(:)';
if any(isnan(arts));
    warning('Some artifact intervals contain NaN.');
    return
end
if isempty(arts)
%      warning('No artifacts found for artcoi.');
    return
end
fs = 2500; %sampling frequency
t = (0:length(coi)-1)/fs;

[n,m] = size (arts);

if m ~= 2
    arts = arts';
    temp = n;
    n = m;
    m = temp;
end

if m ~= 2
    error('Size of artefact vector not OK!');
end


for j = 1:n
    if arts(j,1) > length(t)
        continue
    end
    i = arts(j,1):min([arts(j,2) length(t)]);
    % sqrt(2)/1.03./(t-t(i(1))) according to Torrence & Compo(1998), 
    % but 2/1.03 more conservative
    newCoi = 1./max([1./newCoi;...
        sqrt(2)/1.03./-(t-t(i(1)));...
        sqrt(2)/1.03./(t-t(i(end)))]);
    newCoi(i) = 1e-99;
end