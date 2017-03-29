function [ data, W, coi, Pw, ff, Pf ] = procdata( data, varargin )
%% PROCDATA Process data for spectral analysis:
%  1. Filter time series if desired
%  2. Compute Fourier transform
%  3. Compute wavelet transform
%  4. Compute PSD of both methods
%
%   INPUT:
%           data:    signal of k channels (nt x nch elements)
%
%   OPTIONAL INPUT:
%           dt:     Time increment
%           freq:   Frequency vector in Hz (nf elements) that shall be 
%                   evaluated
%           order:  filter order
%           filter: Structure with high and stop band filter frequencies
%           rect:   1=rectify data
%           inc:    1=increment data
%           art:    Cell containing the artefact-intervals in data points
%
%  	OUTPUT:
%           data:   Processed data (n elements)
%           W:      Wavelet transform of data (nf x nt x nch elements)
%           coi:    Wavelet cone of influence
%           Pw:     Global wavelet spectrum (nf x nch elements)
%           ff:     Frequency vector for Fourier spectrum
%           Pf:     Fourier-based power spectrum (m x k elements)
% 
% Author: Michael von Papen
% 
% Date: 14.10.15


args = struct('filter',[],... % for stop-band e.g. [45 55; 90 110]
            'order',5,...
            'freq',logspace(0,2.7,30),...
            'rect',0,...
            'dt',1/2456,...
            'art',[],...
            'inc',0,...
            'w0',6, ...
	        'dec', 1);
args = parseArgs(varargin,args);
arts = args.art;
w0 = args.w0;
dec = args.dec;


%% Parameters
[n, nch] = size(data);
if nch>n
    data = data';
    [n, nch] = size(data);
end
f = args.freq;
m = length(f);
dt = args.dt;


%% Preset variables
if nargout > 1
    W = NaN(m,n,nch);
    Pw = NaN(m,nch);
end

if nch == 1 && ~iscell(arts)
    arts = {arts};
end


%% Filter time series with two-side butterworth filter order 5
if args.order ~= 0 && ~any(isnan(data(:))) && ~isempty(args.filter);
    F = args.filter;
    for i = 1:nch
        for j = 1:size(F,1)
            data (:,i) = bfilt ( data(:,i), F(j,:)*2*dt, args.order );
        end
    end
end
if args.rect == 1
    data = hilbert(data);
    data = abs(data);
end


%% Prewhiten data if desired
if args.inc == 1
    data = data(2:end,:)-data(1:end-1,:);
    n = n-1;
end


%% Decimate data if desired
if dec ~= 1
    for i = 1:nch
      tmp(:,i) = decimate( data(:,i),dec );
    end
    data = tmp;
    n = length(data);
    dt = dt*dec;
end


%% Spectral Transormation
% Fourier Spectrum
if nargout > 4 && ~any(isnan(data(:)));
    [Pf,ff] = psdf( data, args.dt );
end


% Wavelet Transform
if nargout > 1
    W = NaN(m,n,nch);
    Pw = NaN(m,nch);
    coi = NaN(n,nch);
    for i = 1:nch
        [W(:,:,i),~,~,coi(:,i)] = waveletlin(data(:,i),dt,f,1,'MORLET',w0);
        if ~isempty(arts)
            coi(:,i) = artcoi(coi(:,i),arts{i});
        end
        Pw(:,i) = psdw(squeeze(W(:,:,i)),coi(:,i),f);
    end
end

%% Postdarken spectra if data was prewhitened
if args.inc == 1
    if nargout > 1
        W(:,end+1,:) = NaN;
        coi(end+1,:) = NaN;
        for j = 1:nch
            Pw(:,j) = Pw(:,j)./(4*sin(pi*dt*f).^2)';
            parfor i = 1:n
                W(:,i,j) = W(:,i,j)./(4*sin(pi*dt*f))';
            end
        end
    end
    if nargout > 4
        for i = 1:nch
            Pf(:,i) = Pf(:,i)./(4*sin(pi*dt*ff).^2)';
        end
    end
    data(end+1,:) = NaN;
end