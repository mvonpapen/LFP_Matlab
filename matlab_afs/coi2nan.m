function waveOut=coi2nan(f,waveIn,coi)
%%COI2NAN Set data outside coi to NaN
% waveOut=coi2nan(f,waveIn,coi)
%
%   IN: 
%       f:      frequency vector
%       waveIn: Wavelet coefficient Matrix
%       coi:    cone of influence vector
%   OUT:
%       waveOut: Wavelet coefficient Matrix with NaNs
%
%   Author: Michael von Papen, Date: 13.10.15

[a, b, c] = size( waveIn );

if nargin < 3; coi = gencoi(b); end

[b2, c2] = size( coi );

if b ~= b2
    coi = coi';
    [b2 c2] = size( coi );
    if b ~= b2
        error('Number of data points don''t match!')
    end
end

if c2 == 1
    C = repmat(coi(:)',[a, 1, c]);
else
    error('Only single channels can be processed, sorry.') %TODO????
%     C=repmat(coi,[a 1 1]);
end

warning('off') % switch off warnings because of divide by zero
P = repmat(1./f(:),[1, b, c]);
warning('on')
waveOut = waveIn;

waveOut(C<P) = NaN;