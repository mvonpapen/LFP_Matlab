function c1d = cohere1d ( f, Wx, Wy, coix, coiy )
%% COHERE1D Calculates time-averaged coherence from matrix C
% alternatively calculates coherence for non-volume-conducted signals only
%   IN: ??????
%   OUT: ?????

[nf, nt, nx] = size( Wx );
[nf2, nt2, ny] = size( Wy );

if nf ~= nf2 || nt ~= nt2
    error('Matrices Wx, Wy or coix, coiy do not fit!');
end
clear nf2 nt2;

c1d = NaN(nf, nx, ny);


%% Apply NaNs to matrices
if nargin>3
    for i=1:nx
        Wx(:,:,i) = coi2nan(f, Wx(:,:,i), coix(:,i));
    end
    for j=1:ny
        Wy(:,:,j) = coi2nan(f, Wy(:,:,j), coiy(:,j));
    end
end


%% Calculate Coherence
for i = 1:nx    
    for j = 1:ny        
        c1d(:,i,j) = abs(nanmean(Wx(:,:,i).*conj(Wy(:,:,j)),2)).^2 ./ ...
            ( nanmean(abs(Wx(:,:,i)).^2,2).*nanmean(abs(Wy(:,:,j)).^2,2) );
    end    
end