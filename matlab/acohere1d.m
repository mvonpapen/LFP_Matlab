%% Calculates time-averaged coherence from matrix C
%% alternatively calculates coherence for non-volume-conducted signals only

function [c1d, c1d_novc] = acohere1d ( f, Wx, coix, Wxy, C, phthres, uptri )


if nargin<7; uptri=0; end

[nf, nt, nx] = size ( Wx );

if nargout>1
    [nf2, nt2, nx, nx2] = size ( Wxy );

    if nf~=nf2 || nt~=nt2 || nx~=nx2
        error('Matrices Wx, Wy and Wxy do not fit!')
    end
    clear nf2 nt2 nx2
end


c1d         = NaN (nf, nx, nx);
c1d_novc    = NaN (nf, nx, nx);

if nargin<6
    phthres = 10;
end


%% Apply NaNs to matrices
for i=1:nx
    Wx(:,:,i) = coi2nan(f, Wx(:,:,i), coix(:,i));
    if nargout>1
        for j=i+1:nx
            Wxy(:,:,i,j) = coi2nan(f, Wxy(:,:,i,j), coix(:,i));
            Wxy(:,:,j,i) = Wxy(:,:,i,j);
        end
    end
end


%% Calculate Coherence
for i=1:nx
    
    for j=i+1:nx
        
        c1d(:,i,j) = abs(nanmean(Wx(:,:,i).*conj(Wx(:,:,j)),2)).^2 ./ ...
            ( nanmean(abs(Wx(:,:,i)).^2,2).*nanmean(abs(Wx(:,:,j)).^2,2) );
        if uptri~=0
            c1d(:,j,i) = c1d(:,i,j);
        end
        
        if nargout>1
            x  = Wx(:,:,i);
            y  = Wx(:,:,j);
            c  = C(:,:,i,j);
            ph = abs(angle(Wxy(:,:,i,j))/pi*180);
            vc = (ph<phthres | ph>180-phthres) & c>0.495;
            x(vc) = NaN;
            y(vc) = NaN;
            c1d_novc(:,i,j) = abs(nanmean(x.*conj(y),2)).^2 ./ ...
                ( nanmean(abs(x).^2,2).*nanmean(abs(y).^2,2) );
            if uptri~=0
                c1d_novc(:,j,i) = c1d_novc(:,i,j);
            end
        end
        
    end
    
end