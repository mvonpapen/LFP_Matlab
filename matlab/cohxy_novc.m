%% Eliminates X-X volume conduction effects from X-Y 1d-coherences
%%
%% 

function [cxy_novc, cxy] = cohxy_novc ( Wxx, Cxx, Cxy, f, coix, coiy, phthres )


if nargin<7
    phthres=10;
end


%% Check size of matrix inputs

[nf, nt, nx, nx2] = size(Wxx);
if nx~=nx2
    error('Wxx not symmetric!')
end

[nf, nt, nx2, ny] = size(Cxy);
if nx~=nx2
    error('Wxx and Cxy do not fit!')
end



%% Preset variables
cxy      = NaN (nf, nx, ny);
cxy_novc = NaN (nf, nx, ny);



%% Caluclate 1d coherences

for i=1:nx
    
    for j=1:ny
        
        if nargout>1
            % Simple 1d coherence = temporal mean of coherence matrix
            cxy(:,i,j) = nanmean(coi2nan(f,Cxy(:,:,i,j),min([coix(:,i)'; coiy(:,j)'])),2);
        end
        
        
        %Begin loop over x to determine X-X volume conduction
        c_novc = NaN (nf, nx);
        
        for k=setdiff(1:nx,i)
            % Determine indices where phase near zero and coherence
            % significant
            c  = Cxx(:,:,i,k);
            ph = abs(angle(Wxx(:,:,i,k))/pi*180); 
            vc = (ph<phthres | ph>180-phthres) & c>0.495;
            % Set volcond-indices (vc) to zero and calculate vc-free
            % coherence c_novc
            c  = Cxy(:,:,i,j);
            c(vc) = 0;
            c_novc(:,k) = nanmean(coi2nan(f,c,min([coix(:,i)'; coiy(:,j)'; coix(:,k)'])),2);
        end
        
        % Average c_novc equivalent to subtraction of mean volume
        % conduction
        cxy_novc(:,i,j) = nanmean(c_novc,2);
                
    end
    
end