%% Analyze coherent PSD
i0   = [1  6 11 17 25 31 36 41 50];
iend = [5 10 16 24 30 34 40 49 56];


for n=1:length(i0)
    for i=i0(n):iend(n)
        if nEMG(i)==0
            continue
        end
        [m1 m2 m3] = size (WLFP{i});
        for j=1:nEMG(i)
%             tmp=zeros(m1,m2,m3);
%             c=find(CLE{i}(:,:,:,j)>=0.5);
%             tmp(c)=2/2500*abs(WLFP{i}(c)).^2;
            for k=1:m3
%                 PP{n}(:,i-i0(n)+1,k,j)=mygws(tmp(:,:,k), dt, ...
%                     min([coiLFP{i}(:,k)';coiEMG{i}(:,j)']), f);
                Pn{n}(:,i-i0(n)+1,k,j)=PP{n}(:,i-i0(n)+1,k,j)./PP{n}(:,1,k,j);
%                 D{n}(i-i0(n)+1)=TT(i);
%                 R{n}(i-i0(n)+1)=rigor(i);
            end
        end
    end
end