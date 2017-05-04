function x = forloopbug(N1,N2,N3)

for i = N1:-1:1
    for j = N2:-1:1
        for k = 1:N3
            x(i,j,k) = i+j+k;
        end
    end
end