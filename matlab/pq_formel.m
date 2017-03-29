function x = pq_formel (p, q)

x(1,:) = -p/2 + sqrt( (p/2).^2 -q );
x(2,:) = -p/2 - sqrt( (p/2).^2 -q );