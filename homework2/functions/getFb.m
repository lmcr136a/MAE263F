function [Fb, Jb] = getFb(q, m1, m2)
global kappaBar EI voronoiLength
nv = (length(q)+1) / 4;
Fb = zeros( size(q) );
Jb = zeros( length(q), length(q) );
for c=2:nv-1 % Compute bending force at each internel node
node0 = transpose(q(4*c-7:4*c-5));
node1 = transpose(q(4*c-3:4*c-1));
node2 = transpose(q(4*c+1:4*c+3));
m1e = m1(c-1, :); % m1 vector of c-1-th edge
m2e = m2(c-1, :); % m2 vector of c-1-th edge
m1f = m1(c, :); % m1 vector of c-th edge
m2f = m2(c, :); % m2 vector of c-th edge
[dF, dJ] = ...
gradEb_hessEb(node0, node1, node2, ...
m1e, m2e, m1f, m2f, ...
kappaBar(c,:), voronoiLength(c), EI);
ind = 4*c-7:4*c+3; % 11 numbers
Fb(ind) = Fb(ind) - dF;
Jb(ind, ind) = Jb(ind, ind) - dJ;
end
end
