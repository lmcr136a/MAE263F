function [Fs, Js] = getFs(q)
global EA refLen
nv = (length(q)+1) / 4;
ne = nv - 1;
Fs = zeros( size(q) );
Js = zeros( length(q), length(q) );
for c=1:ne % Each edge
node1 = transpose(q(4*c-3:4*c-1)); % c-th node
node2 = transpose(q(4*c+1:4*c+3)); % c+1-th node
[dF, dJ] = ...
gradEs_hessEs(node1, node2, refLen(c), EA);
ind = [4*c-3, 4*c-2, 4*c-1, 4*c+1, 4*c+2,4*c+3]; % 6 numbers
Fs(ind) = Fs(ind) - dF;
Js(ind, ind) = Js(ind, ind) - dJ;
end
end
