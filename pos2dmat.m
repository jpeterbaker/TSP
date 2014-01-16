function dmat = pos2dmat(pos)
% pos is an n x 2 matrix where pos(i,:) contains the x,y coordinates of node i
% dmat is a symmetric distance matrix: dmat(i,j) is the distance from node i to node j

n = size(pos,1);

a = meshgrid(1:n);
dmat = reshape(sqrt(sum((pos(a,:)-pos(a',:)).^2,2)),n,n);
