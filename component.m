function v = component(g,j0)
% g is the 0-1 adjacency matrix of a 2-regular graph
% returns a logical vector indicating which nodes are in the same component as node j0

[m,n] = size(g);

if m~=n
    error('Adjaceny matrix must be square')
end

if nargin < 2
    j0 = 1;
end

v = false(n,1);
v(j0) = 1;

ja = j0;
% Pick something the first node is connected to
jb = find(g(:,ja),1);

while jb ~= j0
    v(jb) = 1;

    % Find both nodes jb is connected to
    jc = find(g(:,jb),2);

    % The graph is 2-regular: pick the one we didn't just visit:  ja - jb - jc
    jc = jc(jc~=ja);

    ja = jb;
    jb = jc;
end

