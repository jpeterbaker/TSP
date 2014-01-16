function [path,bestDist] = tsp(dmat,type,start,finish)
% This solves the symmetric traveling salesman problem
%
% dmat is a symmetric distance matrix
% type is 0 or 1: 0 for round trip (default)
%                 1 for one-way
% start/finish are optional starting starting/ending points (for one-way case)
% start may be included without finish
%
% path is a vector of point indices giving the shortest tour
%
%
% Uses branch and bound to solve the edge inclusion problem:
% min  z'*c
% s.t. z(i) = 0 or 1
%      z represents edges in a tour

if nargin < 2
    type = 0;
end

if type == 1
    % one-way trip

    % Add a node (with universal free moves)
    dmat(end+1,end+1) = 0;

    if nargin >= 3
        % sum of all edge lengths as a proxy for infinity 
        toomuch = sum(sum(dmat))/2;

        % add a node with free move to start (equal high cost to all other nodes)
        dmat(end,:) = toomuch;
        dmat(:,end) = toomuch;
        dmat(end,[start end]) = 0;
        dmat([start end],end) = 0;

        % give free move to finish as well
        if nargin == 4
            dmat(end,[finish end]) = 0;
            dmat([finish end],end) = 0;
        end
            
    end
    % Else we want the shortest of all one-way trips: the universal teleport node is all we need
end

n = size(dmat,1);
TRIn = n*(n-1)/2;

% A distance vector
c = squareform(dmat);

% Give the linear relaxation some slack on coming up with exactly 0 or 1
tol = 1e-5;

[greedyPath,greedyDist] = greedyTSP(dmat);
bestDist = greedyDist;
%disp(greedyPath)
%disp(bestDist)

options=optimset('Display', 'off');
%options=optimset();

% These are constraints for layer's problem
% Initializing them here lets them persist through all layers
% This persistence grants an enormous speedup
A = [];
b = [];

[z,bestDist] = layer();

% If the greedy solution was optimal, all layers would have terminated early
if bestDist >= greedyDist
    path = greedyPath;
    bestDist = greedyDist;
else
    zmat = squareform(z);
    path = cycleOrder(zmat);
end

% In the one-way case, make path start in correct place
if type == 1
    % remove the loop-closing move
    path = path(1:n);

    % find and remove node "n" (the teleportation node)
    nloc = find(path==n,1);
    path = path([nloc+1:n 1:nloc-1]);

    if nargin > 2
        % make sure we start at the start instead of the end
        if path(1) ~= start
            path=flipud(path);
        end
    end
end

function inds = edgeInds(ii,jj)
% Swap so that i is always larger
[ii,jj] = deal(max(ii,jj),min(ii,jj));

% This is a tricky expression to derive
inds = -n + jj.*(2*n-1-jj)/2 + ii;
end

function [z,dist] = layer(Aeq,beq,LB,UB,depth)
%This is the heart of the algorithm
% Solves:
%    min z'*c
%    s.t. A*z <= b
%         Aeq*z = beq
% If the minimum is greater than the best integral solution in any branch, this branch fails
% Otherwise, if the solution is not integral, modify A and b in two different ways (enforcing 0 and then 1 on some non integer entry) and compare the branches
% If the solution is integral, make sure it represents a tour. If not, add a constraint to join the subtours

% Start with only 0<=z<=1, 2-regularity as constraints
if nargin < 1
    depth = 0;

    Aeq = reg2(n);
    beq = 2*ones(n,1);

    LB = zeros(TRIn,1);
    UB = ones(TRIn,1);

end

if depth < 15
%    disp(depth);
end

[z,dist,flag] = linprog(c,A,b,Aeq,beq,LB,UB,[],options);

if flag ~= 1 && flag ~= -2
    disp '    linprog had trouble'
    disp '    flag'
    disp(flag)
    disp '    equality constraints:'
    disp(size(Aeq,1))
    disp '    inequality constraints:'
    disp(size(A,1))
    dist = inf;
    return
end

if dist >= bestDist
%    disp 'Linear relaxation couldn''t beat best integer solution'
    dist = inf;
    return
end

err = min(abs(z-1),z);

if all(err<tol)
    % We have an integer solution
%    disp 'Integer solution'

    % Eliminate numerical error
    logone = z>.5;
    z(logone)  = 1;
    z(~logone) = 0;

    % Make sure it represents a tour
    comp = component(squareform(z));

    if all(comp)
        % This is a feasible answer
%        disp 'It''s a tour'
        if dist < bestDist
            bestDist = dist;
%            disp(bestDist)
%            disp(cycleOrder(squareform(z)))
        end
        return
    end

%    disp 'It''s not a tour'

    % We need to try again and force the disjoint components to have at least 2 edges between them

    % Identify all the edges between the components
    ii = find(comp);
    jj = find(1-comp);

    [ii,jj] = meshgrid(ii,jj);

    between = edgeInds(ii,jj);
    
    newrow = zeros(1,TRIn);
    newrow(between) = -1;

    % Add the constraint
    A = [A;newrow];
    b = [b;-2];

    [z,dist] = layer(Aeq,beq,LB,UB,depth+1);
else
    % We should consider 0 and 1 branches of the least feasible entry
    [~,i] = max(err);

    % 0 case
    UB(i) = 0;
    [z0,dist0] = layer(Aeq,beq,LB,UB,depth+1);
    UB(i) = 1;
    
    % 1 case
    LB(i) = 1;
    [z1,dist1] = layer(Aeq,beq,LB,UB,depth+1);
    LB(i) = 0;

    % Pick the better one
    if dist0 < dist1
        dist = dist0;
        z    = z0;
    else
        dist = dist1;
        z    = z1;
    end

end

end

function B = reg2(n)
% Returns B such that B*z is the node-degree vector (where z is the adjacency vector)

% Allocate the space
B = zeros(n,TRIn);

rangen = 1:n;

% Consider all i,j combinations
for j = 1:n
    
    % i and j should never match (ii and jj are vectors)
    ii = [rangen(1:j-1) rangen(j+1:n)];

    % Vector-indices of all edges connected to j
    jedges = edgeInds(ii,j);

    % Mark all jedges in B
    B(j,jedges) = 1;

end

end

function cycle = cycleOrder(g)
% g is the adjacency matrix of a cycle
% returns a vector of the node indices in order around the cycle

j0 = 1;

cycle = nan(n+1,1);
cycle(1) = 1;

ja = j0;
% Pick something the first node is connected to
jb = find(g(:,ja),1);

for i=2:n
    cycle(i) = jb;

    % Find both nodes jb is connected to
    jc = find(g(:,jb),2);

    % The graph is 2-regular: pick the one we didn't just visit:  ja -> jb -> jc
    jc = jc(jc~=ja);

    ja = jb;
    jb = jc;
end

cycle(n+1) = cycle(1);

end

function [tourOrder,tourDist] = greedyTSP(dmat,j)
% dmat is a distance matrix
% j is the index of the starting node (j=1 by default)
%
% tourOrder(i) is the index of the ith node to visit
% tourDist is the distance traveled in the tour
%
% Produces a greedy solution to the TSP
% Starts at node j and always chooses the closest unvisited node to visit next

[m,n] = size(dmat);
if m ~= n
    error('Distance matrix must be square');
end

if nargin < 2
    j = 1;
end

tourOrder    = nan(n+1,1);
tourOrder(1) = j;
visited      = false(n,1);

for i = 1:n
    [~,distOrder] = sort(dmat(j,:));

    for next = distOrder
        if ~visited(next)
            visited(next) = 1;
            break
        end
    end
    tourOrder(i) = next;
end

tourOrder(end) = tourOrder(1);

if nargout > 1
    tourDist = 0;
    for i = 1:n
         tourDist = tourDist + dmat(tourOrder(i),tourOrder(i+1));
    end
end

end

end

