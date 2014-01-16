
% This file demonstrates use cases of tsp.m

% Read the contents of the file
filecontents = importdata('EXAMPLE_LOCATIONS.txt',',');

xy = filecontents.data;
names = filecontents.textdata;

% Flip so that longitude is x, latitude is y
xy = fliplr(xy);

dmat = pos2dmat(xy);

% [path,bestDist] = tsp(dmat,1,12,29); % Shortest path with start and finish specified
% [path,bestDist] = tsp(dmat,1,12);    % Shortest path with only the start specified
% [path,bestDist] = tsp(dmat,1);       % Shortest path
[path,bestDist] = tsp(dmat);           % Shortest loop

disp(sprintf('The following tour has length %d:',bestDist));
disp(names(path));

showpath(xy,path);

