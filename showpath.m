function showpath(pos,path)
% pos is an n x 2 matrix where pos(i,:) contains the x,y coordinates of node i
% path(i) is the index of the ith node to visit
%
% This function plots the points in pos connected in the order indicated by path

plot(pos(path,1),pos(path,2),'o-');

