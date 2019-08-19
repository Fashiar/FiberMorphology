function [Y] = FurthestPixel(X)
%This function find the furthest two points (Start and End) of a cluster
%   X is the array of points in a luster

    D = pdist2(X,X);
    UD = unique(D);
%     M3 = maxk(UD,2);
    M2 = findmaxk(UD,2);
    [row,col] = find(D == M2(end));
    Y = [X(col(1),:);X(col(end),:)];
end

function y = findmaxk(A,k)
    A = sort(A,'descend');
    y = A(1:k);
end
