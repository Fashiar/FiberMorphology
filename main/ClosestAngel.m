function [I,dir] = ClosestAngel(B,tol)
    dir = B(1);
    I = find(B >= dir-tol & B <= dir+tol);
end