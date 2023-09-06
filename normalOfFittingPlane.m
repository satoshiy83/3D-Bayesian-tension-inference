% Function to return a normal vector and the closest point to the origin on
% a plane fitted to given points by the least square method.
% Written by Satoshi Yamashita.

function [result1,result2] = normalOfFittingPlane(xyz)
% Function to return a plane by least square fitting for points.
% [result1,result2] = normalOfFittingPlane(xyz)
% Argument xyz is a matrix whose columns represent x, y, and z coordinates
% of points.
% Return value result1 is a normal vector of the fitten plane with length
% 1.
% Return value result2 is a point closest to the origin.

% B is a plane obtained by the least square method for points xyz.
B = [xyz(:,1:2),ones(size(xyz,1),1)] \ xyz(:,3);

p = [1,0,B(1)]; % p is a vevtor on the plane B
q = [0,1,B(2)]; % q is another vector on the plane B
n = [p(2) * q(3) - p(3) * q(2), ...
     p(3) * q(1) - p(1) * q(3), ...
     p(1) * q(2) - p(2) * q(1)];
result1 = n / sqrt(sum(n .^ 2));

o = [0,0,-B(3)]; % o is a vector to the origin from a point on the plane B
d = -dot(result1,o); % d is a distance between the origin and the plane B
result2 = result1 * d;
end
