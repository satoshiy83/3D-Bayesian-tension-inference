% Functon to return lengths of edges around a vertex.
% Written by Satoshi Yamashita.

function result = edgeLengthAround(center,vertices)
% Function to return lengths of edges.
% result = edgeLengthAround(center,vertices)
% Argument center is a row vector representing a vertex at the center.
% Argument vertices is an array of vertices adjacent to the center.
% Return value is an array of angles between the vertices at the center.

[n,c] = normalOfFittingPlane(vertices);

o = lf_projectOnPlane(n,c,center); % o is a projected center.

% projected is an array of projected vertices.
projected = zeros(size(vertices));
for i = 1:size(vertices,1)
    projected(i,:) = lf_projectOnPlane(n,c,vertices(i,:));
end

edges = projected - o;
edges = sqrt(sum(edges .^ 2,2));

result = edges;
end

function result = lf_projectOnPlane(n,c,p)
% Local function to project a point onto a plane.
% result = lf_projectOnPlane(n,c,p)
% Argument n is a normal of the plane.
% Argument c is a point in the plane.
% Argument p is a point to be projected.
v = p - c;
d = -dot(n,v);
result = p + n * d;
end
