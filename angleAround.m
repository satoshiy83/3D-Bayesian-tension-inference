% Function to return angles between edges around a vertex.
% Written by Satoshi Yamashita.

function result = angleAround(center,vertices)
% Function to return angles between incident edges.
% result = angleAround(center,vertices)
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
% normalized is an array of the projected vertices.
normalized = edges ./ sqrt(sum(edges .^ 2,2));

% angles is the angle between the vertices.
% angles(i) represents an angle by projected(i), o, projected(i + 1).
angles = zeros(size(vertices,1),1);
for i = 1:size(vertices,1) - 1
    angles(i) = acos(dot(normalized(i,:),normalized(i + 1,:)));
end
angles(end) = acos(dot(normalized(end,:),normalized(1,:)));

result = angles;
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
