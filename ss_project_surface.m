% Short script to draw depth map of an epithelium surface.
% Written by Satoshi Yamashita.

function result = ss_project_surface(segmented,epithelium)
% Function to draw z depth map of an epithelium surface.
% result = ss_project_surface(segmented,epithelium)
% Argument segmented is an image or name of the image file, and the image 
%   shows cells by none zero value connected pixels separated by 0 value
%   pixels.
% Argument epithelium is an SAEpithelium instance.
% Return value is a grayscale image where brightness indicates the depth.

if isstring(segmented)
    segmented = SYImage(segmented);
elseif isnumeric(segmented)
    segmented = SYImage(SYData(segmented));
end

labeled = IPConnectedComponents.connectedBinaryComponents(segmented,4);

bitmap = zeros(segmented.frameSize);

for i = 1:epithelium.cells.count
    cel = epithelium.cells.objectAtIndex(i);
    if cel.isBoundary
        continue
    end

    mask = labeled == cel.label;

    xyz = [];
    for j = 1:cel.vertices.count
        vertex = cel.vertices.objectAtIndex(j);
        xyz = cat(1,xyz,[vertex.x,vertex.y,vertex.z]);
    end
    B = [xyz(:,1:2),ones(cel.vertices.count,1)] \ xyz(:,3);

    indices = find(mask);
    [r,c] = ind2sub(segmented.frameSize,indices);
    z = [c,r,ones(length(c),1)] * B;
    bitmap(indices) = z;
end

for i = 1:epithelium.edges.count
    edge = epithelium.edges.objectAtIndex(i);
    
    v1 = edge.incidentVertices.objectAtIndex(1);
    v2 = edge.incidentVertices.objectAtIndex(2);
    e_xy = [v2.x - v1.x; v2.y - v1.y];
    l = sqrt(sum(e_xy .^ 2));
    e_xy = e_xy / l;
    dz_xy = (v2.z - v1.z) / l;

    p = double([edge.x,edge.y]) - [v1.x,v1.y];
    z = dz_xy * (p * e_xy) + v1.z;

    indices = sub2ind(segmented.frameSize,edge.y,edge.x);
    bitmap(indices) = z;
end

result = SYImage(SYData(bitmap));
end
