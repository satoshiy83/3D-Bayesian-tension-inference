

function result = ss_tricellulr_junctions(segmented)

if isnumeric(segmented)
    segmented = SYImage(SYData(segmented));
end

% Label cells
labeled = IPConnectedComponents.connectedBinaryComponents(segmented,4);

% Label vertices
vertices = SYArray;

siz = segmented.frameSize;
bitmapRep = segmented.representations.objectAtIndex(1);
bitmap = bitmapRep.bitmap;
indices = (find(bitmap.var == 0))';

for i = indices
    [r,c] = ind2sub(siz,i);
    if r > 1 && r < siz(1) && c > 1 && c < siz(2)
        scope = labeled(r - 1:r + 1,c - 1:c + 1);
        inci = unique(scope(:));
        if length(inci) > 3
            vertex = SAVertex;
            vertex.x = double(c);
            vertex.y = double(r);
            inci(inci == 0) = [];
            vertex.incidentCells = inci;

            vertices.addObject(vertex);
        end
    end
end

result = vertices;
end
