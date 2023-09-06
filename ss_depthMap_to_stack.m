

function result = ss_depthMap_to_stack(depthMap,vertices,depth)

bitmapRep = depthMap.representations.objectAtIndex(1);
bitmapD = bitmapRep.bitmap;

siz = size(bitmapD.var);

stack = zeros(siz(1),siz(2),3,depth,'uint8');

indices = (find(bitmapD.var > 0))';
for i = indices
    [r,c] = ind2sub(siz,i);
    I = depth - bitmapD.var(r,c);
    stack(r,c,2,I) = 255;
end

for i = 1:vertices.count
    vertex = vertices.objectAtIndex(i);
    r = vertex.y;
    c = vertex.x;
    if bitmapD.var(r,c) > 0
        I = depth - bitmapD.var(r,c);
        stack(r - 1:r + 1,c - 1:c + 1,1,I - 1:I + 1) = 0;
        stack(r - 1:r + 1,c - 1:c + 1,2,I - 1:I + 1) = 255;
        stack(r - 1:r + 1,c - 1:c + 1,3,I - 1:I + 1) = 255;
    else
        disp(['Z-position was not assigned to ',num2str(i),'-th vertex.']);
    end
end

image = SYImage(SYData(stack(:,:,:,1)));
for i = 2:depth
    bitmapRep = SYBitmapImageRep(SYData(stack(:,:,:,i)));
    image.addRepresentation(bitmapRep);
end

result = image;
end
