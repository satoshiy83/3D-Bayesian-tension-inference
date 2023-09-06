

function result = ss_mark_vertices(bitmap,vertices,indices,color)

if isnumeric(bitmap)
    bitmap = SYData(bitmap);
elseif isa(bitmap,"SYImage")
    bitmap = SYData(bitmap.drawBitmapRep(nan));
end

d = length(color);

jndices = zeros(vertices.count,1);
count = 0;
for i = 1:vertices.count
    vertex = vertices.objectAtIndex(i);
    if vertex.isBoundary
        continue
    end

    count = count + 1;
    jndices(i) = count;
end

color = permute(color(:),[3,2,1]);
for i = indices
    index = find(jndices == i);
    vertex = vertices.objectAtIndex(index);
    r = round(vertex.y);
    c = round(vertex.x);

    bitmap.var(r - 1:r + 1,c - 1:c + 1,:) = ones(3,3,d) .* color;
end

result = SYImage(bitmap);
end
