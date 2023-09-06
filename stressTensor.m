

function result = stressTensor(epithelium,balanceDict,cells)

projection_edge = balanceDict.objectForKey("projection_edge");
projection_cell = balanceDict.objectForKey("projection_cell");
p = balanceDict.objectForKey("p");

resolution = [epithelium.xResolution, ...
              epithelium.yResolution, ...
              epithelium.zResolution];

totalArea = 0;

pressure = zeros(3);
tension = zeros(3);

edges = SYSet;
for i = cells
    cel = epithelium.cells.objectAtIndex(i);
    edges.addObjectsFromArray(cel.edges);

    P = p(projection_cell(i));
    A = cel.apicalArea(resolution);
    pressure = pressure + eye(3) * P * A;

    totalArea = totalArea + A;
end

edges = edges.allObjects;
for i = 1:edges.count
    edge = edges.objectAtIndex(i);
    index = epithelium.edges.indexOfObject(edge);
    vertex = edge.incidentVertices.objectAtIndex(1);
    wertex = edge.incidentVertices.objectAtIndex(2);
    l = [vertex.x,vertex.y,vertex.z] - [wertex.x,wertex.y,wertex.z];
    l = l .* resolution;

    T = p(projection_edge(index));
    tension = tension + l' * l * T / sqrt(sum(l .^ 2));
end

result = (tension - pressure) / totalArea;
end
