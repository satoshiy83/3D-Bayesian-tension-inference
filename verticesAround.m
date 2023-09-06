

function result = verticesAround(vertex,edges,cells)

% Argument vertex is an index of the vertex.
% Argument edges is an array of edges.
% Argument cells is an array of cells.

% Get adjacent vertices
array = SYArray;
for i = 1:edges.count
    edge = edges.objectAtIndex(i);
    if edge(1) == vertex
        array.addObject(edge(2))
    elseif edge(2) == vertex
        array.addObject(edge(1))
    end
end
vertices = array;

% Get a vertex.
v = vertices.objectAtIndex(1);
v0 = v;

array = SYArray;
array.addObject(v);

% Get a cell containing vertex and v1.
for i = 1:cells.count
    cel = cells.objectAtIndex(i);
    if cel.containsVertex(v)
        break
    end
end

% Enumerate cells and append vertices.
while true
    for i = 1:vertices.count
        w = vertices.objectAtIndex(i);
        if w ~= v && cel.containsVertex(w)
            break
        end
    end
    if w == v0
        break
    end

    array.addObject(w);

    for i = 1:cells.count
        del = cells.objectAtIndex(i);
        if ~del.isEqual(cel) && del.containsVertex(w)
            break
        end
    end
    cel = del;
    v = w;
end

result = array;
end
