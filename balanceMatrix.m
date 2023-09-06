% Function to return a balance matrix for a geometrical tension inference
% with an epithelium.
% Written by Satoshi Yamashita.

function result = balanceMatrix(epithelium)
% Function to return a balance matrix from an epithelium.
% result = balanceMatrix(epithelium)
% Argument epithelium is an SAEpithelium instance.
% Return value is an SYDictionary instance containing:
% % "balanceMatrix": the balance matrix
% % "projection_edge": indexing vector from index in epithelium.edges to
% index in the balance matrix column.
% % "projection_cell": indexing vector from index in epithelium.cells to
% index in the balance matrix column.

projection_edge = zeros(epithelium.edges.count,1);
count_m = 0;
for i = 1:epithelium.edges.count
    edge = epithelium.edges.objectAtIndex(i);
    if ~edge.isBoundary
        count_m = count_m + 1;
        projection_edge(i) = count_m;
    end
end
projection_cell = zeros(epithelium.cells.count,1);
for i = 1:epithelium.cells.count
    cel = epithelium.cells.objectAtIndex(i);
    if ~cel.isBoundary
        count_m = count_m + 1;
        projection_cell(i) = count_m;
    end
end

count_n = 0;
for i = 1:epithelium.vertices.count
    vertex = epithelium.vertices.objectAtIndex(i);
    if ~vertex.isBoundary
        count_n = count_n + 1;
    end
end
A = zeros(count_n * 2,count_m);

count = 1:2;
resolution = [epithelium.xResolution, ...
              epithelium.yResolution, ...
              epithelium.zResolution];
for i = 1:epithelium.vertices.count
    vertex = epithelium.vertices.objectAtIndex(i);
    if vertex.isBoundary
        continue
    end

    center = [vertex.x,vertex.y,vertex.z] .* resolution;
    around = [];
    dict = epithelium.verticesAround(vertex);
    vertices = dict.objectForKey("vertices");
    for j = 1:vertices.count
        wertex = vertices.objectAtIndex(j);
        around = cat(1,around,[wertex.x,wertex.y,wertex.z] .* resolution);
    end

    [balance_t,balance_p] = balanceAround(center,around);

    edges = dict.objectForKey("edges");
    for j = 1:edges.count
        edge = edges.objectAtIndex(j);
        A(count,projection_edge(edge.index)) = balance_t(:,j);
    end

    cells = dict.objectForKey("cells");
    for j = 1:cells.count
        cel = cells.objectAtIndex(j);
        A(count,projection_cell(cel.label)) = balance_p(:,j);
    end

    count = count + 2;
end

dict = SYDictionary;
dict.setObjectForKey("balanceMatrix",A);
dict.setObjectForKey("projection_edge",projection_edge);
dict.setObjectForKey("projection_cell",projection_cell);

result = dict;
end
