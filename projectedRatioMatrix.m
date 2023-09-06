

function result = projectedRatioMatrix(epithelium)

%% Inner parameter
threshold_ratio = 0.9;

%% Implementation
projection_edge = zeros(epithelium.edges.count,1);
count_m = 0;
for i = 1:epithelium.edges.count
    edge = epithelium.edges.objectAtIndex(i);
    if ~edge.isBoundary
        count_m = count_m + 1;
        projection_edge(i) = count_m;
    end
end

count_n = 0;
for i = 1:epithelium.vertices.count
    vertex = epithelium.vertices.objectAtIndex(i);
    if ~vertex.isBoundary
        count_n = count_n + 1;
    end
end
A = zeros(count_n,count_m);

count = 1;
resolution = [epithelium.xResolution, ...
              epithelium.yResolution, ...
              epithelium.zResolution];
array = SYArray;
for i = 1:epithelium.vertices.count
    vertex = epithelium.vertices.objectAtIndex(i);
    if vertex.isBoundary
        continue
    end

    center = [vertex.x,vertex.y,vertex.z] .* resolution;
    around = [];
    dict = epithelium.verticesAround(vertex);
    vertices = dict.objectForKey('vertices');
    for j = 1:vertices.count
        wertex = vertices.objectAtIndex(j);
        around = cat(1,around,[wertex.x,wertex.y,wertex.z] .* resolution);
    end

    ratio = projectionRatioAround(center,around);

    if any(ratio < threshold_ratio)
        disp('It found an edge angled to projected plane.')
        eict = SYDictionary;
        eict.setObjectForKey('vertex',vertex);
        eict.setObjectForKey('center',center);
        eict.setObjectForKey('around',around);
        array.addObject(eict);
    end

    edges = dict.objectForKey('edges');
    for j = 1:edges.count
        edge = edges.objectAtIndex(j);
        A(count,projection_edge(edge.index)) = ratio(j);
    end

    count = count + 1;
end

mean_dev = zeros(2,count_m);
for i = 1:count_m
    ratios = A(:,i);
    ratios = ratios(ratios > 0);
    if length(ratios) > 1
        mean_dev(1,i) = mean(ratios);
        mean_dev(2,i) = std(ratios);
    else
        mean_dev(1,i) = ratios;
    end
end

dict = SYDictionary;
dict.setObjectForKey('ratioMatrix',A);
dict.setObjectForKey('projection_edge',projection_edge);
dict.setObjectForKey('steep_edges',array);
dict.setObjectForKey('mean_dev',mean_dev);

result = dict;
end
