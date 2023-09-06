

function result = list_edgeLength_tension(epithelium,balanceDict)

projection_edge = balanceDict.objectForKey('projection_edge');
p = balanceDict.objectForKey('p');
resolution = [epithelium.xResolution, ...
              epithelium.yResolution, ...
              epithelium.zResolution];

length_tension = nan(epithelium.edges.count,2);

for i = 1:epithelium.edges.count
    if projection_edge(i) < 1
        continue
    end
    edge = epithelium.edges.objectAtIndex(i);
    if edge.isBoundary
        continue
    end

    T = p(projection_edge(i));

    length_tension(i,1) = edge.edgeLength(resolution);
    length_tension(i,2) = T;
end

indices = any(isnan(length_tension),2);
length_tension(indices,:) = [];

result = length_tension;
end
