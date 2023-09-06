

function result = draw_tension(epithelium,balanceDict)

frameSize = epithelium.frameSize;
bitmap = nan(frameSize);

projection_edge = balanceDict.objectForKey("projection_edge");
p = balanceDict.objectForKey("p");

mask = SYData(zeros(frameSize,'logical'));

if ~balanceDict.isNanForKey("range_pressure")
    range = balanceDict.objectForKey("range_pressure");
else
    n_edge = sum(projection_edge > 0);
    tensions = p(1:n_edge);
    range = [min(tensions),max(tensions)];
end

for i = 1:epithelium.edges.count
    if projection_edge(i) < 1
        continue
    end
    edge = epithelium.edges.objectAtIndex(i);
    if edge.isBoundary
        continue
    end

    T = p(projection_edge(i));

    mask.var(:) = false;
    v = edge.incidentVertices.objectAtIndex(1);
    w = edge.incidentVertices.objectAtIndex(2);
    ss_draw_line(mask,[v.y,v.x],[w.y,w.x],1);

    bitmap(mask.var) = T;
end

color = SYColor;
color.lut = 'violet-orange';
context = SYGraphicsContext(SYData(bitmap),frameSize(2),frameSize(1), ...
    8,SYGraphicsContext.CompositeModeOver, ...
    SYGraphicsContext.ColorSpaceIndexed,color.lut,range);
image = SYImage(context);

result = image;
end
