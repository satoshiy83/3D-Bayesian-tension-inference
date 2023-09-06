

function result = draw_cell_graph_2D(epithelium)

frameSize = epithelium.frameSize;
bitmap = zeros([frameSize,3],'uint8');
bitmap(:) = 255;
context = SYGraphicsContext(SYData(bitmap),frameSize(2),frameSize(1), ...
    8,SYGraphicsContext.CompositeModeOver, ...
    SYGraphicsContext.ColorSpaceRGB,nan,nan);
painter = SYPainter(context);
painter.lineWidth = 1.0;
painter.isAliased = true;

mask = SYData(zeros(frameSize,'logical'));

for i = 1:epithelium.edges.count
    edge = epithelium.edges.objectAtIndex(i);

    mask.var(:) = false;
    v = edge.incidentVertices.objectAtIndex(1);
    w = edge.incidentVertices.objectAtIndex(2);
    
    painter.removeAllPoints;
    startPoint.x = v.x; startPoint.y = v.y;
    endPoint.x = w.x; endPoint.y = w.y;
    painter.move(startPoint);
    painter.addLine(endPoint);

    painter.stroke([0,0,0]);
end

result = SYImage(context);
end
