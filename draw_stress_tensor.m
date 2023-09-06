

function result = draw_stress_tensor(context, epithelium, balanceDict, ...
    clusters, scale_factor, lineWidth, color)

painter = SYPainter(context);
painter.lineWidth = lineWidth;
painter.isAliased = false;

for i = 1:clusters.count
    clus = clusters.objectAtIndex(i);

    % get center.
    c_list = [];
    for j = clus
        cel = epithelium.cells.objectAtIndex(j);
        v_list = [];
        for k = 1:cel.vertices.count
            v = cel.vertices.objectAtIndex(k);
            v_list = cat(1,v_list,[v.y, v.x]);
        end
        c_list = cat(1,c_list,mean(v_list,1));
    end
    center = mean(c_list,1);

    % plot tensor.
    T = stressTensor(epithelium,balanceDict,clus);
    [V,D] = eig(T);
    [m,indices] = sort(diag(D),'descend');
    v = V(:,indices(1));

    painter.removeAllPoints;
    startPoint.x = center(2) - v(1) * m(1) * scale_factor;
    startPoint.y = center(1) - v(2) * m(1) * scale_factor;
    painter.move(startPoint);
    endPoint.x = center(2) + v(1) * m(1) * scale_factor;
    endPoint.y = center(1) + v(2) * m(1) * scale_factor;
    painter.addLine(endPoint);

    painter.stroke(color);
end

result = SYImage(context.copy);
end
