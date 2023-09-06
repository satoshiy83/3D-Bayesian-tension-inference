

function result = draw_surfaceEnergy(epithelium,balanceDict)

frameSize = epithelium.frameSize;
bitmap = nan(frameSize);
resolution = [epithelium.xResolution, ...
              epithelium.yResolution, ...
              epithelium.zResolution];

mask = SYData(zeros(frameSize,'uint8'));
nask = SYData(zeros(frameSize,'logical'));
neighborhood = [-frameSize(1),-1,1,frameSize(1)];
m = length(mask.var(:));

surfaceEnergies = nan([epithelium.cells.count,1]);
for i = 1:epithelium.cells.count
    cel = epithelium.cells.objectAtIndex(i);
    surfaceEnergies(i) = cel.surfaceEnergy(balanceDict,resolution);
end

if ~balanceDict.isNanForKey("range_meanTension")
    range = balanceDict.objectForKey("range_meanTension");
else
    range = [min(surfaceEnergies),max(surfaceEnergies)];
end

for i = 1:epithelium.cells.count
    if isnan(surfaceEnergies(i))
        continue
    end

    cel = epithelium.cells.objectAtIndex(i);
    
    mask.var(:) = 0;

    vertices = cel.encircleVertices;
    v0 = vertices.objectAtIndex(1);
    v_pre = v0;
    v = [v_pre.y,v_pre.x];
    v_list = round(v);
    for j = 2:vertices.count
        v_cur = vertices.objectAtIndex(j);

        w = [v_cur.y,v_cur.x];
        ss_draw_line(mask,v,w,1);

        v = w;
        v_list = cat(1,v_list,round(v));
    end
    w = [v0.y,v0.x];
    ss_draw_line(mask,v,w,1);

    c = round(mean(v_list,1));
    x_min = min(v_list(:,2));
    x_max = max(v_list(:,2));
    y_min = min(v_list(:,1));
    y_max = max(v_list(:,1));
    index = sub2ind(frameSize,c(1),c(2));
    indicese = zeros((x_max - x_min)*(y_max - y_min),1);
    indicese(1) = index;
    a = 1;
    b = 2;
    mask.var(index) = 2;
    while a < b
        index = indicese(a);
        a = a + 1;

        neig = index + neighborhood;
        neig(neig <= 0 | neig > m) = [];
        neig(mask.var(neig) > 0) = [];
        if ~isempty(neig)
            mask.var(neig) = 2;
            indicese(b:b + length(neig) - 1) = neig;
            b = b + length(neig);

            % temporal treatment
            if b > (x_max - x_min) * (y_max - y_min)
                disp('flooding');
                break
            end
        end
    end
    nask.var(:) = mask.var(:) == 2;

    bitmap(nask.var) = surfaceEnergies(i);
end

color = SYColor;
color.lut = 'violet-orange';
context = SYGraphicsContext(SYData(bitmap),frameSize(2),frameSize(1), ...
    8,SYGraphicsContext.CompositeModeOver, ...
    SYGraphicsContext.ColorSpaceIndexed,color.lut,range);
image = SYImage(context);

result = image;
end
