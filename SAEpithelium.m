% SA Class SAEpithelium < SYObject.
% Written by Satoshi Yamashita.
% Data class representing an epithelium.

classdef SAEpithelium < SYObject
properties
    frameSize = [0,0]; % [int,int].
    depth = 0; % int.

    xResolution = nan; % float.
    yResolution = nan; % float.
    zResolution = nan; % float.

    cells = nan; % SYArray<SACell>.
    vertices = nan; % SYArray<SAVertex>.
    edges = nan; % SYArray<SAEdge>.

    % Adjustable parameters.
    threshold_smooth = 3.0;
    threshold_n_edge = 10;
    threshold_merge = 2.0;
end

methods
function obj = SAEpithelium
% Data class representing an epithelium.
% obj = SAEpithelium
    obj.cells = SYArray;
    obj.vertices = SYArray;
    obj.edges = SYArray;
end

function obj = initWithData(obj,data)
% Initialization method with SYData instance.
% obj = initWithData(obj,data)
    obj.frameSize = data.var.frameSize;
    obj.depth = data.var.depth;
    obj.xResolution = data.var.xResolution;
    obj.yResolution = data.var.yResolution;
    obj.zResolution = data.var.zResolution;
    obj.threshold_smooth = data.var.threshold_smooth;
    obj.threshold_n_edge = data.var.threshold_n_edge;
    obj.threshold_merge = data.var.threshold_merge;

    obj.cells.initWithData(SYData(data.var.cells));
    obj.vertices.initWithData(SYData(data.var.vertices));
    obj.edges.initWithData(SYData(data.var.edges));

    data_ = SYData(data.var.cells_vertices_indices);
    array = SYArray;
    array.initWithData(data_);
    data_ = SYData(data.var.cells_edges_indices);
    brray = SYArray;
    brray.initWithData(data_);
    for i = 1:obj.cells.count
        cel = obj.cells.objectAtIndex(i);
        indices = array.objectAtIndex(i);
        for j = indices
            cel.vertices.addObject(obj.vertices.objectAtIndex(j));
        end
        indices = brray.objectAtIndex(i);
        for j = indices
            cel.edges.addObject(obj.edges.objectAtIndex(j));
        end
    end

    data_ = SYData(data.var.vertices_cells_indices);
    array = SYArray;
    array.initWithData(data_);
    data_ = SYData(data.var.vertices_edges_indices);
    brray = SYArray;
    brray.initWithData(data_);
    for i = 1:obj.vertices.count
        vertex = obj.vertices.objectAtIndex(i);
        indices = array.objectAtIndex(i);
        for j = indices
            vertex.incidentCells.addObject(obj.cells.objectAtIndex(j));
        end
        indices = brray.objectAtIndex(i);
        for j = indices
            vertex.incidentEdges.addObject(obj.edges.objectAtIndex(j));
        end
    end

    data_ = SYData(data.var.edges_cells_indices);
    array = SYArray;
    array.initWithData(data_);
    data_ = SYData(data.var.edges_vertices_indices);
    brray = SYArray;
    brray.initWithData(data_);
    for i = 1:obj.edges.count
        edge = obj.edges.objectAtIndex(i);
        indices = array.objectAtIndex(i);
        for j = indices
            edge.incidentCells.addObject(obj.cells.objectAtIndex(j));
        end
        indices = brray.objectAtIndex(i);
        for j = indices
            edge.incidentVertices.addObject(obj.vertices.objectAtIndex(j));
        end
    end
end
function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = SAEpithelium;
    end
    copy@SYObject(obj,dest);

    dest.frameSize = obj.frameSize;
    dest.depth = obj.depth;
    dest.xResolution = obj.xResolution;
    dest.yResolution = obj.yResolution;
    dest.zResolution = obj.zResolution;
    dest.threshold_smooth = obj.threshold_smooth;
    dest.threshold_n_edge = obj.threshold_n_edge;
    dest.threshold_merge = obj.threshold_merge;

    dest.cells = obj.cells.copy;
    dest.vertices = obj.vertices.copy;
    dest.edges = obj.edges.copy;
end

function result = data(obj)
% Method to convert the instance to an SYData instance.
% result = data(obj)
    s.frameSize = obj.frameSize;
    s.depth = obj.depth;
    s.xResolution = obj.xResolution;
    s.yResolution = obj.yResolution;
    s.zResolution = obj.zResolution;
    s.threshold_smooth = obj.threshold_smooth;
    s.threshold_n_edge = obj.threshold_n_edge;
    s.threshold_merge = obj.threshold_merge;

    s.cells = obj.cells.data.var;
    s.vertices = obj.vertices.data.var;
    s.edges = obj.edges.data.var;
    
    array = SYArray;
    brray = SYArray;
    for i = 1:obj.cells.count
        cel = obj.cells.objectAtIndex(i);
        array.addObject(cel.indicesOfVertices);
        brray.addObject(cel.indicesOfEdges);
    end
    s.cells_vertices_indices = array.data.var;
    s.cells_edges_indices = brray.data.var;

    array = SYArray;
    brray = SYArray;
    for i = 1:obj.vertices.count
        vertex = obj.vertices.objectAtIndex(i);
        array.addObject(vertex.indicesOfCells);
        brray.addObject(vertex.indicesOfEdges);
    end
    s.vertices_cells_indices = array.data.var;
    s.vertices_edges_indices = brray.data.var;

    array = SYArray;
    brray = SYArray;
    for i = 1:obj.edges.count
        edge = obj.edges.objectAtIndex(i);
        array.addObject(edge.indicesOfCells);
        brray.addObject(edge.indicesOfVertices);
    end
    s.edges_cells_indices = array.data.var;
    s.edges_vertices_indices = brray.data.var;

    result = SYData(s);
end

function readSegmentedSlice(obj,segmented)
% Method to prepare an epithelium from a segmented image.
% readSegmentedSlice(obj,segmented)
% Argument segmented is an image or name of the image file, and the image 
%   shows cells by none zero value connected pixels separated by 0 value
%   pixels.
    if isstring(segmented)
        segmented = SYImage(segmented);
    elseif isnumeric(segmented)
        segmented = SYImage(SYData(segmented));
    end

    obj.cells.removeAllObjects;
    obj.vertices.removeAllObjects;
    obj.edges.removeAllObjects;
    
    % Prepare cells.
    labeled = IPConnectedComponents.connectedBinaryComponents(segmented,4);

    b = [labeled(1,:),labeled(end,:),labeled(:,1)',labeled(:,end)'];
    b = unique(b);

    for i = 1:max(labeled(:))
        cel = SACell;
        cel.label = i;

        if any(b == i)
            cel.isBoundary = true;
        else
            cel.isBoundary = false;
        end
        
        obj.cells.addObject(cel);
    end

    % Prepare vertices.
    siz = segmented.frameSize;
    obj.frameSize = siz;
    bitmapRep = segmented.representations.objectAtIndex(1);
    bitmap = bitmapRep.bitmap;
    indices = (find(bitmap.var == 0))';

%     vitmap = zeros(siz,'logical');

    label = 1;
    for i = indices
        [r,c] = ind2sub(siz,i);
        if r > 1 && r < siz(1) && c > 1 && c < siz(2)
            scope = labeled(r - 1:r + 1,c - 1:c + 1);
            inci = unique(scope(:));
            if length(inci) > 3
                vertex = SAVertex();
                vertex.label = label;
                vertex.x = double(c);
                vertex.y = double(r);
                vertex.isBoundary = false;
                inci(inci == 0) = [];
                for j = inci'
                    cel = obj.cells.objectAtIndex(j);
                    if cel.isBoundary
                        vertex.isBoundary = true;
                    end
                    vertex.incidentCells.addObject(cel);
                    cel.vertices.addObject(vertex);
                end

%                 vitmap(r,c) = true;

                obj.vertices.addObject(vertex);
                label = label + 1;
            elseif bitmap.var(r + 1,c) == 0 && ...
                    bitmap.var(r,c + 1) == 0 && ...
                    bitmap.var(r + 1,c + 1) == 0
                vertex = SAVertex();
                vertex.label = label;
                vertex.x = double(c) + 0.5;
                vertex.y = double(r) + 0.5;
                vertex.isBoundary = false;
                scope = labeled(r - 1:r + 2,c - 1:c + 2);
                inci = unique(scope(:));
                inci(inci == 0) = [];
                for j = inci'
                    cel = obj.cells.objectAtIndex(j);
                    if cel.isBoundary
                        vertex.isBoundary = true;
                    end
                    vertex.incidentCells.addObject(cel);
                    cel.vertices.addObject(vertex);
                end

                obj.vertices.addObject(vertex);
                label = label + 1;
            end
        end
    end

    % Prepare edges.
    label = 1;
    for i = 1:(obj.cells.count - 1)
        cel1 = obj.cells.objectAtIndex(i);

        citmap = labeled == cel1.label;
        citmap1 = citmap;
        citmap1(1:end - 1,:) = citmap1(1:end - 1,:) | citmap(2:end,:);
        citmap1(2:end,:) = citmap1(2:end,:) | citmap(1:end - 1,:);
        citmap1(:,1:end - 1) = citmap1(:,1:end - 1) | citmap(:,2:end);
        citmap1(:,2:end) = citmap1(:,2:end) | citmap(:,1:end - 1);

        for j = (i + 1):obj.cells.count
            cel2 = obj.cells.objectAtIndex(j);

            array = obj.verticesBetweenCells(cel1,cel2);
            if array.count < 2
                continue
            end

            citmap = labeled == cel2.label;
            citmap2 = citmap;
            citmap2(1:end - 1,:) = citmap2(1:end - 1,:) | citmap(2:end,:);
            citmap2(2:end,:) = citmap2(2:end,:) | citmap(1:end - 1,:);
            citmap2(:,1:end - 1) = citmap2(:,1:end - 1) | citmap(:,2:end);
            citmap2(:,2:end) = citmap2(:,2:end) | citmap(:,1:end - 1);

            eitmap = citmap1 & citmap2;

            for k = 1:array.count
                vertex = array.objectAtIndex(k);
                eitmap(round(vertex.y),round(vertex.x)) = true;
            end

            if array.count > 2
                for k = array.count:-1:1
                    vertex = array.objectAtIndex(k);
                    r = round(vertex.y);
                    c = round(vertex.x);
                    scope = eitmap(r - 1:r + 1,c - 1:c + 1);
                    scope(5) = false;
                    if ~any(scope(:))
                        array.removeObjectAtIndex(k);
                    end
                end
                if array.count ~= 2
                    disp('More than 2 incident vertices for an edge.');
                end
            end

            vertex1 = array.objectAtIndex(1);
            eitmap(round(vertex1.y),round(vertex1.x)) = true;
            vertex2 = array.objectAtIndex(2);
            eitmap(round(vertex2.y),round(vertex2.x)) = true;

            indices = find(eitmap);

            edge = SAEdge;
            edge.index = label;
            [r,c] = ind2sub(siz,indices);
            edge.x = c;
            edge.y = r;

            if cel1.isBoundary || cel2.isBoundary
                edge.isBoundary = true;
            else
                edge.isBoundary = false;
            end

            edge.incidentCells.addObjects(cel1,cel2);
            cel1.edges.addObject(edge);
            cel2.edges.addObject(edge);

            edge.incidentVertices.addObjects(vertex1,vertex2);
            vertex1.incidentEdges.addObject(edge);
            vertex2.incidentEdges.addObject(edge);

            obj.edges.addObject(edge);
            label = label + 1;
        end
    end
end
function result = drawSegmentedSlice(obj)
% Method to draw edges and vertices.
% result = drawSegmentedSlice(obj)
% Return value is a 2D image showing edges in red and vertices in cyan.
    bitmap = zeros([obj.frameSize,3],'uint8');

    for i = 1:obj.edges.count
        edge = obj.edges.objectAtIndex(i);
        indices = sub2ind(obj.frameSize,edge.y,edge.x);
        bitmap(indices) = 255;
    end

    for i = 1:obj.vertices.count
        vertex = obj.vertices.objectAtIndex(i);
        bitmap(round(vertex.y),round(vertex.x),:) = [0,255,255];
    end

    result = SYImage(SYData(bitmap));
end
function markCell(obj,index,image)
% Method to draw a cell into an image.
% markCell(obj,index,image)
% Argument index indicates the cell to be drawn.
% Argument image is an SYImage instance into which the cell is drawn.
% The cell is drawn with purple edges and white vertices.
    bitmapRep = image.representations.objectAtIndex(1);
    bitmap = bitmapRep.bitmap;

    mask = zeros(obj.frameSize,'logical');

    cel = obj.cells.objectAtIndex(index);
    for i = 1:cel.edges.count
        edge = cel.edges.objectAtIndex(i);
        indices = sub2ind(obj.frameSize,edge.y,edge.x);
        mask(indices) = true;
    end
    bitmap.var(mask(:,:,[1,1,1])) = repmat([255,0,255],sum(mask(:)),1);

    mask(:) = false;
    for i = 1:cel.vertices.count
        vertex = cel.vertices.objectAtIndex(i);
        mask(round(vertex.y),round(vertex.x)) = true;
    end
    bitmap.var(mask(:,:,[1,1,1])) = repmat([255,255,255],sum(mask(:)),1);
end

function result = vertexAt(obj,x,y)
% Method to get a vertex at given position.
% result = vertexAt(obj,x,y)
% Argument x and y are numbers specifying the position.
% Return value is an SAVertex instance or nan if not found.
    for i = 1:obj.vertices.count
        vertex = obj.vertices.objectAtIndex(i);
        if vertex.x == double(x) && vertex.y == double(y)
            result = vertex;
            return
        end
    end
    
    result = nan;
end
function result = verticesBetweenCells(~,cel1,cel2)
% Method to get vertices incident with two adjacent cells.
% result = verticesBetweenCells(~,cel1,cel2)
% Argument cel1 and cel2 are SACell instances adjacent to each other.
% Return value is an array of SAVertex instances incident with both of cel1
%   and cel2.
    array = SYArray;
    for i = 1:cel1.vertices.count
        vertex = cel1.vertices.objectAtIndex(i);
        if cel2.containsVertex(vertex)
            array.addObject(vertex);
        end
    end

    result = array;
end

function readDepthMap(obj,depthMap)
% Method to read a depth map.
% readDepthMap(obj,depthMap)
% Argument depthMap is the depth map image or its file name.
% Instance variable threshold_smooth evaluates a smoothness of the depth
%   map, where a difference between the depth of an edge pixel and average
%   depth of its adjacent edge pixels must be smaller than the threshold,
%   and otherwise the depth will be masked out.
    if isstring(depthMap)
        depthMap = SYImage(depthMap);
    elseif isnumeric(depthMap)
        depthMap = SYImage(SYData(depthMap));
    end

    bitmapRep = depthMap.representations.objectAtIndex(1);
    bitmap = bitmapRep.bitmap;

    for i = 1:obj.edges.count
        edge = obj.edges.objectAtIndex(i);

        indices = sub2ind(obj.frameSize,edge.y,edge.x);
        array = bitmap.var(indices);
        for j = 1:length(indices)
            if ~isSmooth(edge.x(j),edge.y(j))
                array(j) = 0;
            end
        end

        edge.z = array;
    end

    function result = isSmooth(x,y)
        scope = bitmap.var(y - 1:y + 1,x - 1:x + 1);
        scope(5) = 0;
        scope = scope(scope > 0);
        if length(scope) < 2
            result = false;
            return
        end
        m = mean(scope);
        result = abs(double(bitmap.var(y,x)) - m) <= obj.threshold_smooth;
    end
end
function inferVerticesZ(obj)
% Method to assign vertices with z depth.
% inferVerticesZ(obj)
% Instance variable threshold_n_edge defines the minimum number of pixels
%   of edges incident to a vertex and assigned depth, and a vertex with
%   incident edges without enough number of pixels with depth will be
%   assigned a depth averaged from adjacent vertices.
    array = SYArray;
    for i = 1:obj.vertices.count
        vertex = obj.vertices.objectAtIndex(i);
        
        x = [];
        y = [];
        z = [];
        for j = 1:vertex.incidentEdges.count
            edge = vertex.incidentEdges.objectAtIndex(j);
            indices = edge.z > 0;
            x = cat(1,x,edge.x(indices));
            y = cat(1,y,edge.y(indices));
            z = cat(1,z,edge.z(indices));
        end

        if length(z) < obj.threshold_n_edge
            array.addObject(vertex);
            continue
        end

        B = double([x,y,ones(length(x),1)]) \ double(z);
        vertex.z = [vertex.x,vertex.y,1] * B;
    end

    f = false;
    while true
        if array.count < 1 || f
            break
        end

        f = true;
        for i = array.count:-1:1
            vertex = array.objectAtIndex(i);
            z = [];
            for j = 1:vertex.incidentEdges.count
                edge = vertex.incidentEdges.objectAtIndex(j);
                wertex = edge.vertexOppositeFrom(vertex);
                if ~isnan(wertex.z)
                    z = cat(1,z,wertex.z);
                end
            end
            if ~isempty(z)
                vertex.z = mean(z);
                array.removeObjectAtIndex(i);
                f = false;
            end
        end
    end
end
function result = drawDepthMap(obj)
% Method to draw vertices and edges in a 3D stack.
% result = drawDepthMap(obj)
% Return value is an image in which the vertices are drawn with cyan and
%   the edges are drawn with green.
    stack = zeros([obj.frameSize,3,obj.depth],'uint8');

    for i = 1:obj.edges.count
        edge = obj.edges.objectAtIndex(i);
        for j = find(edge.z' > 0)
            I = obj.depth - edge.z(j);
            stack(edge.y(j),edge.x(j),2,I) = 255;
        end
    end

    for i = 1:obj.vertices.count
        vertex = obj.vertices.objectAtIndex(i);
        if ~isnan(vertex.z)
            r = round(vertex.y);
            c = round(vertex.x);
            I = obj.depth - round(vertex.z);
            if I < 2
                I = 1:2;
            elseif I == obj.depth
                I = (obj.depth - 1):obj.depth;
            else
                I = I - 1:I + 1;
            end
            stack(r - 1:r + 1,c - 1:c + 1,1,I) = 0;
            stack(r - 1:r + 1,c - 1:c + 1,2,I) = 255;
            stack(r - 1:r + 1,c - 1:c + 1,3,I) = 255;
        end
    end

    image = SYImage(SYData(stack(:,:,:,1)));
    for i = 2:obj.depth
        bitmapRep = SYBitmapImageRep(SYData(stack(:,:,:,i)));
        image.addRepresentation(bitmapRep);
    end

    result = image;
end

function mergeCloseVertices(obj)
% Method to merge vertices in close distance.
% mergeCloseVertices(obj)
% Instance variable threshold_merge defines the distance, where vertices
%   closer to each other than the threshold will be merged.
    x = zeros(obj.vertices.count,1);
    y = zeros(obj.vertices.count,1);
    z = zeros(obj.vertices.count,1);
    for i = 1:obj.vertices.count
        vertex = obj.vertices.objectAtIndex(i);
        x(i) = vertex.x;
        y(i) = vertex.y;
        z(i) = vertex.z;
    end

    dx = x - x';
    dy = y - y';
    dz = z - z';
    d = sqrt(dx .^ 2 + dy .^ 2 + dz .^ 2);
    c = d < obj.threshold_merge;
    c(logical(diag(ones(obj.vertices.count,1)))) = false;

    array = SYArray;
    for i = 1:obj.vertices.count
        cluster = zeros(obj.vertices.count,1,'logical');
        if any(c(i,:))
            cluster(i) = true;
            while true
                neighbor = any(c(:,cluster),2);
                neighbor(cluster) = false;
                if ~any(neighbor)
                    break
                end
                cluster(neighbor) = true;
            end
            c(cluster,cluster) = false;

            array.addObject(obj.vertices.objectsAtIndexes(find(cluster)));
        end
    end

    for i = 1:array.count
        cluster = array.objectAtIndex(i);

        % Make a vertex.
        vertex = SAVertex;

        x = 0.0; y = 0.0; z = 0.0;
        incidentCells = SYSet;
        incidentEdges = SYSet;
        for j = 1:cluster.count
            wertex = cluster.objectAtIndex(j);
            x = x + wertex.x;
            y = y + wertex.y;
            z = z + wertex.z;
            incidentCells.addObjectsFromArray(wertex.incidentCells);
            incidentEdges.addObjectsFromArray(wertex.incidentEdges);
        end
        incidentCells = incidentCells.allObjects;
        incidentEdges = incidentEdges.allObjects;

        includedEdges = SYArray;
        for j = 1:incidentEdges.count
            edge = incidentEdges.objectAtIndex(j);
            wertex = edge.incidentVertices.objectAtIndex(1);
            xertex = edge.incidentVertices.objectAtIndex(2);
            if cluster.containsObject(wertex) && ...
                    cluster.containsObject(xertex)
                includedEdges.addObject(edge);
            end
        end
        incidentEdges.removeObjectsInArray(includedEdges);

        vertex.x = x / cluster.count;
        vertex.y = y / cluster.count;
        vertex.z = z / cluster.count;
        vertex.isBoundary = false;
        vertex.incidentCells.addObjectsFromArray(incidentCells);
        vertex.incidentEdges.addObjectsFromArray(incidentEdges);

        % Update incident cells.
        for j = 1:incidentCells.count
            cel = incidentCells.objectAtIndex(j);
            cel.vertices.removeObjectsInArray(cluster);
            cel.edges.removeObjectsInArray(includedEdges);
            cel.vertices.addObject(vertex);

            if cel.isBoundary
                vertex.isBoundary = true;
            end
        end

        % Update incident edges.
        obj.edges.removeObjectsInArray(includedEdges);

        for j = 1:incidentEdges.count
            edge = incidentEdges.objectAtIndex(j);
            edge.incidentVertices.removeObjectsInArray(cluster);
            edge.incidentVertices.addObject(vertex);
        end

        % Update vertices array.
        obj.vertices.removeObjectsInArray(cluster);
        obj.vertices.addObject(vertex);
    end

    % Update index.
    for i = 1:obj.vertices.count
        vertex = obj.vertices.objectAtIndex(i);
        vertex.label = i;
    end
    for i = 1:obj.edges.count
        edge = obj.edges.objectAtIndex(i);
        edge.index = i;
    end
end
function result = drawGraph(obj)
% Method to draw vertices and edges in a 3D stack, where the edges are
%   drawn straight.
% result = drawGraph(obj)
% Return value is an image in which the vertices are drawn with cyan and
%   the edges are drawn with green.
    rStack = SYData(zeros([obj.frameSize,obj.depth],'uint8'));
    gStack = SYData(zeros([obj.frameSize,obj.depth],'uint8'));
    bStack = SYData(zeros([obj.frameSize,obj.depth],'uint8'));

    for i = 1:obj.edges.count
        edge = obj.edges.objectAtIndex(i);
        vertex = edge.incidentVertices.objectAtIndex(1);
        z = double(obj.depth) - vertex.z;
        if z < 1
            z = 1;
        elseif z > obj.depth
            z = double(obj.depth);
        end
        p = [vertex.y,vertex.x,z];
        vertex = edge.incidentVertices.objectAtIndex(2);
        z = double(obj.depth) - vertex.z;
        if z < 1
            z = 1;
        elseif z > obj.depth
            z = double(obj.depth);
        end
        q = [vertex.y,vertex.x,z];

        ss_draw_line(gStack,p,q,255);
    end

    for i = 1:obj.vertices.count
        vertex = obj.vertices.objectAtIndex(i);
        r = round(vertex.y);
        c = round(vertex.x);
        I = obj.depth - round(vertex.z);
        if I < 2
            I = 1:2;
        elseif I == obj.depth
            I = (obj.depth - 1):obj.depth;
        else
            I = I - 1:I + 1;
        end
        gStack.var(r - 1:r + 1,c - 1:c + 1,I) = 255;
        bStack.var(r - 1:r + 1,c - 1:c + 1,I) = 255;
    end

    bitmap = cat(3,rStack.var(:,:,1),gStack.var(:,:,1),bStack.var(:,:,1));
    image = SYImage(SYData(bitmap));
    for i = 2:obj.depth
        bitmap = cat(3,rStack.var(:,:,i),gStack.var(:,:,i),bStack.var(:,:,i));
        bitmapRep = SYBitmapImageRep(SYData(bitmap));
        image.addRepresentation(bitmapRep);
    end

    result = image;
end
function result = drawCellSegmentation(obj)
% Method to draw cells in grayscale, where their edges are drawn straight.
% result = drawCellSegmentation(obj)
    bitmap = zeros(obj.frameSize);

    mask = SYData(zeros(obj.frameSize,'logical'));

    for i = 1:obj.cells.count
        cel = obj.cells.objectAtIndex(i);
        if cel.isBoundary
            continue
        end
        
        obj.maskCell(mask,cel);

        bitmap(mask.var) = i;
    end

    result = SYImage(SYData(bitmap));
end
function maskCell(~,mask,cel)
% Method to make a mask covering a cell with straight edges.
% maskCell(~,mask,cel)
% Argument cel is an SACell instance to be masked.
    siz = size(mask.var);
    neighborhood = [-siz(1),-1,1,siz(1)];
    m = length(mask.var(:));

    mask.var(:) = false;

    bitmap = SYData(zeros(siz));
    circ = cel.encircleVertices;
    v0 = circ.objectAtIndex(1);
    v_pre = v0;
    v = [v_pre.y,v_pre.x];
    v_list = round(v);
    for j = 2:circ.count
        v_cur = circ.objectAtIndex(j);

        w = [v_cur.y,v_cur.x];
        ss_draw_line(bitmap,v,w,1);

        v = w;
        v_list = cat(1,v_list,round(v));
    end
    w = [v0.y,v0.x];
    ss_draw_line(bitmap,v,w,1);

    c = round(mean(v_list,1));
    x_min = min(v_list(:,2));
    x_max = max(v_list(:,2));
    y_min = min(v_list(:,1));
    y_max = max(v_list(:,1));
    index = sub2ind(siz,c(1),c(2));
    indicese = zeros((x_max - x_min)*(y_max - y_min),1);
    indicese(1) = index;
    a = 1;
    b = 2;
    bitmap.var(index) = 2;
    while a < b
        index = indicese(a);
        a = a + 1;

        neig = index + neighborhood;
        neig(neig <= 0 | neig > m) = [];
        neig(bitmap.var(neig) > 0) = [];
        if ~isempty(neig)
            bitmap.var(neig) = 2;
            indicese(b:b + length(neig) - 1) = neig;
            b = b + length(neig);

            % temporal treatment
            if b > (x_max - x_min) * (y_max - y_min)
                disp('flooding');
                break
            end
        end
    end
    mask.var(:) = bitmap.var(:) == 2;
end

function result = verticesAround(~,vertex)
% Method to get adjacent vertices for a vertex, where the vertices are
%   aligned in clockwise or couterclockwise.
% result = verticesAround(~,vertex)
% Argument vertex is an SAVertex instance for the center vertex.
    % Set initial vertex.
    edge = vertex.incidentEdges.objectAtIndex(1);
    v0 = edge.vertexOppositeFrom(vertex);
    arrayV = SYArray(v0);
    arrayE = SYArray(edge);

    % Find a cell containing the initial vertex.
    for i = 1:vertex.incidentCells.count
        cel = vertex.incidentCells.objectAtIndex(i);
        if cel.containsVertex(v0)
            break
        end
    end
    arrayC = SYArray(cel);

    % Enumerate incident cells and sort the adjacent vertices.
    v = v0;
    while true
        for i = 1:vertex.incidentEdges.count
            edge = vertex.incidentEdges.objectAtIndex(i);
            w = edge.vertexOppositeFrom(vertex);
            if v ~= w && cel.containsVertex(w)
                break
            end
        end
        if w == v0
            break
        end

        arrayV.addObject(w);
        arrayE.addObject(edge);

        for i = 1:vertex.incidentCells.count
            del = vertex.incidentCells.objectAtIndex(i);
            if del ~= cel && del.containsVertex(w)
                break
            end
        end
        arrayC.addObject(del);

        cel = del;
        v = w;
    end

    dict = SYDictionary;
    dict.setObjectForKey("vertices",arrayV);
    dict.setObjectForKey("edges",arrayE);
    dict.setObjectForKey("cells",arrayC);

    result = dict;
end

end
end
