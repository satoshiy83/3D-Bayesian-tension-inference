% SA Class SACell < SYObject.
% Written by Satoshi Yamashita.
% Data class representing a cell.

classdef SACell < SYObject
properties
    label = nan; % int.

    isBoundary = nan; % bool.

    vertices = nan; % SYArray<SAVertex>.
    edges = nan; % SYArray<SAEdge>.
end

methods
function obj = SACell
% Data class representing a cell.
% obj = SACell
    obj.label = nan;
    obj.isBoundary = nan;

    obj.vertices = SYArray;
    obj.edges = SYArray;
end

function obj = initWithData(obj,data)
% Initialization method with SYData instance.
% obj = initWithData(obj,data)
    s = data.var;

    obj.label = s.label;
    obj.isBoundary = s.isBoundary;
end
function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = SACell;
    end
    copy@SYObject(obj,dest);

    dest.label = obj.label;
    dest.isBoundary = obj.isBoundary;
    dest.vertices = obj.vertices.copy;
    dest.edges = obj.edges.copy;
end

function result = data(obj)
% Method to convert the instance to an SYData instance.
% result = data(obj)
    s.label = obj.label;
    s.isBoundary = obj.isBoundary;

    result = SYData(s);
end

function result = indicesOfVertices(obj)
% Method to get an array of incident vertices labels.
% result = indicesOfVertices(obj)
% Return value is a number array containing labels of the incident
%   vertices.
    array = zeros(1,obj.vertices.count);
    for i = 1:obj.vertices.count
        vertex = obj.vertices.objectAtIndex(i);
        array(i) = vertex.label;
    end
    result = array;
end
function result = indicesOfEdges(obj)
% Method to get an array of incident edges labels.
% result = indicesOfEdges(obj)
% Return value is a number array containing labels of the incident edges.
    array = zeros(1,obj.edges.count);
    for i = 1:obj.edges.count
        edge = obj.edges.objectAtIndex(i);
        array(i) = edge.index;
    end
    result = array;
end

function result = containsVertex(obj,vertex)
% Method to ask if the cell is incident with a vertex.
% result = containsVertex(obj,vertex)
% Argument vertex is the SAVertex instance.
% Return value is a boolean.
    result = obj.vertices.containsObject(vertex);
end
function result = isAdjacentWithCell(obj,cel)
% Method to ask if the cell is adjacent with a cell.
% result = isAdjacentWithCell(obj,cel)
% Argument cell is the SACell instance.
% Return value is a boolean.
    b = false;
    for i = 1:obj.edges.count
        edge = obj.edges.objectAtIndex(i);
        del = edge.cellOppositeFrom(obj);
        if cel == del
            b = true;
            break
        end
    end
    result = b;
end
function result = isAdjacentWithBoundaryCell(obj)
% Method to ask if the cell is adjacent with a boundary cell.
% result = isAdjacentWithBoundaryCell(obj)
% Return value is a boolean.
    b = false;
    for i = 1:obj.edges.count
        edge = obj.edges.objectAtIndex(i);
        cel = edge.cellOppositeFrom(obj);
        if cel.isBoundary
            b = true;
            break
        end
    end
    result = b;
end

function result = encircleVertices(obj)
% Method to get an aligned vertices incident with the cell.
% result = encircleVertices(obj)
% Return value is an array of incident vertice aligned in clockwise or
%   counterclockwise.
    if obj.vertices.count < 4
        result = obj.vertices.copy;
        return
    end

    v0 = obj.vertices.objectAtIndex(1);
    for i = 1:v0.incidentEdges.count
        edge = v0.incidentEdges.objectAtIndex(i);
        v1 = edge.vertexOppositeFrom(v0);
        if obj.containsVertex(v1)
            break
        end
    end

    v_pre = v0;
    v_cur = v1;
    array = SYArray(v0,v1);
    for c = 3:obj.vertices.count
        for i = 1:v_cur.incidentEdges.count
            edge = v_cur.incidentEdges.objectAtIndex(i);
            v = edge.vertexOppositeFrom(v_cur);
            if v ~= v_pre && obj.containsVertex(v)
                break
            end
        end
        array.addObject(v);
        v_pre = v_cur;
        v_cur = v;
    end

    result = array;
end

function result = apicalArea(obj,resolution)
% Method to get an apical area of the cell.
% result = apicalArea(obj,resolution)
% Argument resolution is a number array indicating x, y, and z resolution.
% Return value is a number representing the area.
    array = obj.encircleVertices;

    xyz = [];
    for i = 1:array.count
        vertex = array.objectAtIndex(i);
        xyz = cat(1,xyz,[vertex.x,vertex.y,vertex.z] .* resolution);
    end
    [n,c] = normalOfFittingPlane(xyz);

    v = xyz - c;
    d = -dot(repmat(n,[array.count,1]),v,2);
    xyz_projected = xyz + d * n - c;

    v_a = cross(xyz_projected,xyz_projected([2:array.count,1],:),2);
    v_a = sum(v_a,1);
    
    result = sqrt(sum(v_a .^ 2)) / 2;
end
function result = perimeterLength(obj,resolution)
% Method to get a perimeter length of the cell.
% result = perimeterLength(obj,resolution)
% Argument resolution is a number array indicating x, y, and z resolution.
% Return value is a number representing the length.
    l = 0;
    for i = 1:obj.edges.count
        edge = obj.edges.objectAtIndex(i);
        v = edge.incidentVertices.objectAtIndex(1);
        w = edge.incidentVertices.objectAtIndex(2);

        xyz = ([v.x,v.y,v.z] - [w.x,w.y,w.z]) .* resolution;
        d = sqrt(sum(xyz .^ 2));

        l = l + d;
    end
    
    result = l;
end
function result = surfaceEnergy(obj,balanceDict,resolution)
% Method to get a sum of junctional tension weighted by junction length.
% result = surfaceEnergy(obj,balanceDict,resolution)
% Argument balanceDict is a dictionary made by balanceMatrix().
% Argument resolution is a number array indicating x, y, and z resolution.
% Return value is a number representing the surface energy.
    if obj.isBoundary || obj.isAdjacentWithBoundaryCell
        result = nan;
        return
    end

    e = 0;

    projection_edge = balanceDict.objectForKey("projection_edge");
    p = balanceDict.objectForKey("p");
    for i = 1:obj.edges.count
        edge = obj.edges.objectAtIndex(i);
        v = edge.incidentVertices.objectAtIndex(1);
        w = edge.incidentVertices.objectAtIndex(2);
        T = p(projection_edge(edge.index));

        xyz = ([v.x,v.y,v.z] - [w.x,w.y,w.z]) .* resolution;
        d = sqrt(sum(xyz .^ 2));

        e = e + d * T;
    end

    result = e;
end
function result = meanSurfaceTension(obj,balanceDict,resolution)
% Method to get a mean junctional tension weighted by junction length.
% result = meanSurfaceTension(obj,balanceDict,resolution)
% Argument balanceDict is a dictionary made by balanceMatrix().
% Argument resolution is a number array indicating x, y, and z resolution.
% Return value is a number representing the mean tension.
    if obj.isBoundary || obj.isAdjacentWithBoundaryCell
        result = nan;
        return
    end

    e = obj.surfaceEnergy(balanceDict,resolution);
    l = obj.perimeterLength(resolution);

    result = e / l;
end

end
end
