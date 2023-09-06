% SA Class SAEdge < SYObject.
% Written by  Satoshi Yamashita.
% Data class representing an edge.

classdef SAEdge < SYObject
properties
    index = nan; % int.

    x = []; % [int].
    y = []; % [int].
    z = []; % [int].

    isBoundary = nan; % bool.

    incidentCells = nan; % SYArray<SACell>.
    incidentVertices = nan; % SYArray<SAVertex>.
end

methods
function obj = SAEdge
% Data class representing an edge.
% obj = SAEdge
    obj.x = [];
    obj.y = [];
    obj.z = [];

    obj.isBoundary = nan;

    obj.incidentCells = SYArray;
    obj.incidentVertices = SYArray;
end

function obj = initWithData(obj,data)
% Initialization method with SYData instance.
% obj = initWithData(obj,data)
    s = data.var;

    obj.index = s.index;
    obj.x = s.x;
    obj.y = s.y;
    obj.z = s.z;
    obj.isBoundary = s.isBoundary;
end
function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = SAEdge;
    end
    copy@SYObject(obj,dest);

    dest.index = obj.index;
    dest.x = obj.x;
    dest.y = obj.y;
    dest.z = obj.z;
    dest.isBoundary = obj.isBoundary;

    dest.incidentCells = obj.incidentCells.copy;
    dest.incidentVertices = obj.incidentVertices.copy;
end

function result = data(obj)
% Method to convert the instance to an SYData instance.
% result = data(obj)
    s.index = obj.index;
    s.x = obj.x;
    s.y = obj.y;
    s.z = obj.z;
    s.isBoundary = obj.isBoundary;

    result = SYData(s);
end

function result = indicesOfCells(obj)
% Method to get an array of incident cells labels.
% result = indicesOfCells(obj)
% Return value is a number array containing labels of the incident cells.
    array = zeros(1,obj.incidentCells.count);
    for i = 1:obj.incidentCells.count
        cel = obj.incidentCells.objectAtIndex(i);
        array(i) = cel.label;
    end
    result = array;
end
function result = indicesOfVertices(obj)
% Method to get an array of labels of vertices linked by the edge.
% result = indicesOfVertices(obj)
% Return value is a number array containing labels of the incident
%   vertices.
    array = zeros(1,obj.incidentVertices.count);
    for i = 1:obj.incidentVertices.count
        vertex = obj.incidentVertices.objectAtIndex(i);
        array(i) = vertex.label;
    end
    result = array;
end

function result = vertexOppositeFrom(obj,vertex)
% Method to get an incident vertex other than a given vertex.
% vertexOppositeFrom(obj,vertex)
% Argument vertex is an SAVertex instance.
% Return value is an SAVertex instance opposite from the argeument or nan
%   if the argument vertex was not incident with the edge.
    if vertex == obj.incidentVertices.objectAtIndex(1)
        result = obj.incidentVertices.objectAtIndex(2);
    elseif vertex == obj.incidentVertices.objectAtIndex(2)
        result = obj.incidentVertices.objectAtIndex(1);
    else
        result = nan;
    end
end
function result = cellOppositeFrom(obj,cel)
% Method to get an incident cell other than a gicen cell.
% result = cellOppositeFrom(obj,cel)
% Argument cel is an SACell instance opposite from the argument or nan if
%   the argument cell was not incident with the edge.
% Return value is an SACell instance.
    if cel == obj.incidentCells.objectAtIndex(1)
        result = obj.incidentCells.objectAtIndex(2);
    elseif cel == obj.incidentCells.objectAtIndex(2)
        result = obj.incidentCells.objectAtIndex(1);
    else
        result = nan;
    end
end

function result = edgeLength(obj,resolution)
% Method to get a length of the edge.
% result = edgeLength(obj,resolution)
% Argument resolution is a number array indicating x, y, and z resolution.
% Return value is a number.
    vertex = obj.incidentVertices.objectAtIndex(1);
    wertex = obj.incidentVertices.objectAtIndex(2);
    v = [vertex.x,vertex.y,vertex.z] - [wertex.x,wertex.y,wertex.z];
    v = v .* resolution;

    result = sqrt(sum(v .^ 2));
end

end
end
