% SA Class SAVertex < SYObject.
% Written by Satoshi Yamashita.
% Data class representing a vertex.

classdef SAVertex < SYObject
properties
    label = nan; % int.

    x = nan; % float.
    y = nan; % float.
    z = nan; % float.

    isBoundary = nan; % bool.

    incidentCells = nan; % SYArray<SACell>.
    incidentEdges = nan; % SYArray<SAEdge>.
end

methods
function obj = SAVertex
% Data class representing a vertex
% obj = SAVertex
    obj.label = nan;

    obj.x = nan;
    obj.y = nan;
    obj.z = nan;

    obj.isBoundary = nan;

    obj.incidentCells = SYArray;
    obj.incidentEdges = SYArray;
end

function obj = initWithData(obj,data)
% Initialization method with SYData instance.
% obj = initWithData(obj,data)
    s = data.var;

    obj.label = s.label;
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
        dest = SAVertex;
    end
    copy@SYObject(obj,dest);

    dest.label = obj.label;
    dest.x = obj.x;
    dest.y = obj.y;
    dest.z = obj.z;
    dest.isBoundary = obj.isBoundary;

    dest.incidentCells = obj.incidentCells.copy;
    dest.incidentEdges = obj.incidentEdges.copy;
end

function result = data(obj)
% Method to convert the instance to an SYData instance.
% result = data(obj)
    s.label = obj.label;
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
function result = indicesOfEdges(obj)
% Method to get an array of incident edges labels.
% result = indicesOfEdges(obj)
% Return value is a number array containing labels of the incident edges.
    array = zeros(1,obj.incidentEdges.count);
    for i = 1:obj.incidentEdges.count
        edge = obj.incidentEdges.objectAtIndex(i);
        array(i) = edge.index;
    end
    result = array;
end

end
end
