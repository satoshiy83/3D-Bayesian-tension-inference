% SA Class SATrack < SYObject.
% Written by Satoshi Yamashita.
% Data class for cell tracking.

classdef SATrack < SYObject
properties
    track_table = nan; % table.
    t_length = nan; % double.
    num_id = nan; % double.

    epithelia = nan; % SYDictionary.
    balanceDicts = nan; % SYDictionary.

    cell_track_match = nan; % double[n,m].
end

methods
function obj = initWithTable(obj,newTable)
% Initialization method with table.
% obj = initWithTable(obj,newTable)
% Argument newTable is a table made by readtable() from csv file.
    obj.track_table = newTable;
    obj.epithelia = SYDictionary;
    obj.balanceDicts = SYDictionary;
    
    obj.t_length = max(newTable.POSITION_T) + 1;
    obj.num_id = max(newTable.TRACK_ID) + 1;
    obj.cell_track_match = nan(obj.t_length,obj.num_id);
end
function obj = initWithData(obj,data)
% Initialization method with SYData instance.
% obj = initWithData(obj,data)
    obj.track_table = data.var.track_table;
    
    dict = SYDictionary;
    dict.initWithData(SYData(data.var.epithelia));
    obj.epithelia = dict;
    obj.t_length = data.var.t_length;
    obj.num_id = data.var.num_id;
    dict = SYDictionary;
    dict.initWithData(SYData(data.var.balanceDicts));
    obj.balanceDicts = dict;

    obj.cell_track_match = data.var.cell_track_match;
end
function dest = copy(obj,dest)
% Method to make a shallow copy.
% dest = copy(obj,dest)
% Return value is an instance of the class.
    if nargin < 2
        dest = SATrack;
    end
    copy@SYObject(obj,dest);

    dest.track_table = obj.track_table;
    dest.t_length = obj.t_length;
    dest.num_id = obj.num_id;
    dest.epithelia = obj.epithelia.copy;
    dest.balanceDicts = obj.balanceDicts.copy;
    dest.cell_track_match = obj.cell_track_match;
end

function result = data(obj)
% Method to convert the instance to an SYData instance.
% result = data(obj)
    s.track_table = obj.track_table;
    s.t_length = obj.t_length;
    s.num_id = obj.num_id;
    s.epithelia = obj.epithelia.data.var;
    s.balanceDicts = obj.balanceDicts.data.var;
    s.cell_track_match = obj.cell_track_match;

    result = SYData(s);
end

function mask_cell(~,mask,cel)
% Method to prepare a mask covering a cell.
% mask_cell(~,mask,cel)
% Argument mask is an SYData instance containing a raw data of the mask.
% Argument cel is the cell.
    frameSize = size(mask.var);
    neighborhood = [-frameSize(1),-1,1,frameSize(1)];
    m = length(mask.var(:));

    bitmap = SYData(zeros(frameSize));
    vertices = cel.encircleVertices;
    v0 = vertices.objectAtIndex(1);
    v_pre = v0;
    v = [v_pre.y,v_pre.x];
    v_list = round(v);
    for j = 2:vertices.count
        v_cur = vertices.objectAtIndex(j);

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
    index = sub2ind(frameSize,c(1),c(2));
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

function insertEpitheliumAt(obj,epithelium,T)
% Method to set epithelium data at a time point.
% insertEpitheliumAt(obj,epithelium,T)
% Argument epithelium is an SAEpithelium instance.
% Argument T is a number corresponding to the time point in table.
    obj.epithelia.setObjectForKey(string(T),epithelium);

    obj.matchCellAndTrack(T);
end
function matchCellAndTrack(obj,T)
% Method to match cells in an epithelium at a time point to the track.
% matchCellAndTrack(obj,T)
% Argument T is a number corresponding to the time point in table.
    epithelium = obj.epithelia.objectForKey(string(T));
    bitmap = obj.draw_cell_segment(epithelium);

    indices_time = obj.track_table.POSITION_T(:) == T - 1;
    for i = 1:obj.num_id
        indices_track = obj.track_table.TRACK_ID(:) == i - 1;
        indices = indices_time & indices_track;
        x = round(obj.track_table.POSITION_X(indices));
        y = round(obj.track_table.POSITION_Y(indices));

        if ~isnumeric(x) || length(x) ~= 1 || isnan(x) || ...
                ~isnumeric(y) || length(y) ~= 1 || isnan(y)
            disp('Invalid x and y.');
            disp('x: ');
            disp(x);
            disp('y: ');
            disp(y);
            continue
        end
        label = bitmap(y,x);
        if isnan(label)
            disp(['No cell was matched with track ',num2str(i), ...
                ' at x:',num2str(x),' y:',num2str(y)]);
            continue
        end

        obj.cell_track_match(T,i) = label;
    end
end
function result = draw_cell_segment(obj,epithelium)
% Method to prepare cell segmentation map.
% result = draw_cell_segment(obj,epithelium)
% Argument epithelium is an SAEpithelium instance.
% Return value is a 2D number array of map raw data.
    frameSize = epithelium.frameSize;
    bitmap = nan(frameSize);

    mask = SYData(zeros(frameSize,'logical'));

    for i = 1:epithelium.cells.count
        cel = epithelium.cells.objectAtIndex(i);
        if cel.vertices.count < 3
            continue
        end

        obj.mask_cell(mask,cel);

        bitmap(mask.var) = i;
    end

    result = bitmap;
end

function insertBalanceDictAt(obj,balanceDict,T)
% Method to set a dictionary holding infered tension data.
% insertBalanceDictAt(obj,balanceDict,T)
% Argument balanceDict is a dictionary holding infered tension data.
% Argument T is a number corresponding to the time point in table.
    obj.balanceDicts.setObjectForKey(string(T),balanceDict);
end

function result = drawExpansionRatio(obj)
% Method to draw a stack of heatmaps showing cell apical surface expansion
%   ratio.
% result = drawExpansionRatio(obj)
% Return value is a stack of heatmaps.
    Ts = find(any(~isnan(obj.cell_track_match),2))';

    if isempty(Ts) || length(Ts) < 2
        disp('Two epithelia must be matched to track.');
        result = nan;
        return
    end

    epithelium = obj.epithelia.objectForKey(string(Ts(1)));
    indices_cell = obj.cell_track_match(Ts(1),:);
    area_i = get_apical_area();

    epithelium = obj.epithelia.objectForKey(string(Ts(end)));
    indices_cell = obj.cell_track_match(Ts(end),:);
    area_t = get_apical_area();

    expansion_ratio = area_t ./ area_i;
    color = SYColor;
    color.lut = 'blue-red';
    m = max(abs(expansion_ratio - 1.0));
    range = [1 - m,1 + m];
    mask = SYData(zeros(epithelium.frameSize,'logical'));
    epithelium = obj.epithelia.objectForKey(string(Ts(1)));
    indices_cell = obj.cell_track_match(Ts(1),:);
    slice = draw_slice();
    image = SYImage(slice);
    for T = Ts(2:end)
        epithelium = obj.epithelia.objectForKey(string(T));
        indices_cell = obj.cell_track_match(T,:);
        slice = draw_slice();
        image.addRepresentation(slice);
    end

    function result = get_apical_area()
        array = nan(size(obj.cell_track_match,2),1);

        resolution = [epithelium.xResolution, ...
                      epithelium.yResolution, ...
                      epithelium.zResolution];
        for i = 1:length(indices_cell)
            index = indices_cell(i);
            if isnan(index)
                continue
            end

            cel = epithelium.cells.objectAtIndex(index);
            array(i) = cel.apicalArea(resolution);
        end

        result = array;
    end

    function result = draw_slice()
        frameSize = epithelium.frameSize;
        bitmap = nan(frameSize);

        for i = 1:length(indices_cell)
            index = indices_cell(i);
            if isnan(index)
                continue
            end

            cel = epithelium.cells.objectAtIndex(index);
            obj.mask_cell(mask,cel);

            bitmap(mask.var) = expansion_ratio(i);
        end

        context = SYGraphicsContext(SYData(bitmap), ...
            frameSize(2),frameSize(1), ...
            8,SYGraphicsContext.CompositeModeOver, ...
            SYGraphicsContext.ColorSpaceIndexed,color.lut,range);
        jmage = SYImage(context);
        bitmap = jmage.drawBitmapRep();

        mask.var(:) = false;
        
        for i = 1:epithelium.edges.count
            edge = epithelium.edges.objectAtIndex(i);
            v = edge.incidentVertices.objectAtIndex(1);
            w = edge.incidentVertices.objectAtIndex(2);
            ss_draw_line(mask,[v.y,v.x],[w.y,w.x],1);
        end

        bitmap(mask.var(:,:,[1,1,1])) = 0;

        result = SYBitmapImageRep(SYData(bitmap));
    end

    result = image;
end

function result = cellsToTracksAtTime(obj,cells,T)
% Method to convert cells labels to track.
% result = cellsToTracksAtTime(obj,cells,T)
% Argument cells is a number array of the cells labels.
% Argument T is a number corresponding to the time point in table.
% Return value is a number array of the track labels.
    labels = obj.cell_track_match(T,:);
    selected = any(labels == cells(:),1);

    result = find(selected);
end

function result = trackCellsArea(obj,cells,T)
% Method to get an array representing cells area time evolution.
% result = trackCellsArea(obj,cells,T)
% Argument cells is a number array of the cells labels.
% Argument T is a number corresponding to the time point in table.
% Return value is a 2D number array showing cell apical area, where its row
%   corresponds to time and its column corresponds to the cell.
    tracks = obj.cellsToTracksAtTime(cells,T);

    areas = nan(obj.t_length,length(tracks));

    for t = 1:obj.t_length
        if obj.epithelia.isNanForKey(string(t))
            continue
        end
        epithelium = obj.epithelia.objectForKey(string(t));
        resolution = [epithelium.xResolution, ...
                      epithelium.yResolution, ...
                      epithelium.zResolution];
        for i = 1:length(tracks)
            index = obj.cell_track_match(t,tracks(i));
            if isnan(index)
                continue
            end
            cel = epithelium.cells.objectAtIndex(index);
            areas(t,i) = cel.apicalArea(resolution);
        end
    end

    result = areas;
end
function result = averageCellArea(obj,excluded,T)
% Method to get a time evolution of average cell apical area.
% result = averageCellArea(obj,excluded,T)
% Argument excluded is a number array of cells labels to be excluded from
%   the calculation.
% Argument T is a number corresponding to the time point in table.
% Return value is a 2D array whose row corresponds to time and column
%   corresponds to average and standard deviation respectively.
    tracks = obj.cellsToTracksAtTime(excluded,T);

    areas = nan(obj.t_length,2);

    for t = 1:obj.t_length
        if obj.epithelia.isNanForKey(string(t))
            continue
        end
        epithelium = obj.epithelia.objectForKey(string(t));
        resolution = [epithelium.xResolution, ...
                      epithelium.yResolution, ...
                      epithelium.zResolution];
        array = zeros(epithelium.cells.count,1);
        count = 0;
        for i = 1:epithelium.cells.count
            if any(i == obj.cell_track_match(t,tracks))
                continue
            end
            cel = epithelium.cells.objectAtIndex(i);
            if cel.isBoundary
                continue
            end

            count = count + 1;
            array(count) = cel.apicalArea(resolution);
        end

        areas(t,1) = mean(array(1:count));
        areas(t,2) = std(array(1:count));
    end

    result = areas;
end
function result = trackCellsPressure(obj,cells,T)
% Method to get an array representing cells pressure time evolution.
% result = trackCellsPressure(obj,cells,T)
% Argument cells is a number array of the cells labels.
% Argument T is a number corresponding to the time point in table.
% Return value is a 2D number array showing cell pressure, where its row
%   corresponds to time and its column corresponds to the cell.
    tracks = obj.cellsToTracksAtTime(cells,T);

    pressures = nan(obj.t_length,length(tracks));

    for t = 1:obj.t_length
        if obj.balanceDicts.isNanForKey(string(t))
            continue
        end
        balanceDict = obj.balanceDicts.objectForKey(string(t));
        projection_cell = balanceDict.objectForKey("projection_cell");
        p = balanceDict.objectForKey("p");
        for i = 1:length(tracks)
            index = obj.cell_track_match(t,tracks(i));
            if isnan(index)
                continue
            end
            jndex = projection_cell(index);
            if jndex < 1
                continue
            end
            pressures(t,i) = p(jndex);
        end
    end

    result = pressures;
end
function result = averageCellPressure(obj,excluded,T)
% Method to get a time evolution of average cell pressure.
% result = averageCellPressure(obj,excluded,T)
% Argument excluded is a number array of cells labels to be excluded from
%   the calculation.
% Argument T is a number corresponding to the time point in table.
% Return value is a 2D array whose row corresponds to time and column
%   corresponds to average and standard deviation respectively.
    tracks = obj.cellsToTracksAtTime(excluded,T);

    pressures = nan(obj.t_length,2);

    for t = 1:obj.t_length
        if obj.epithelia.isNanForKey(string(t))
            continue
        end
        if obj.balanceDicts.isNanForKey(string(t))
            continue
        end
        balanceDict = obj.balanceDicts.objectForKey(string(t));
        projection_cell = balanceDict.objectForKey("projection_cell");
        p = balanceDict.objectForKey("p");

        array = zeros(length(projection_cell),1);
        count = 0;
        for i = 1:length(projection_cell)
            if any(i == obj.cell_track_match(t,tracks))
                continue
            end
            if projection_cell(i) == 0
                continue
            end

            count = count + 1;
            array(count) = p(projection_cell(i));
        end

        pressures(t,1) = mean(array(1:count));
        pressures(t,2) = std(array(1:count));
    end

    result = pressures;
end
function result = trackCellsSurfaceEnergy(obj,cells,T)
% Method to get an array representing cells surface energy time evolution.
% result = trackCellsSurfaceEnergy(obj,cells,T)
% Argument cells is a number array of the cells labels.
% Argument T is a number corresponding to the time point in table.
% Return value is a 2D number array showing cell junctional energy, where
%   its row corresponds to time and its column corresponds to the cell.
    tracks = obj.cellsToTracksAtTime(cells,T);

    surfaceEnergy = nan(obj.t_length,length(tracks));

    for t = 1:obj.t_length
        if obj.epithelia.isNanForKey(string(t))
            continue
        end
        balanceDict = obj.balanceDicts.objectForKey(string(t));
        epithelium = obj.epithelia.objectForKey(string(t));
        resolution = [epithelium.xResolution, ...
                      epithelium.yResolution, ...
                      epithelium.zResolution];
        for i = 1:length(tracks)
            index = obj.cell_track_match(t,tracks(i));
            if isnan(index)
                continue
            end
            cel = epithelium.cells.objectAtIndex(index);
            surfaceEnergy(t,i) = cel.surfaceEnergy(balanceDict,resolution);
        end
    end

    result = surfaceEnergy;
end
function result = averageCellSurfaceEnergy(obj,excluded,T)
% Method to get a time evolution of average cell surface energy.
% result = averageCellSurfaceEnergy(obj,excluded,T)
% Argument excluded is a number array of cells labels to be excluded from
%   the calculation.
% Argument T is a number corresponding to the time point in table.
% Return value is a 2D array whose row corresponds to time and column
%   corresponds to average and standard deviation respectively.
    tracks = obj.cellsToTracksAtTime(excluded,T);

    surfaceEnergies = nan(obj.t_length,2);

    for t = 1:obj.t_length
        if obj.epithelia.isNanForKey(string(t))
            continue
        end
        balanceDict = obj.balanceDicts.objectForKey(string(t));
        epithelium = obj.epithelia.objectForKey(string(t));
        resolution = [epithelium.xResolution, ...
                      epithelium.yResolution, ...
                      epithelium.zResolution];
        array = zeros(epithelium.cells.count,1);
        count = 0;
        for i = 1:epithelium.cells.count
            if any(i == obj.cell_track_match(t,tracks))
                continue
            end
            cel = epithelium.cells.objectAtIndex(i);
%             if cel.isBoundary
%                 continue
%             end
% 
%             count = count + 1;
%             array(count) = cel.surfaceEnergy(balanceDict,resolution);
            se = cel.surfaceEnergy(balanceDict,resolution);
            if isnan(se)
                continue
            end
            count = count + 1;
            array(count) = se;
        end

        surfaceEnergies(t,1) = mean(array(1:count));
        surfaceEnergies(t,2) = std(array(1:count));
    end

    result = surfaceEnergies;
end
function result = trackCellsMeanSurfaceTension(obj,cells,T)
% Method to get an array representing cells mean junctional tension time
%   evolution.
% result = trackCellsMeanSurfaceTension(obj,cells,T)
% Argument cells is a number array of the cells labels.
% Argument T is a number corresponding to the time point in table.
% Return value is a 2D number array showing cell mean junctional tension,
%   where its row corresponds to time and its column corresponds to the
%   cell.
    tracks = obj.cellsToTracksAtTime(cells,T);

    meanSurfaceTension = nan(obj.t_length,length(tracks));

    for t = 1:obj.t_length
        if obj.epithelia.isNanForKey(string(t))
            continue
        end
        balanceDict = obj.balanceDicts.objectForKey(string(t));
        epithelium = obj.epithelia.objectForKey(string(t));
        resolution = [epithelium.xResolution, ...
                      epithelium.yResolution, ...
                      epithelium.zResolution];
        for i = 1:length(tracks)
            index = obj.cell_track_match(t,tracks(i));
            if isnan(index)
                continue
            end
            cel = epithelium.cells.objectAtIndex(index);
            meanSurfaceTension(t,i) = ...
                cel.meanSurfaceTension(balanceDict,resolution);
        end
    end

    result = meanSurfaceTension;
end
function result = averageCellMeanSurfaceTension(obj,excluded,T)
% Method to get a time evolution of average cell mean junctional tension.
% result = averageCellMeanSurfaceTension(obj,excluded,T)
% Argument excluded is a number array of cells labels to be excluded from
%   the calculation.
% Argument T is a number corresponding to the time point in table.
% Return value is a 2D array whose row corresponds to time and column
%   corresponds to average and standard deviation respectively.
    tracks = obj.cellsToTracksAtTime(excluded,T);

    meanSurfaceTension = nan(obj.t_length,2);

    for t = 1:obj.t_length
        if obj.epithelia.isNanForKey(string(t))
            continue
        end
        balanceDict = obj.balanceDicts.objectForKey(string(t));
        epithelium = obj.epithelia.objectForKey(string(t));
        resolution = [epithelium.xResolution, ...
                      epithelium.yResolution, ...
                      epithelium.zResolution];
        array = zeros(epithelium.cells.count,1);
        count = 0;
        for i = 1:epithelium.cells.count
            if any(i == obj.cell_track_match(t,tracks))
                continue
            end
            cel = epithelium.cells.objectAtIndex(i);
%             if cel.isBoundary
%                 continue
%             end
% 
%             count = count + 1;
%             array(count) = cel.meanSurfaceTension(balanceDict,resolution);
            mst = cel.meanSurfaceTension(balanceDict,resolution);
            if isnan(mst)
                continue
            end
            count = count + 1;
            array(count) = mst;
        end

        meanSurfaceTension(t,1) = mean(array(1:count));
        meanSurfaceTension(t,2) = std(array(1:count));
    end

    result = meanSurfaceTension;
end

end
end
