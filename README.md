# 3D Bayesian tension inference
This project provides tools to implement the Bayesian tension inference with a 3D microscopic image labeled for adherens junction.

## Installation
Download files and put them in a folder with a suitable name. Go to Matlab command line, enter “addpath" + a full path of the folder, and enter “savepath”.

## Requirement
This project requires no Matlab toolbox, but a custom framework of objective classes SYObject family which is available at [![DOI](https://zenodo.org/badge/235579182.svg)](https://zenodo.org/badge/latestdoi/235579182).

## Example
Below is an example of the tension inference.
Assume that you have a 3D microscopic image labeled for adherens junction. The 3D image is first converted to x-y plane with maximum z projection, and cells are segmented by some means. The segmented image shows the cells with positive values and the cells were separated by 4-connected 1-pixel width boundary of 0 value. Let a name of the 3D image be smaple.tif and a name of the segmented image be segmented.tif.

### Converting image to vertices and edges representation
The image is first converted to vertices and edges representation, and a detected apical surface of the tissue is projected onto x-y plane. The projected surface will be clearer than the initial maximum z projection, and the cell segmentation can be done again with higher accuracy.
```
%% parameters.
image = "sample.tif";
segmented = "sample.seg.tif";
name = "sample";

threshold_brightness = 1000;

xResolution = 0.125;
yResolution = 0.125;
zResolution = 0.7;

threshold_merge = 2.0;
threshold_smooth = 3.0;

%% Initialization
image = SYImage(image);
segmented = SYImage(segmented);

depthMap = ss_mark_surface(image,segmented,threshold_brightness);

epithelium = SAEpithelium;
epithelium.depth = image.countStack;
epithelium.xResolution = xResolution;
epithelium.yResolution = yResolution;
epithelium.zResolution = zResolution;
epithelium.threshold_merge = threshold_merge;
epithelium.threshold_smooth = threshold_smooth;

%% Read data
epithelium.readSegmentedSlice(segmented);
epithelium.readDepthMap(depthMap);

epithelium.inferVerticesZ;
epithelium.mergeCloseVertices;
data = epithelium.data;
str = name + ".tb" + string(threshold_brightness) + ".tm" + ...
    string(threshold_merge) + ".ts" + string(threshold_smooth);
data.writeToFile(str + ".epit.mat");
gmage = epithelium.drawGraph;
gmage.writeToFile(str + ".epit.tif", false);

%% Draw apical surface
surface = ss_project_surface(segmented,epithelium);
amage = ss_draw_surface(image,surface,1);
amage.writeToFile(str + ".surf.tif", false);
```

The epithelium data can be later loaded from the .mat file.
```
data = SYData;
data.initWithContentsOfFile(str + “.epit.mat”);
epithelium = SAEpithelium;
epithelium.initWithData(data);

```

### Bayesian tension inference
The Bayesian tension inference is implemented as described in [![DOI](http://dx.doi.org/10.1016/j.jtbi.2012.08.017)].
The inferred tensions and pressures are expressed as a vector p and written into a dictionary returned by balanceMatrix() with a key “p”. Number arrays projection_edge and projection_cell map indices of edges and cells to the inferred result, where a tension on i-th edge is stored in p(projection_edge(i)) and a pressure in j-th cell is stored in p(projection_cell(j)). Both projection_edge and projection_cell are in the dictionary with keys of their names.
```
%% get balance matrix
dict = balanceMatrix(epithelium);

A = dict.objectForKey("balanceMatrix");
projection_edge = dict.objectForKey("projection_edge");
projection_cell = dict.objectForKey("projection_cell");

k_edge = sum(projection_edge > 0);
k_cell = sum(projection_cell > 0);

%% Bayesian inference
param_num = 2;
un = k_edge;
nkm = size(A,1) - k_cell + 1;
sSa = cat(2,A,zeros(size(A,1),1));
g = cat(1,ones(k_edge,1),zeros(k_cell,1));
B = diag(g);
sSb0 = cat(2,B,g);

mu = get_mu_for_min_ABIC(param_num,un,nkm,sSa,sSb0);

tau = sqrt(mu);
sSb = tau*sSb0;
S = cat(1,sSa,sSb);

[~,R] = qr(S);
H = R(1:end - 1,1:end - 1);
h = R(1:end - 1,end);

p = pinv(H)*h;
dict.setObjectForKey("mu",mu);
dict.setObjectForKey("p",p);

data = dict.data;
data.writeToFile(str + ".bm.p.mat");

%% drawing pressure
image = draw_tension_pressure(epithelium,dict);
image.writeToFile(str + ".tp.tif", true);
```

### Tracking cell apical area and pressure
Assume that you have a cell tracking data with Fiji plugin TrackMate and exported csv file, and epithelium and inferred tension data corresponding to 1st, 6th, and 11th time points in the tracking. Let files names be sample01.epit.mat, sample01.bm.p.mat, sample06.epit.mat, sample06.bm.p.mat, sample11.epit.mat, and sample11.bm.p.mat. Also you know labels of cells to be tracked at time = 1.
Then a tracking of area and pressure are exported as below.
```
%% parameters
tableName = “sample.csv”;
nameArray = SYArray(“sample01”, “sample06”, “sample11”);
TArray = [1, 6, 11];

Cells = [88,81,78,101,91,85,69,103,72,64,80];
t1 = 1;

%% initialization
tabl = readtable(tableName);
track = SATrack;
track.initWithTable(tabl);

for i = 1:nameArray.count
    str = nameArray.objectAtIndex(i);
    data = SYData;
    data.initWithContentsOfFile(str + “.epit.mat”);
    epit = SAEpithelium;
    epit.initWithData(data);

    data.initWithContentsOfFile(str + ”.bm.p.mat”);
    dict = SYDictionary;
    dict.initWithData(data);

    T = TArray(i);
    track.insertEpitheliumAt(epit,T);
    track.insertBalanceDictAt(dict,T);
end

%% export areas and pressures
areas = track.trackCellsArea(cells,t1);
pressures = track.trackCellsPressure(cells,t1);
areas_control = track.averageCellArea(cells,t1);
pressures_control = track.averageCellPressure(cells,t1);

dict = SYDictionary;
dict.setObjectForKey("track",track.data);
dict.setObjectForKey("cells",cells);
dict.setObjectForKey("t_cells",t);
dict.setObjectForKey("areas",areas);
dict.setObjectForKey("pressures",pressures);
dict.setObjectForKey("areas_control",areas_control);
dict.setObjectForKey("pressures_control",pressures_control);
data = dict.data;

data.writeToFile([name,'.area.pressure.mat']);
```
