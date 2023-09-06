% Short script to project an epithelium surface.
% Written by Satoshi Yamashita.

function result = ss_draw_surface(image,depthMap,zWidth)
% Function to project an epithelium surface onto x-y plane.
% result = ss_draw_surface(image,depthMap,zWidth)
% Argument image is an image or name of the image file of the epithelium.
% Argument depthMap is a grayscale image where brightness indicates the
%   depth.
% Argument zWidth is a length to dilate the projecting space in z
%   direction.
% Return value is an image of the projected surface.

if isstring(image)
    image = SYImage(image);
elseif isnumeric(image)
    image = SYImage(SYData(image));
    image.splitChannels;
end

if isstring(depthMap)
    depthMap = SYImage(depthMap);
elseif isnumeric(depthMap)
    depthMap = SYImage(SYData(depthMap));
end

context = image.graphicsContext;

if context.colorSpace ~= SYGraphicsContext.ColorSpaceGrayscale
    disp('Image must be gray scale.');
end

bitmap = [];
for i = 1:image.representations.count
    bitmapRep = image.representations.objectAtIndex(i);
    if bitmapRep.repType == SYBitmapImageRep.RepTypeImage
        bitmap = cat(3,bitmap,bitmapRep.bitmap.var);
    end
end

bitsPerComponent = context.bitsPerComponent;
switch bitsPerComponent
    case 8
        d = 'uint8';
    case 16
        d = 'uint16';
    case 32
        d = 'single';
    case 64
        d = 'double';
    otherwise
        d = 'uint8';
end
citmap = zeros(image.frameSize,d);

h = size(bitmap,3);
bitmapRep = depthMap.representations.objectAtIndex(1);
ditmap = bitmapRep.bitmap.var;
indices = (find(ditmap > 0))';
for i = indices
    [r,c] = ind2sub(image.frameSize,i);
    
    z = h - round(ditmap(i));
    if z < 1
        z = 1;
    elseif z > h
        z = h;
    end
    z_i = max([1,z - zWidth]);
    z_t = min([h,z + zWidth]);
    if z_i > z_t
        disp('invalid z range.');
    end
    citmap(i) = max(bitmap(r,c,z_i:z_t));
end
indices = (find(ditmap == 0))';
for i = indices
    [r,c] = ind2sub(image.frameSize,i);
    citmap(i) = max(bitmap(r,c,:));
end

result = SYImage(SYData(citmap));
end
