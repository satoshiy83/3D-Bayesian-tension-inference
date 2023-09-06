% Short script to mark z dpeth of Ecdh signal.
% Written by Satoshi Yamashita.

function result = ss_mark_surface(image,mask,threshold)
% Function to mark z depth of fluorescence signal peak.
% result = ss_mark_surface(image,mask,threshold)
% Argument image is an image labeled for Ecdh or other junctional proteins.
% Argument mask is to exclude cell faces, and points of mask value equal 0
%   are measured for the depth.
% Argument threshold is a number for the minimum of signal intensity to
%   mark the depth.
% Return value is an SYImage instance of the depth map.

% Get bitmap data of image and mask.
bitmap = [];
for i = 1:image.representations.count
    bitmapRep = image.representations.objectAtIndex(i);
    if bitmapRep.repType == SYBitmapImageRep.RepTypeImage
        bitmap = cat(3,bitmap,bitmapRep.bitmap.var);
    end
end
bitmap = SYData(bitmap);

bitmapRep = mask.representations.objectAtIndex(1);
mask_ = bitmapRep.bitmap;

% Z depth.
h = size(bitmap.var,3);

% Draw z position of peak if it is brighter than the threshold.
siz = size(mask_.var);
depthMap = zeros(siz,'uint8');

indices = (find(mask_.var == 0))';
for i = indices
    [r,c] = ind2sub(siz,i);
    l = bitmap.var(r,c,:);
    [M,I] = max(l);
    
    if M > threshold
        depthMap(r,c) = h - I;
    end
end

result = SYImage(SYData(depthMap));
end
