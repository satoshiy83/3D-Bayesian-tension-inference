

function result = draw_tension_pressure(epithelium,balanceDict)

image_T = draw_tension(epithelium,balanceDict);
image_P = draw_pressure(epithelium,balanceDict);

bitmap = image_P.drawBitmapRep();

citmap = image_T.drawBitmapRep();
bitmapRep = image_T.representations.objectAtIndex(1);
mask = ~isnan(bitmapRep.bitmap.var);

mask = repmat(mask,1,1,3);
bitmap(mask) = citmap(mask);

result = SYImage(SYData(bitmap));
end
