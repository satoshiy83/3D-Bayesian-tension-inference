

function describe_inferred_tension(dict,fileName)

A = dict.objectForKey("balanceMatrix");
projection_edge = dict.objectForKey("projection_edge");
projection_cell = dict.objectForKey("projection_cell");
p = dict.objectForKey("p");

%% condition number
s = svd(A);
cn = s(1) / s(end);
txt = ['Condition number: ',num2str(cn),'\n'];

%% trivial solution residue
num_edge = sum(projection_edge > 0);
num_cell = sum(projection_cell > 0);
t = cat(1,zeros(num_edge,1),ones(num_cell,1));
residue_t = A * t;
txt = [txt,'Residue of trivial solution: ', ...
    num2str(mean(abs(residue_t))),' ± ',num2str(std(abs(residue_t))),'\n'];

%% inferred tension residue
residue_p = A * p;
txt = [txt,'Residue of inferred solution: ', ...
    num2str(mean(abs(residue_p))),' ± ',num2str(std(abs(residue_p))),'\n'];

%% inferred tension
txt = [txt,'Inferred tension and pressure: ', ...
    num2str(mean(p)),' ± ',num2str(std(p)),'\n'];

%% export
fid = fopen([fileName,'.txt'],'w');
fprintf(fid,txt);
fclose(fid);

end
