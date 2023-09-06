

function result = histogram_residue(dict)

A = dict.objectForKey("balanceMatrix");
p = dict.objectForKey("p");
projection_edge = dict.objectForKey("projection_edge");

num_edge = sum(projection_edge > 0);

residue = A * p;
residue = sqrt(residue(1:2:end - 1) .^ 2 + residue(2:2:end) .^ 2);

B = A .* p';
B = sqrt(B(1:2:end - 1,:) .^ 2 + B(2:2:end,:) .^ 2);
B_T = B(:,1:num_edge);
B_T = B_T(B_T > 0);
B_P = B(:,num_edge + 1:end);
B_P = B_P(B_P > 0);

range = [0,max([residue',B_T',B_P'])];

hint = SYDictionary;
hint.setObjectForKey("hg_range",range);

hist_residue = histogram_scores(residue,hint);
hist_T = histogram_scores(B_T,hint);
hist_P = histogram_scores(B_P,hint);

result = cat(1,hist_residue(2:3,:),hist_T(3,:),hist_P(3,:));
end
