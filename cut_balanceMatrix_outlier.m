

function result = cut_balanceMatrix_outlier(A,projection_edge, ...
    projection_cell,removed_vertices)

indices = ones(size(A,1),1,'logical');
indices(removed_vertices * 2 - 1) = false;
indices(removed_vertices * 2) = false;
B = A(indices,:);

remains = any(B,1);
B = B(:,remains);

indices = find(remains);

new_projection_edge = zeros(size(projection_edge));
new_projection_cell = zeros(size(projection_cell));
for i = 1:length(indices)
    index = projection_edge == indices(i);
    if any(index)
        new_projection_edge(index) = i;
        continue
    end

    index = projection_cell == indices(i);
    if any(index)
        new_projection_cell(index) = i;
    end
end

dict = SYDictionary;
dict.setObjectForKey("balanceMatrix",B);
dict.setObjectForKey("projection_edge",new_projection_edge);
dict.setObjectForKey("projection_cell",new_projection_cell);

result = dict;
end
