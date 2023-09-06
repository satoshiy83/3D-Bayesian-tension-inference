

function result = track_condition_number(A,n)

vertices = ones(size(A,1)/2,1,'logical');
indices_cut = ones(size(A,1),1,'logical');
decrements = nan(length(vertices),n);
track = nan(1,n);
vertex_removed = zeros(1,n);
vertex_second = zeros(1,n);

for i = 1:n
    B = A(indices_cut,:);
    cn_cur = cond(B);
    for j = (find(vertices))'
        indices = indices_cut;
        indices([j * 2 - 1,j * 2]) = false;
        B = A(indices,:);

        cn = cond(B);
        decrements(j,i) = cn_cur - cn;
    end

    [d_cn,indices] = sort(decrements(vertices,i),'descend');
    if d_cn(1) > 0
        jndices = find(vertices);
        index = jndices(indices(1));
        vertex_removed(i) = index;
        vertices(index) = false;
        indices_cut([index * 2 - 1,index * 2]) = false;

        track(i) = cn_cur;

        vertex_second(i) = jndices(indices(2));
    else
        break
    end
end

txt = 'Track';
for i = 1:n
    if isnan(track(i))
        break
    end

    txt = cat(2,txt,'\t',num2str(track(i)));
end
txt = [txt,'\n'];

txt = [txt,'Removed'];
for i = 1:n
    if isnan(track(i))
        break
    end

    txt = cat(2,txt,'\t',num2str(vertex_removed(i)));
end
txt = [txt,'\n'];

txt = [txt,'Next'];
for i = 1:n
    if isnan(track(i))
        break
    end

    txt = cat(2,txt,'\t',num2str(vertex_second(i)));
end
txt = [txt,'\n'];

txt = [txt,'Decrements\n'];
for j = 1:size(decrements,1)
    txt = cat(2,txt,'Vertex_',num2str(j));
    for i = 1:n
        if isnan(track(i))
            break
        end

        txt = cat(2,txt,'\t');
        if ~isnan(decrements(j,i))
            txt = cat(2,txt,num2str(decrements(j,i)));
        end
    end
    txt = cat(2,txt,'\n');
end

dict = SYDictionary;
dict.setObjectForKey("Track",track);
dict.setObjectForKey("Vertex_removed",vertex_removed);
dict.setObjectForKey("Vertex_second",vertex_second);
dict.setObjectForKey("Decrements",decrements);
dict.setObjectForKey("Text",txt);

result = dict;
end
