function obj = write_gephi_edges_file(obj)
A = General.load_var(fullfile(obj.path,'adjacency_matrix.mat'));
Source = [];Target = []; Type = {}; Weight = [];
if isequal(A,A')
    typeStr = 'Undirected';
else
    typeStr = 'Directed';
end
for i =1:size(A,1)
    for j = 1:size(A,2)
        if A(i,j)
            Source(end+1,1) = i;
            Target(end+1,1) = j;
            Type{end+1,1} = typeStr;
            Weight(end+1,1) = A(i,j);
        end
    end
end
T = table(Source,Target,Type,Weight);
writetable(T,fullfile(obj.path,'gephi_edges.csv'));
end