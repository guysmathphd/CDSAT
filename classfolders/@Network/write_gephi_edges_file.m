function obj = write_gephi_edges_file(obj,DT,filenamesufstr)

Source = [];Target = []; Type = {}; Weight = [];
if isempty(DT)
    A = General.load_var(fullfile(obj.path,'adjacency_matrix.mat'));
else
    pu_datetime = General.load_var(fullfile(obj.source_data_path,[obj.source_data_name '_filtered_pu_datetime']));
    do_datetime = General.load_var(fullfile(obj.source_data_path,[obj.source_data_name '_filtered_do_datetime']));
    pu_locID = General.load_var(fullfile(obj.source_data_path,[obj.source_data_name '_filtered_pu_locID']));
    do_locID = General.load_var(fullfile(obj.source_data_path,[obj.source_data_name '_filtered_do_locID']));
    A = zeros(max([pu_locID;do_locID]));
    n = size(pu_datetime,1);
    for i1 = 1:n
        if pu_datetime(i1) < DT && do_datetime(i1) > DT
            A(pu_locID(i1),do_locID(i1)) = A(pu_locID(i1),do_locID(i1)) + 1;
        end
    end
end
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
writetable(T,fullfile(obj.path,['gephi_edges_' filenamesufstr '.csv']));
end