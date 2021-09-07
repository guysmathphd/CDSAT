function batchFunction(set_obj,function_handle,inds)
for i = inds
    disp(['batch function i = ' num2str(i)]);
    cur_obj_path = set_obj.Paths{i};
    cur_obj = load(cur_obj_path);
    function_handle(cur_obj.obj);
%     clear eng;
end
end