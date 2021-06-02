function obj = batchFunction(obj,function_handle,inds)
for i = inds
    disp(['batch function i = ' num2str(i)]);
    eng_path = obj.enginePaths{i};
    eng = load(eng_path);
    function_handle(eng.obj);
%     clear eng;
end
end