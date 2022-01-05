function obj = addengine(obj,foldernamestr,legendnamestr)
obj.folderNames{end+1} = foldernamestr;
obj.enginePaths{end+1} = fullfile('tests',foldernamestr,'results',[foldernamestr 'Obj']);
obj.objPropertiesPaths{end+1} = fullfile('tests',foldernamestr,'results','obj_properties');
obj.resultsPaths{end+1} = fullfile('tests',foldernamestr,'results');
obj.numEngines = obj.numEngines + 1;
obj.legendNames{end+1} = legendnamestr;
end