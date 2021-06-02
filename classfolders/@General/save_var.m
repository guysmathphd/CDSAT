function save_var(var,path,name)
varname = inputname(1);
if ~isfolder(path)
    mkdir(path);
end
if ~isempty(var)
    save(fullfile(path,name),'var','-v7.3');
end
end