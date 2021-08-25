function var = load_var(filepath)
mydata = load(filepath);
if isfield(mydata,'var')
    var = mydata.var;
elseif isfield(mydata,'obj')
    var = mydata.obj;
end
end