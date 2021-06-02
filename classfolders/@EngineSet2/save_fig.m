function obj = save_fig(obj,f,name)
if ~isfolder(obj.figs_path)
    mkdir(obj.figs_path);
end
try
    saveas(f,fullfile(obj.path,'setfigs',name),'fig');
catch exception
    disp(exception.message);
end
try
    saveas(f,fullfile(obj.path,'setfigs',name),'png');
catch exception
    disp(exception.message);
end
end