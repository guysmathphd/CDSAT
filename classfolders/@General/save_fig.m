function obj = save_fig(f,name,path)
if ~isfolder(path)
    mkdir(path);
end
try
    saveas(f,fullfile(path,name),'fig');
catch exception
    disp(exception.message);
end
try
    saveas(f,fullfile(path,name),'png');
catch exception
    disp(exception.message);
end
end