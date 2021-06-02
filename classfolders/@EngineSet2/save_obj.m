function obj = save_obj(obj)
save(fullfile(obj.path,[obj.name 'Obj.mat']),'obj','-v7.3');
end