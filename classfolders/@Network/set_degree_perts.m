function obj = set_degree_perts(obj)
powers = [1,3,10];
path = fullfile(obj.path,'perts');
mydata = load(fullfile(obj.path,'degree_vector.mat'));k = mydata.var;
for i = powers
    ki = k.^i;
    kinorm = norm(ki);
    pert = ki/kinorm;
    pertname = ['pert_k_power_' num2str(i)];
    General.save_var(pert,path,pertname);
end
end
