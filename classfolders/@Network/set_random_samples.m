function obj = set_random_samples(obj)
path = fullfile(obj.path,'random_samples');
num_nodes = [1 10 100 1000];
for n = num_nodes
    nodes = datasample(1:obj.N,n,'Replace',false);
    sample_name = ['n_' num2str(n)];
    General.save_var(nodes,path,sample_name);
end
end