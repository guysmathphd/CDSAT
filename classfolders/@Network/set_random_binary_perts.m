function obj = set_random_binary_perts(obj)
num_perts = 1e5;
% perts = zeros(obj.N,num_perts);
perts = logical(randi([0, 1], [obj.N, num_perts]));
General.save_var(perts,obj.path,'rand_binary_perts');



end