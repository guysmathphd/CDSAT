function obj = set_logarithmic_bins(obj)
numbins = 15;
values_vec = General.load_var(fullfile(obj.path,'degree_vector.mat'));
tol = 1e-13;cond_vec = true(obj.N,1);

[ind_bins_var,bins_edges,ind_bins_var_sizes] = EngineClass.set_bins_generic(numbins,values_vec,tol,cond_vec);

General.save_var(ind_bins_var,obj.path,'ind_bins_var');
General.save_var(bins_edges,obj.path,'bins_edges');
General.save_var(ind_bins_var_sizes,obj.path,'ind_bins_var_sizes');



end