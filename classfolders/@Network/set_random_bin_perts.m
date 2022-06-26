function obj = set_random_bin_perts(obj)
disp('Setting random bin perts...');
ind_bins_var = General.load_var(fullfile(obj.path,'ind_bins_var'));

pert = zeros(obj.N,1);
numbins = size(ind_bins_var,1);
for i1 = 1:numbins
    bin = ind_bins_var{i1};
    node = bin(1);
    rand_pert = rand + .5;
    pert(node) = rand_pert;
end
General.save_var(pert,fullfile(obj.path,'random_bin_perts'),'pert1');
disp('Done');
end