function obj = calc_eigs(obj)
A = General.load_var(fullfile(obj.path,'adjacency_matrix.mat'));
[v,d,flag] = eigs(A,obj.N,'largestreal');
if flag
    str = 'not';
else
    str = '';
end
disp(['Network.calc_eigs: eigenvalues did ' str ' converge']);
General.save_var(v,obj.path,'adj_mat_eigenvectors');
General.save_var(diag(d),obj.path,'adj_mat_eigenvalues');
end