function obj = write_gephi_nodes_table_PEV(obj)
nodes = 1:obj.N;nodes=nodes';
v = General.load_var(fullfile(obj.path,'adj_mat_eigenvectors.mat'));
k = General.load_var(fullfile(obj.path,'degree_vector.mat'));
v1 = v(:,1);
values = abs(v1)/sum(abs(v1));
filepath = obj.path;
filename = 'gephi_nodes_PEV';
General.degree_radius_gephi_nodes_table(nodes,values,k,filepath,filename)
% General.general_gephi_nodes_table(nodes,values,x,y,filepath,filename);
end