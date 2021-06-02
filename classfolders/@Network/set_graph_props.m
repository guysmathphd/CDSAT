function obj = set_graph_props(obj)
A = General.load_var(fullfile(obj.path,'adjacency_matrix'));
G = graph(A);
d = distances(G);
k = sum(A,2);
k_inn = (A*k)./k;
General.save_var(G,obj.path,'graph_object');
General.save_var(d,obj.path,'shortest_distances');
General.save_var(k,obj.path,'degree_vector');
General.save_var(k_inn,obj.path,'nearest_neighbor_degree_vector');

end