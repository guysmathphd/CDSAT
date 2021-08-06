function obj = set_graph_props(obj)
A = General.load_var(fullfile(obj.path,'adjacency_matrix'));
k = sum(A,2);
N = size(A,1);
obj.N = N;
if ~isequal(A,A')
    G = digraph(A);
    k_in = sum(A,1);
    General.save_var(k_in,obj.path,'degree_vector_in');
    num_triangles_max = k.*k_in;
    L = sum(k) + sum(k_in);
    General.save_var(L,obj.path,'L_directed');
    k_avg = L/N;
    General.save_var(k_avg,obj.path,'k_avg_directed');
else
    G = graph(A);
    num_triangles_max = k.*(k-1)/2;
    L = .5 * sum(k);
    General.save_var(L,obj.path,'L_undirected');
    k_avg = sum(k)/N;
    General.save_var(k_avg,obj.path,'k_avg_undirected');
end

d = distances(G);
k_inn = (A*k)./k;
d_i_avg = mean(d,2);
d_avg = mean(d,'all');

if ~(any(A~=0 & A~=1,'all'))
    S = 2*sum(A,1)/(N-1);
    General.save_var(S,obj.path,'density');
    num_triangles = diag(A^3);
    C_i = num_triangles./num_triangles_max;
    General.save_var(C_i,obj.path,'clustering_coefficient');
    C = sum(C_i)/N;
    General.save_var(C,obj.path,'clustering_coefficient_avg');
end
if any(diag(A)~=0)
    A_noselfloops = A - diag(diag(A));
    General.save_var(A_noselfloops,obj.path,'adjacency_matrix_noselfloops');
end

sig_sq = var(k);
closeness_i = (sum(d,2)/N).^(-1);

General.save_var(G,obj.path,'graph_object');
General.save_var(d,obj.path,'shortest_distances');
General.save_var(k,obj.path,'degree_vector');
General.save_var(k_inn,obj.path,'nearest_neighbor_degree_vector');
General.save_var(d_i_avg,obj.path,'shortest_distance_i_avg');
General.save_var(d_avg,obj.path,'shortest_distances_avg');
General.save_var(sig_sq,obj.path,'degree_variance');
General.save_var(closeness_i,obj.path,'closeness_i');

end