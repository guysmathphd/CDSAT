function obj = create_BA(obj)
m = obj.BA_m;
num_nodes = 5000;
tmax = num_nodes - m;
A = zeros(5000);
a = ones(m) - eye(m);
A(1:m,1:m) = a;
for t = 1:tmax
    k = sum(A,2);
    cur_i = m + t;
    num_links = 0;
%     A(cur_i,1) = 0;
%     A(1,cur_i) = 0;
    while num_links < m
        for j = 1:cur_i - 1
            k_j = k(j);
            sum_k = sum(k);
            p_j = k_j/sum_k;
            if rand < p_j && ~A(cur_i,j)
                A(cur_i,j) = 1;
                A(j,cur_i) = 1;
                num_links = num_links + 1;
                if num_links == m
                    break
                end
            end
        end
    end
end
obj.N = size(A,2);
General.save_var(A,obj.path,'adjacency_matrix');
end