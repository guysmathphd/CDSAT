function obj = create_ER(obj)
A = rand(obj.N);
A = A < obj.ER_p;
A = A - diag(diag(A));
A = triu(A);
A = A + A';
General.save_var(A,obj.path,'adjacency_matrix');
end