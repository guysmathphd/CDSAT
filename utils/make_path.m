function A = make_path(n_vert)

A = diag(ones(n_vert-1,1),1);
A = A + A';

end

