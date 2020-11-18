function A = make_circuit(n_vert)

A = make_path(n_vert);
A(1,n_vert) = 1;
A(n_vert,1) = 1;

end
