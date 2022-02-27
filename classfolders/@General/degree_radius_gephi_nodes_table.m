function degree_radius_gephi_nodes_table(nodes,values,degree_vector,filepath,filename)
k = degree_vector;
k_norm = k./min(k);k_norm = k_norm+1;
r = 1./log(k_norm);r = r - .5*min(r);
[s,i] = sort(r);p = floor(.99*length(i));
r(i(p:end)) = s(p);
% r = 1./k_norm;
rng(1);
thetas = rand(size(k))*2*pi;
x = r.*cos(thetas);
y = r.*sin(thetas);
General.general_gephi_nodes_table(nodes,values,x,y,filepath,filename);

end
