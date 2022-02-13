function obj = plot_eigvecs_mass(obj)

M = General.load_var(fullfile(obj.path,'adj_mat_eigenvectors.mat'));
k = General.load_var(fullfile(obj.path,'degree_vector.mat'));
[~,ind] = sort(k,'descend');
M_sorted = M(ind,:);
name = 'net08a';namef = 'net08a';path = fullfile(obj.path);
mytitle = 'Adjacency Matrix Eigenvectors Absolute Value';
xlabelstr = '$v^{(i)}$';
ylabelstr = 'Nodes (sorted by degree)';
EngineClass.plot_image_static(abs(M_sorted),name,namef,path,mytitle,xlabelstr,ylabelstr);

sums = sum(abs(M_sorted));
mass_prop = abs(M_sorted)./sums;
name = 'net08b';namef = 'net08b';
mytitle = 'Adjacency Matrix Eigenvectors Mass Proportion';
EngineClass.plot_image_static(mass_prop,name,namef,path,mytitle,xlabelstr,ylabelstr);

pev_mass_prop = mass_prop(:,1);
name = 'net08c';namef = 'net08b';
mytitle = 'Adjacency Matrix PEV Mass Proportion';
EngineClass.plot_image_static(pev_mass_prop,name,namef,path,mytitle,xlabelstr,ylabelstr);

end