function plotNet03(obj)
subnetsize = 1000;
fname = 'net10a';
k = General.load_var(fullfile(obj.path, 'degree_vector.mat'));
[K,I] = sort(k,'descend');
k_sub = k(I(1:subnetsize));
f = General.ThreeDim_network_layout_visualization(k_sub,zeros(size(k_sub)),subnetsize,fname);
General.save_fig(f,fname,fullfile(obj.path,'figs'));

end