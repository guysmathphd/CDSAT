function obj = plot_node_degrees_vs_time(obj)

times_out = General.load_var(fullfile(obj.path,'times_out'));
k_outs = General.load_var(fullfile(obj.path,'k_outs_vs_times_out'));
node_names_out = General.load_var(fullfile(obj.path,'node_names_out'));
times_in = General.load_var(fullfile(obj.path,'times_in'));
k_ins = General.load_var(fullfile(obj.path,'k_ins_vs_times_in'));
node_names_in = General.load_var(fullfile(obj.path,'node_names_in'));

num_nodes_out = length(node_names_out);
num_nodes_in = length(node_names_in);

name = 'net03a';figdesc = 'Number of active outgoing rides vs time';
f = figure('Name',name,'NumberTitle','off');hold on;
for i1 = 1:num_nodes_out
    plot(times_out{1,i1},k_outs{1,i1});
end
legend(node_names_out);
xlabel('t');ylabel('k_{out}');
title({[name ' ' obj.name ' ' obj.desc];figdesc});
General.save_fig(f,name,fullfile(obj.path,'netfigs'));

name = 'net03b';figdesc = 'Number of active incoming rides vs time';
f = figure('Name',name,'NumberTitle','off');hold on;
for i1 = 1:num_nodes_in
    plot(times_in{1,i1},k_ins{1,i1});
end
legend(node_names_in);
xlabel('t');ylabel('k_{in}');
title({[name ' ' obj.name ' ' obj.desc];figdesc});
General.save_fig(f,name,fullfile(obj.path,'netfigs'));


end