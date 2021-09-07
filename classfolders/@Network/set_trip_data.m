function obj = set_trip_data(obj)
csvpath_nodes = fullfile(obj.source_data_path,'nodes');
data_nodes = readtable(csvpath_nodes);



% N = size(data,1);
% find day of week
% dow = weekday(data.tpep_pickup_datetime);

% set number of active rides counter
pu_times = General.load_var(fullfile(obj.source_data_path,[obj.source_data_name '_filtered_pu_datetime']));
do_times = General.load_var(fullfile(obj.source_data_path,[obj.source_data_name '_filtered_do_datetime']));

weight = General.load_var(fullfile(obj.source_data_path,[obj.source_data_name '_filtered_' obj.weight_type]));
pu_counter = weight.*ones(size(pu_times));do_counter = -1*weight.*ones(size(do_times));
all_times = [pu_times;do_times];[all_times_sorted,all_times_sorted_inds] = sort(all_times);

counter = [pu_counter; do_counter];
counter_sorted = counter(all_times_sorted_inds);
total_weight_vs_all_times_sorted = cumsum(counter_sorted);
General.save_var(all_times_sorted,obj.path,'all_times_sorted');
General.save_var(total_weight_vs_all_times_sorted,obj.path,'total_weight_vs_all_times_sorted');

% set k_out and k_in over time
nodes_pu =  General.load_var(fullfile(obj.source_data_path,[obj.source_data_name '_filtered_pu_locID']));
nodes_do = General.load_var(fullfile(obj.source_data_path,[obj.source_data_name '_filtered_do_locID']));
nodes_pu_unique_sorted = sort(unique(nodes_pu));
num_nodes_out = size(nodes_pu_unique_sorted,1);
% k_out = zeros(size(counter,1),num_nodes_out);
loc_ids_out = [nodes_pu;nodes_pu];
loc_ids_out_sorted = loc_ids_out(all_times_sorted_inds);
times_out = cell(1,num_nodes_out);
k_outs = cell(1,num_nodes_out);
node_names_out = cell(1,num_nodes_out);
for i1 = 1:num_nodes_out
    cur_loc = nodes_pu_unique_sorted(i1);
    node_names_out{i1} = data_nodes.Name{data_nodes.ID == cur_loc};
    cur_inds = find(loc_ids_out_sorted == cur_loc);
    times_out{1,i1} = all_times_sorted(cur_inds);
    counter_cur_loc = counter_sorted(cur_inds);
    k_outs{1,i1} = cumsum(counter_cur_loc);
end
General.save_var(k_outs,obj.path,'k_outs_vs_times_out');
General.save_var(times_out,obj.path,'times_out');
General.save_var(node_names_out,obj.path,'node_names_out');


nodes_do_unique_sorted = sort(unique(nodes_do));
num_nodes_in = size(nodes_do_unique_sorted,1);
loc_ids_in = [nodes_do;nodes_do];
loc_ids_in_sorted = loc_ids_in(all_times_sorted_inds);
times_in = cell(1,num_nodes_in);
k_ins = cell(1,num_nodes_in);
node_names_in = cell(1,num_nodes_in);
for i1 = 1:num_nodes_in
    cur_loc = nodes_do_unique_sorted(i1);
    node_names_in{i1} = data_nodes.Name{data_nodes.ID == cur_loc};
    cur_inds = find(loc_ids_in_sorted == cur_loc);
    times_in{1,i1} = all_times_sorted(cur_inds);
    counter_cur_loc = counter_sorted(cur_inds);
    k_ins{1,i1} = cumsum(counter_cur_loc);
end
General.save_var(k_ins,obj.path,'k_ins_vs_times_in');
General.save_var(times_in,obj.path,'times_in');
General.save_var(node_names_in,obj.path,'node_names_in');

name = 'net03a';figdesc = 'Number of active outgoing rides vs time';
f = figure('Name',name,'NumberTitle','off');hold on;
for i1 = 1:num_nodes_out
    plot(times_out{1,i1},k_outs{1,i1});
end
legend(node_names_out);
xlabel('t');ylabel('k_out');
title({[name ' ' obj.name ' ' obj.desc];figdesc});
General.save_fig(f,name,fullfile(obj.path,'netfigs'));



name = 'net03b';figdesc = 'Number of active incoming rides vs time';
f = figure('Name',name,'NumberTitle','off');hold on;
for i1 = 1:num_nodes_in
    plot(times_in{1,i1},k_ins{1,i1});
end
legend(node_names_in);
xlabel('t');ylabel('k_in');
title({[name ' ' obj.name ' ' obj.desc];figdesc});
General.save_fig(f,name,fullfile(obj.path,'netfigs'));
    
end




