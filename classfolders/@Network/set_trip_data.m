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
pu_counter_unweighted = ones(size(pu_times));do_counter_unweighted = -1*ones(size(do_times));
all_times = [pu_times;do_times];[all_times_sorted,all_times_sorted_inds] = sort(all_times);
d2 = '2018-01-31 23:59:59'; %last pick up cutoff
t2 = datetime(d2,'InputFormat','yyyy-MM-dd HH:mm:ss');
d1 = '2018-01-01 01:22:30'; %first peak
t1 = datetime(d1,'InputFormat','yyyy-MM-dd HH:mm:ss');
max_time_ind = find(all_times_sorted >= t2,1);
min_time_ind = find(all_times_sorted >= t1,1);



counter = [pu_counter; do_counter];
counter_unweighted = [pu_counter_unweighted; do_counter_unweighted];
counter_sorted = counter(all_times_sorted_inds);
counter_sorted_unweighted = counter_unweighted(all_times_sorted_inds);

total_weight_vs_all_times_sorted = cumsum(counter_sorted);
% total_rides_vs_all_times_sorted = cumsum(counter_sorted_unweighted);
% average_weight_vs_all_times_sorted = total_weight_vs_all_times_sorted./total_rides_vs_all_times_sorted;
if obj.average_weight_denominator_type
    average_weight_denominator = General.load_var(fullfile('networks',...
        ['yellow taxi weighted-' obj.average_weight_denominator_type],...
        'total_weight_vs_all_times_sorted'));
    average_k_out_denominator = General.load_var(fullfile('networks',...
        ['yellow taxi weighted-' obj.average_weight_denominator_type],...
        'k_outs_vs_times_out'));
else
    average_weight_denominator = [];
    average_k_out_denominator = [];
end

all_times_sorted_cut = all_times_sorted(min_time_ind:max_time_ind);
all_times_sorted_inds_cut = all_times_sorted_inds(min_time_ind:max_time_ind);
total_weight_vs_all_times_sorted_cut = total_weight_vs_all_times_sorted(min_time_ind:max_time_ind);
% total_rides_vs_all_times_sorted_cut = total_rides_vs_all_times_sorted(min_time_ind:max_time_ind);
if average_weight_denominator
    average_weight_vs_all_times_sorted_cut = total_weight_vs_all_times_sorted_cut./average_weight_denominator;
%     average_weight_vs_all_times_sorted_cut = average_weight_vs_all_times_sorted(min_time_ind:max_time_ind);
end

General.save_var(all_times_sorted_cut,obj.path,'all_times_sorted');
General.save_var(total_weight_vs_all_times_sorted_cut,obj.path,'total_weight_vs_all_times_sorted');
if average_weight_denominator
    General.save_var(average_weight_vs_all_times_sorted_cut,obj.path,'average_weight_vs_all_times_sorted');
end
% General.save_var(total_rides_vs_all_times_sorted_cut,obj.path,'total_rides_vs_all_times_sorted_cut');

[M,I] = sort(total_weight_vs_all_times_sorted_cut);
General.save_var(M,obj.path,'total_weight_sorted_min_to_max');
General.save_var(I,obj.path,'total_weight_sorted_min_to_max_inds');



% set k_out and k_in per node over time
nodes_pu =  General.load_var(fullfile(obj.source_data_path,[obj.source_data_name '_filtered_pu_locID']));
nodes_do = General.load_var(fullfile(obj.source_data_path,[obj.source_data_name '_filtered_do_locID']));
nodes_pu_unique_sorted = sort(unique(nodes_pu));
num_nodes_out = size(nodes_pu_unique_sorted,1);
% k_out = zeros(size(counter,1),num_nodes_out);
loc_ids_out = [nodes_pu;nodes_pu];
loc_ids_out_sorted = loc_ids_out(all_times_sorted_inds);
times_out = cell(1,num_nodes_out);
k_outs = cell(1,num_nodes_out);
k_outs_num_rides = cell(1,num_nodes_out);
k_outs_average = cell(1,num_nodes_out);
node_names_out = cell(1,num_nodes_out);
max_k_out = zeros(1,num_nodes_out);max_k_out_inds = zeros(1,num_nodes_out);
min_k_out = zeros(1,num_nodes_out);min_k_out_inds = zeros(1,num_nodes_out);
for i1 = 1:num_nodes_out
 %   disp(num2str(i1));
    cur_loc = nodes_pu_unique_sorted(i1);
    node_names_out{1,i1} = data_nodes.Name{data_nodes.ID == cur_loc};
    cur_inds = find(loc_ids_out_sorted == cur_loc);
    times_out{1,i1} = all_times_sorted(cur_inds);
    counter_cur_loc = counter_sorted(cur_inds);
    counter_unweighted_cur_loc = counter_sorted_unweighted(cur_inds);
    k_outs{1,i1} = cumsum(counter_cur_loc);
    %     k_outs_num_rides{1,i1} = cumsum(counter_unweighted_cur_loc);

    max_time_ind = find(times_out{1,i1} <= t2,1,'last');
    min_time_ind = find(times_out{1,i1} >= t1,1,'first');
    times_out{1,i1} = times_out{1,i1}(min_time_ind:max_time_ind);
    k_outs{1,i1} = k_outs{1,i1}(min_time_ind:max_time_ind);
    if ~isempty(average_k_out_denominator)
        k_outs_average{1,i1} = k_outs{1,i1}./average_k_out_denominator{1,i1};
    end
    [max_k_out(1,i1),max_k_out_inds(1,i1)] = max(k_outs{1,i1});
    [min_k_out(1,i1),min_k_out_inds(1,i1)] = min(k_outs{1,i1});
end
General.save_var(nodes_pu_unique_sorted,obj.path,'nodes_pu_unique_sorted');
General.save_var(k_outs,obj.path,'k_outs_vs_times_out');
if ~isempty(average_k_out_denominator)
    General.save_var(k_outs_average,obj.path,'k_outs_average');
end
General.save_var(times_out,obj.path,'times_out');
General.save_var(node_names_out,obj.path,'node_names_out');
General.save_var(max_k_out,obj.path,'max_k_out');
General.save_var(max_k_out_inds,obj.path,'max_k_out_inds');
General.save_var(min_k_out,obj.path,'min_k_out');
General.save_var(min_k_out_inds,obj.path,'min_k_out_inds');


nodes_do_unique_sorted = sort(unique(nodes_do));
num_nodes_in = size(nodes_do_unique_sorted,1);
loc_ids_in = [nodes_do;nodes_do];
loc_ids_in_sorted = loc_ids_in(all_times_sorted_inds);
times_in = cell(1,num_nodes_in);
k_ins = cell(1,num_nodes_in);
node_names_in = cell(1,num_nodes_in);
min_k_in = zeros(1,num_nodes_in);min_k_in_inds = zeros(1,num_nodes_in);

for i1 = 1:num_nodes_in
    cur_loc = nodes_do_unique_sorted(i1);
    node_names_in{i1} = data_nodes.Name{data_nodes.ID == cur_loc};
    cur_inds = find(loc_ids_in_sorted == cur_loc);
    times_in{1,i1} = all_times_sorted(cur_inds);
    counter_cur_loc = counter_sorted(cur_inds);
    k_ins{1,i1} = cumsum(counter_cur_loc);
    [max_k_in(1,i1),max_k_in_inds(1,i1)] = max(k_ins{1,i1});
    [min_k_in(1,i1),min_k_in_inds(1,i1)] = min(k_ins{1,i1});
end
General.save_var(nodes_do_unique_sorted,obj.path,'nodes_do_unique_sorted');
General.save_var(k_ins,obj.path,'k_ins_vs_times_in');
General.save_var(times_in,obj.path,'times_in');
General.save_var(node_names_in,obj.path,'node_names_in');
General.save_var(max_k_in,obj.path,'max_k_in');
General.save_var(max_k_in_inds,obj.path,'max_k_in_inds');
General.save_var(min_k_in,obj.path,'min_k_in');
General.save_var(min_k_in_inds,obj.path,'min_k_in_inds');
end




