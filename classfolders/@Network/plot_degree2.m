function obj = plot_degree2(obj)
% net06*
total_weight_sorted_min_to_max = General.load_var(fullfile(obj.path,'total_weight_sorted_min_to_max'));
total_weight_sorted_min_to_max_inds = General.load_var(fullfile(obj.path,'total_weight_sorted_min_to_max_inds'));

all_times_sorted = General.load_var(fullfile(obj.path,'all_times_sorted'));

max_weight_time = all_times_sorted(total_weight_sorted_min_to_max_inds(end));
min_weight_time = all_times_sorted(total_weight_sorted_min_to_max_inds(1));

times_out = General.load_var(fullfile(obj.path,'times_out'));

k_outs_vs_times_out = General.load_var(fullfile(obj.path,'k_outs_vs_times_out'));

nodes_pu_unique_sorted = General.load_var(fullfile(obj.path, 'nodes_pu_unique_sorted'));
node_names_out = General.load_var(fullfile(obj.path,'node_names_out'));
N = length(nodes_pu_unique_sorted);
k_out_max_weight_time = zeros(1,N);time_diffs_max = zeros(1,N);
k_out_min_weight_time = zeros(1,N);time_diffs_min = zeros(1,N);
for i1 = 1:N
    cur_times = times_out{1,i1};
    [time_diff,ind] = min(abs(cur_times - max_weight_time));
    if time_diff > 1800
        continue
    end
    time_diffs_max(1,i1) = seconds(time_diff);
    k_out_max_weight_time(1,i1) = k_outs_vs_times_out{1,i1}(ind);
    [time_diff,ind] = min(abs(cur_times - min_weight_time));
    time_diffs_min(1,i1) = seconds(time_diff);
    k_out_min_weight_time(1,i1) = k_outs_vs_times_out{1,i1}(ind);
end

[k_out_max_weight_time_sorted,inds_max] = sort(k_out_max_weight_time,'descend');
[k_out_min_weight_time_sorted,inds_min] = sort(k_out_min_weight_time,'descend');

node_names_out_sorted = node_names_out(inds_max);
X = categorical(node_names_out_sorted);
X = reordercats(X,node_names_out_sorted);
Y = k_out_max_weight_time_sorted;
name = 'net06a';figdesc1 = ['k_{out} when k_{tot} is greatest, at t = ' datestr(max_weight_time)];
f = figure('Name',name,'NumberTitle','off');
bar(X,Y);
xlabel('Location');
ylabel('k_i(t_{k_{totmax}}) ');
title({[name ' ' obj.name ' ' obj.desc];[figdesc1]});
General.save_fig(f,name,fullfile(obj.path,'netfigs'));


node_names_out_sorted = node_names_out(inds_min);
X = categorical(node_names_out_sorted);
X = reordercats(X,node_names_out_sorted);
Y = k_out_min_weight_time_sorted;
name = 'net06b';figdesc1 = ['k_{out} when k_{tot} is least, at t = ' datestr(min_weight_time)];
f = figure('Name',name,'NumberTitle','off');
bar(X,Y);
xlabel('Location');
ylabel('k_i(t_{k_{totmin}}) ');
title({[name ' ' obj.name ' ' obj.desc];[figdesc1]});
General.save_fig(f,name,fullfile(obj.path,'netfigs'));
end