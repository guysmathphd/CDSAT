function obj = set_trip_data(obj)
csvpath = fullfile(obj.path,'yellow_tripdata_2018-01');
data = readtable(csvpath);
d1 = '2018-01-01 00:00:00';
t1 = datetime(d1,'InputFormat','yyyy-MM-dd HH:mm:ss');
d2 = '2018-01-31 23:59:59';
t2 = datetime(d2,'InputFormat','yyyy-MM-dd HH:mm:ss');

% clean data
t_pu = data.tpep_pickup_datetime;
t_do = data.tpep_dropoff_datetime;
hms_do = hms(t_do);
midnight = duration(0,0,0);
max_dur = hours(3);
durs = t_do - t_pu;
data(t_pu<t1 | t_pu>t2 | tod_do == midnight | durs>max_dur | t_pu>t_do,:) = [];
N = size(data,1);
% find day of week
dow = weekday(data.tpep_pickup_datetime);

% set number of active rides counter
pu_times = data.tpep_pickup_datetime;do_times = data.tpep_dropoff_datetime;
pu_counter = ones(size(pu_times));do_counter = -1*ones(size(do_times));
all_times = [pu_times;do_times];counter = [pu_counter; do_counter];
[B,I] = sort(all_times);counter_sorted = counter(I);
num_active_rides = cumsum(counter_sorted);

% set k_out and k_in over time
nodes_pu = data.PULocationID;nodes_do = data.DOLocationID;
nodes_pu_unique_sorted = sort(unique(nodes_pu));
num_nodes_out = size(nodes_pu_unique_sorted,1);
% k_out = zeros(size(counter,1),num_nodes_out);
loc_ids_out = [data.PULocationID;data.PULocationID];
loc_ids_out_sorted = loc_ids_out(I);
times_out = cell(1,num_nodes_out);
k_outs = cell(1,num_nodes_out);
for i1 = 1:num_nodes_out
    cur_loc = nodes_pu_unique_sorted(i1);
    cur_inds = find(loc_ids_out_sorted == cur_loc);
    times_out{1,i1} = B(cur_inds);
    counter_cur_loc = counter_sorted(cur_inds);
    k_outs{1,i1} = cumsum(counter_cur_loc);
end
figure;hold on;
for i1 = 1:num_nodes_out
    plot(times_out{1,i1},k_outs{1,i1});
end
nodes_do_unique_sorted = sort(unique(nodes_do));
num_nodes_in = size(nodes_do_unique_sorted,1);
loc_ids_in = [data.DOLocationID;data.DOLocationID];
loc_ids_in_sorted = loc_ids_in(I);
times_in = cell(1,num_nodes_in);
k_ins = cell(1,num_nodes_in);
for i1 = 1:num_nodes_in
    cur_loc = nodes_do_unique_sorted(i1);
    cur_inds = find(loc_ids_in_sorted == cur_loc);
    times_in{1,i1} = B(cur_inds);
    counter_cur_loc = counter_sorted(cur_inds);
    k_ins{1,i1} = cumsum(counter_cur_loc);
end
figure;hold on;
for i1 = 1:num_nodes_in
    plot(times_in{1,i1},k_ins{1,i1});
end
    
    
end




