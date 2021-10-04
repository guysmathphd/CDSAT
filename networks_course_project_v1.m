%%
%C:\16 - PhD\15 - courses tashpa 2\01 - networks\Project\yellow_taxi_trip_data
csvpath = fullfile('yellow_taxi_trip_data/','yellow_tripdata_2018-01.csv');
data = readtable(csvpath);
%%
weight_types = {'num-rides','trip-distance','fare-amount','tip-amount'};% 1 : for number of active rides
average_weight_denominator_types = {'','num-rides','trip-distance','fare-amount'};
file_name_pre = 'yellow_tripdata_2018-01';
pu_locID = General.load_var(fullfile('source_data',file_name_pre,[file_name_pre '_filtered_pu_locID']));
do_locID = General.load_var(fullfile('source_data',file_name_pre,[file_name_pre '_filtered_do_locID']));
num_nodes = 265; A = zeros(num_nodes);
num_rows = size(pu_locID,1);ind = 1;
for weight_type = weight_types
    weight_type = weight_type{1};
    average_weight_denominator_type = average_weight_denominator_types{ind};
    ind = ind + 1;
    if isequal(weight_type,'num-rides')
        weights = ones(num_rows,1);
        General.save_var(weights,fullfile('source_data',file_name_pre),[file_name_pre '_filtered_' weight_type]);
    else
        weights = General.load_var(fullfile('source_data',file_name_pre,[file_name_pre '_filtered_' weight_type]));
    end
    for i = 1:size(do_locID,1)
        if mod(i,100000) == 0
            disp(['i = ' num2str(i)]);
        end
        A(pu_locID(i),do_locID(i)) = A(pu_locID(i),do_locID(i)) + weights(i);
    end
    Network(['yellow taxi weighted-' weight_type],A,0,0,'yellow_tripdata_2018-01',weight_type,...
    average_weight_denominator_type);
end

% obj1 = Network('yellow taxi unweighted',A>0);
%%
k_out = sum(A,2);
k_in = sum(A,1);

%%
figure;
histogram(k_out);

%%
G = graph(A,'upper');
%%
figure;
plot(G,'Layout','force');
%%
obj.set_trip_data();
obj.plot_degree_vs_time();
obj.plot_node_degrees_vs_time();
obj.plot_Pk();
obj.plot_degree();
obj.plot_degree2();
obj = plot_node_degrees_vs_time2(obj);


%%
folderNames = {'yellow taxi weighted-num-rides',...
    'yellow taxi weighted-trip-distance',...
    'yellow taxi weighted-fare-amount',...
    'yellow taxi weighted-tip-amount'};
name = 'yellow taxi weighted set';
legendNames = {'Number of rides','Average Trip Distance',...
    'Average fare amount',...
    'Average tip amount'};
obj = NetworkSet(folderNames,name,legendNames);

General.batchFunction(obj,@set_trip_data,[3:4]);
General.batchFunction(obj,@plot_degree,[3:4]);
General.batchFunction(obj,@plot_degree_vs_time,[3:4]);


%%

total_weight_sorted_min_to_max_inds = General.load_var(fullfile(obj.path,'total_weight_sorted_min_to_max_inds'));
all_times_sorted = General.load_var(fullfile(obj.path,'all_times_sorted'));

DTmin = all_times_sorted(total_weight_sorted_min_to_max_inds(end));

obj.write_gephi_edges_file(DTmin,'max');




