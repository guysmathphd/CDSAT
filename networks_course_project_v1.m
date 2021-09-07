%%
%C:\16 - PhD\15 - courses tashpa 2\01 - networks\Project\yellow_taxi_trip_data
csvpath = fullfile('yellow_taxi_trip_data/','yellow_tripdata_2018-01.csv');
data = readtable(csvpath);
%%
weight_types = {'num_rides','trip_distance','fare_amount','tip_amount'};% 1 : for number of active rides
file_name_pre = 'yellow_tripdata_2018-01';
pu_locID = General.load_var(fullfile('source_data',file_name_pre,[file_name_pre '_filtered_pu_locID']));
do_locID = General.load_var(fullfile('source_data',file_name_pre,[file_name_pre '_filtered_do_locID']));
num_nodes = 265; A = zeros(num_nodes);
num_rows = size(pu_locID,1);
for weight_type = weight_types
    weight_type = weight_type{1};
    if isequal(weight_type,'num_rides')
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
    Network(['yellow taxi weighted-' weight_type],A,0,0,'yellow_tripdata_2018-01',weight_type);
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
