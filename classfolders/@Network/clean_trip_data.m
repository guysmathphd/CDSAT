function obj = clean_trip_data(obj)
csvpath = fullfile(obj.path,'yellow_tripdata_2018-01');
data = readtable(csvpath);
d1 = '2018-01-01 00:00:00'; % first pick up cutoff
t1 = datetime(d1,'InputFormat','yyyy-MM-dd HH:mm:ss');
d2 = '2018-01-31 23:59:59'; %last pick up cutoff
t2 = datetime(d2,'InputFormat','yyyy-MM-dd HH:mm:ss');

% clean data
t_pu = data.tpep_pickup_datetime;
t_do = data.tpep_dropoff_datetime;
tod_do = timeofday(t_do);
midnight = duration(0,0,0);
max_dur = hours(3);
durs = t_do - t_pu;
data(t_pu<t1 | t_pu>t2 | tod_do == midnight | durs>max_dur | t_pu>t_do,:) = [];
General.save_var(data.tpep_pickup_datetime,obj.path,'yellow_tripdata_2018-01_filtered_pu_datetime');
General.save_var(data.tpep_dropoff_datetime ,obj.path,'yellow_tripdata_2018-01_filtered_do_datetime');
General.save_var(data.PULocationID,obj.path,'yellow_tripdata_2018-01_filtered_pu_locID');
General.save_var(data.DOLocationID,obj.path,'yellow_tripdata_2018-01_filtered_do_locID');
General.save_var(data.trip_distance,obj.path,'yellow_tripdata_2018-01_filtered_trip_distance');
General.save_var(data.fare_amount,obj.path,'yellow_tripdata_2018-01_filtered_fare_amount');
General.save_var(data.tip_amount,obj.path,'yellow_tripdata_2018-01_filtered_tip_amount');

end