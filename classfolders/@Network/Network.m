classdef Network
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        path
        N
        desc
        ER_p
        BA_m
        source_data_name
        source_data_path
        weight_type
    end
    
    methods
        function obj = Network(name,adjacencyMatrix,type,ER_p_BA_m,source_data_name,weight_type)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.name = name;
            obj.source_data_name = source_data_name;
            obj.weight_type = weight_type;
            obj.source_data_path = fullfile('source_data',source_data_name);
            obj.N = size(adjacencyMatrix,1);
            obj.path = fullfile('networks',name);
            if ~isfolder(obj.path)
                mkdir(obj.path);
            end
            switch type
                case 0 % from file, do nothing
                    General.save_var(adjacencyMatrix,obj.path,'adjacency_matrix');
                case 1 % ER with p = ER_p_BA_m
                    obj.ER_p = ER_p_BA_m;
                    obj.create_ER();
                case 2 % BA with m = ER_P_BA_m
                    obj.BA_m = ER_p_BA_m;
                    obj.create_BA();
            end
            
            obj = obj.set_graph_props();
            General.save_obj(obj,obj.path,name);
        end
        obj = set_graph_props(obj);
        obj = set_degree_perts(obj);
        obj = plotNet01(obj);
        obj = set_random_samples(obj);
        obj = create_ER(obj);
        obj = create_BA(obj);
        obj = plot_Pk(obj);
        obj = plot_degree(obj);
        obj = plot_degree_vs_time(obj);
        obj = write_gephi_edges_file(obj);
        obj = calc_geo_distances(obj);
        obj = set_trip_data(obj);
        obj = clean_trip_data(obj);
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

