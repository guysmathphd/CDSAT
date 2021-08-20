classdef Network
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        path
        N
        desc
        ER_p
    end
    
    methods
        function obj = Network(name,adjacencyMatrix,type,ER_p)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.name = name;
            obj.N = size(adjacencyMatrix,1);
            obj.path = fullfile('networks',name);
            if ~isfolder(obj.path)
                mkdir(obj.path);
            end
            switch type
                case 0 % from file, do nothing
                    General.save_var(adjacencyMatrix,obj.path,'adjacency_matrix');
                case 1 % ER with p = ER_p
                    obj.ER_p = ER_p;
                    obj.create_ER();
            end
            
            obj = obj.set_graph_props();
            General.save_obj(obj,obj.path,name);
        end
        obj = set_graph_props(obj);
        obj = set_degree_perts(obj);
        obj = plotNet01(obj);
        obj = set_random_samples(obj);
        obj = create_ER(ojb);
        obj = plot_Pk(obj);
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

