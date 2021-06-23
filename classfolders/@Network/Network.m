classdef Network
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        path
    end
    
    methods
        function obj = Network(name,adjacencyMatrix)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.name = name;
            obj.path = fullfile('networks',name);
            if ~isfolder(obj.path)
                mkdir(obj.path);
            end
            General.save_var(adjacencyMatrix,obj.path,'adjacency_matrix');
            obj.set_graph_props();
            General.save_obj(obj,obj.path,name);
        end
        obj = set_graph_props(obj);
        obj = set_degree_perts(obj);
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

