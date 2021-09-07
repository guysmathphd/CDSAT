classdef NetworkSet
    %NETWORKSET Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        folderNames
        networkPaths % cell array with paths of all engine folders included in this engine set
        numNetworks
%         objPropertiesPaths
%         resultsPaths
        name
        path
        figs_path
        legendNames
    end
    
    methods
        function obj = NetworkSet(folderNames,name,legendNames)
            %NETWORKSET Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.folderNames = folderNames;
            obj.name = name;
            obj.legendNames = legendNames;
            obj.path = fullfile('networks',obj.name);
            if ~isfolder(obj.path)
                mkdir(obj.path)
            end
            obj.figs_path = fullfile(obj.path,'setfigs');
            if ~isfolder(obj.figs_path)
                mkdir(obj.figs_path);
            end
            obj.numNetworks = length(folderNames);
            for i = 1:obj.numNetworks
                obj.networkPaths{i} = fullfile('networks',obj.folderNames{i},[obj.folderNames{i} 'Obj']);
            end
            obj.save_obj();
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

