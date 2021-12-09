classdef EngineSet2
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        folderNames
        enginePaths % cell array with paths of all engine folders included in this engine set
        numEngines
        objPropertiesPaths
        resultsPaths
        name
        path
        figs_path
        legendNames
    end
    
    methods
        function obj = EngineSet2(folderNames,name,legendNames)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.folderNames = folderNames;
            obj.name = name;
            obj.legendNames = legendNames;
            obj.path = fullfile('tests',obj.name);
            if ~isfolder(obj.path)
                mkdir(obj.path)
            end
            obj.figs_path = fullfile(obj.path,'setfigs');
            if ~isfolder(obj.figs_path)
                mkdir(obj.figs_path);
            end
            obj.numEngines = length(folderNames);
            for i = 1:obj.numEngines
                obj.resultsPaths{i} = fullfile('tests',obj.folderNames{i},'results');
                obj.objPropertiesPaths{i} = fullfile(obj.resultsPaths{i},'obj_properties');
                obj.enginePaths{i} = fullfile(obj.resultsPaths{i},[obj.folderNames{i} 'Obj']);                
            end
            obj.save_obj();
        end
        obj = batchFunction(obj,function_handle,inds)
        obj = plotLocalization1(obj)
        obj = plotTauVsQ(obj)
        obj = save_obj(obj)
        obj = save_fig(obj,f,name)

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
        out1 = testfun1(obj,inputArg1)

        
    end
end

