classdef EngineSet < handle
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        engines
    end
    
    methods
        function obj = EngineSet(enginesPropertiesMaps)
            ind = 1;
            for p = enginesPropertiesMaps
            obj.engines{ind} = EngineClass(p{1});
            ind = ind+1;
            end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

