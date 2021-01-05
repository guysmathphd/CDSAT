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
                disp(['EngienSet ind = ' num2str(ind)]);
                p1 = p{1};
                p1('scenarioName') = [p1('scenarioName') '-' num2str(ind)];
                p1('resultsPath') = fullfile(p1('resultsPath'),p1('scenarioName'),'results');
            obj.engines{ind} = EngineClass(p1);
            ind = ind+1;
            end
        end
        function obj = setSolve(obj)
            for engcel = obj.engines
                eng = engcel{1};
                eng.solve(false,false,true,1,1,false)
            end
        end
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

