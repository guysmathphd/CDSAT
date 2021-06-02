classdef General
    %GENERAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        function obj = General()
            %GENERAL Construct an instance of this class
            %   Detailed explanation goes here
        end
    end
    methods (Static)
        save_var(var,path,name)
        var = load_var(path)
        [ bytes ] = getSize( variable )
    end
end

