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
        save_obj(obj,path,name)
        save_fig(f,name,path)
        [p_0, p_1, yfit, rsq] = lin_reg(x_vals,y_vals)
    end
end

