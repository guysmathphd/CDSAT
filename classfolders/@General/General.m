classdef General
    %GENERAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        colors1 = [0 0.4470 0.7410;
            0.8500 0.3250 0.0980;
            0.9290 0.6940 0.1250;
            0.4940 0.1840 0.5560;
            0.4660 0.6740 0.1880;
            0.3010 0.7450 0.9330;
            0.6350 0.0780 0.1840];
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
        batchFunction(set_obj,function_handle,inds)
        general_gephi_nodes_table(nodes,values,x,y,filepath,filename);
        degree_radius_gephi_nodes_table(nodes,values,degree_vector,filepath,filename);
        make_fig_bold(folderpath);
        f=ThreeDim_network_layout_visualization(d,v,n,fname);
    end
end

