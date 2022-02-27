function general_gephi_nodes_table(nodes,values,x,y,filepath,filename)
if isempty(x)
    x = zeros(size(values));
    y = x;
end
T = array2table([nodes,nodes,values,x,y],'VariableNames',{'ID','Label','value','x','y'});
writetable(T,fullfile(filepath,[filename '.csv']));

end