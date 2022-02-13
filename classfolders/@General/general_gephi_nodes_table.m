function general_gephi_nodes_table(nodes,values,filepath,filename)

T = array2table([nodes,nodes,values],'VariableNames',{'ID','Label','value'});
writetable(T,fullfile(filepath,[filename '.csv']));

end