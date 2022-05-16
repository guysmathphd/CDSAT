function obj = plotNet01(obj)
k_out = General.load_var(fullfile(obj.path,'degree_vector'));
density =[];
% % try
% % density = General.load_var(fullfile(obj.path,'density'));
% % catch
% end
k_out_avg = sum(k_out)/obj.N;
issymmetric = true;
try
    k_in = General.load_var(fullfile(obj.path,'degree_vector_in'));
    issymmetric = false;
    k_in_avg = sum(k_in)/obj.N;
catch
end
k_out_prop = k_out;
%k_in_prop = k_in;

% name = 'Net01a';
% f = figure('Name',name,'NumberTitle','off');
% histogram(k_out_prop);
% 
% name = 'Net01b';
% f = figure('Name',name,'NumberTitle','off');
% histogram(k_in_prop);

% numbins = 15;values_vec = k_out;tol = 1e-13;cond_vec = true(obj.N,1);
% [ind_bins_var,bins_edges,ind_bins_var_sizes] = EngineClass.set_bins_generic(numbins,values_vec,tol,cond_vec);

num_bins = 20;values_vec = k_out;
[ind_bins_var,bins_edges,ind_bins_var_sizes] = EngineClass.set_bins_nonlog(num_bins,values_vec);
name = 'Net01c';figdesc = 'Degree distributions';
figdesc1 = ['<k_{out}> = ' num2str(k_out_avg) ];
if ~isempty(density)
    figdesc1 = [figdesc1 ', density = ' num2str(density)];
end
f = figure('Name',name,'NumberTitle','off');legendStr = {};
loglog(bins_edges(2:end),ind_bins_var_sizes/obj.N,'o-');
legendStr{end+1} = 'k';hold on;
if ~issymmetric
    values_vec = k_in;
    [ind_bins_var,bins_edges,ind_bins_var_sizes] = EngineClass.set_bins_nonlog(num_bins,values_vec);
    loglog(bins_edges(2:end),ind_bins_var_sizes/obj.N,'*-');
    legendStr{end} = 'k_{out}';legendStr{end+1} = 'k_{in}';
end
legend(legendStr);
xlabel('Degree');ylabel('P(k)','interpreter','latex');
title({[name ' ' obj.name];figdesc;figdesc1});
General.save_fig(f,name,fullfile(obj.path,'netfigs'));
end