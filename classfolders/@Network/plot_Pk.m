function obj = plot_Pk(obj)
% fig net02*
myplot = {@plot,@loglog};
numplots = length(myplot);
namesuf = {'','-loglog'};
for i1 = 1:numplots
    name = ['net02a' namesuf{i1}];legendStr = {};
    figdesc = 'Degree distribution';
    f = figure('Name',name,'NumberTitle','off');
    Pk = General.load_var(fullfile(obj.path,'Pk'));
    U = General.load_var(fullfile(obj.path,'k_sorted_unique'));
    myplot{i1}(U,Pk,'-o');
    legendStr{end+1} = 'Emprical';
    hold on;
    if obj.ER_p
        lambda = obj.ER_p * obj.N;
        a = lambda.^U;
        b = factorial(U);
        c = exp(-lambda);
        Pk_poisson = c*a./b;
        % Pk_poisson = (exp(-lambda)*lambda.^U')/factorial(U');
        myplot{i1}(U,Pk_poisson,'-');
        legendStr{end+1} = ['Poisson $\lambda = ' num2str(lambda) '$'];
    elseif obj.BA_m
        lambda = 3;
        Pk_th = U.^-lambda;
        myplot{i1}(U,Pk_th,'-');
        legendStr{end+1} = ['$k^{-\lambda}$, $\lambda = ' num2str(lambda) '$'];
    end
    xlabel('k');
    ylabel('P(k)');
    legend(legendStr,'Interpreter','latex');
    title({[name ' ' obj.name ' ' obj.desc];figdesc});
    General.save_fig(f,name,fullfile(obj.path,'netfigs'));
end


%% net02b P(k) binned
k = General.load_var(fullfile(obj.path,'degree_vector.mat'));
tol = 1e-10;numbins=15;values_vec = k;cond_vec = true(obj.N,1);
[ind_bins_var,bins_edges,ind_bins_var_sizes] = EngineClass.set_bins_generic(numbins,values_vec,tol,cond_vec);
bin_lengths = diff(bins_edges);
P_k = ind_bins_var_sizes./bin_lengths;
name = 'net02b';figdesc = 'Degree distribution Density Binned';
f = figure('Name',name,'NumberTitle','off');
loglog(bins_edges(2:end),P_k,'-o');
legend('Empirical');
xlabel('k');
ylabel('$P(k) = N_n/\mid B_n \mid$','Interpreter','latex');
title({[name ' ' obj.name ' ' obj.desc];figdesc});
General.save_fig(f,name,fullfile(obj.path,'netfigs'));

end





