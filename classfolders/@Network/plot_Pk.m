function obj = plot_Pk(obj)
% fig net02*
myplot = {@plot,@loglog};
numplots = length(myplot);
namesuf = {'','-loglog'};
figdesc = 'Degree distribution';
Pk_out = General.load_var(fullfile(obj.path,'Pk'));
U_out = General.load_var(fullfile(obj.path,'k_sorted_unique'));
k_out = General.load_var(fullfile(obj.path,'degree_vector.mat'));
try
    Pk_in = General.load_var(fullfile(obj.path,'Pk_in'));
    U_in = General.load_var(fullfile(obj.path,'k_in_sorted_unique'));
    k_in = General.load_var(fullfile(obj.path,'degree_vector_in.mat'))';
    namesuf1 = {'k_{out}','k_{in}'};
    Pks = {Pk_out, Pk_in};Us = {U_out,U_in};
    ks = {k_out,k_in};
catch
    namesuf1 = {'k_{out}'};
    Pk = {Pk_out};Us = {U_out};
    ks = {k_out};
end

for i2 = 1:length(namesuf1)
    Pk = Pks{i2};U = Us{i2};
    for i1 = 1:numplots
        name = ['net02a' '-' namesuf1{i2} namesuf{i1}];legendStr = {};
        f = figure('Name',name,'NumberTitle','off');
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
        ylabel(['P(' namesuf1{i2} ')']);
        legend(legendStr,'Interpreter','latex');
        title({[name ' ' obj.name ' ' obj.desc];figdesc});
        General.save_fig(f,name,fullfile(obj.path,'netfigs'));
    end
    
    
    %% net02b P(k) binned
    
    tol = 1e-10;numbins=15;values_vec = ks{i2};cond_vec = true(obj.N,1);
    [ind_bins_var,bins_edges,ind_bins_var_sizes] = EngineClass.set_bins_generic(numbins,values_vec,tol,cond_vec);
    bin_lengths = diff(bins_edges);
    P_k = ind_bins_var_sizes./bin_lengths;
    name = ['net02b' '-' namesuf1{i2}];
    figdesc = 'Degree distribution Density Binned';
    f = figure('Name',name,'NumberTitle','off');
    x = bins_edges(2:end);
    loglog(x,P_k,'-o');hold on;
    loglog(x,x.^(-1));
    legend('Empirical', 'P(k) ~ k^{-1}');
    xlabel(namesuf1{i2});
    ylabel(['$P(' namesuf1{i2} ') = N_n/\mid B_n \mid$'],'Interpreter','latex');
    title({[name ' ' obj.name ' ' obj.desc];figdesc});
    General.save_fig(f,name,fullfile(obj.path,'netfigs'));
end
end





