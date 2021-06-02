function obj = plotLocalization1(obj)
vars = cell(1,obj.numEngines); means = cell(1,obj.numEngines); sol_ts = cell(1,obj.numEngines);
mydatadist = load(fullfile('networks',['SF1shortest_distances']));

for i=1:obj.numEngines
    disp(['plotLocalization i = ' num2str(i)]);
    mydata = load(fullfile(obj.objPropertiesPaths{i},'eigvec_pert_max_hub','sol_t_hub'));
    sol_ts{i} = mydata.var; clear mydata;
    mydata = load(fullfile(obj.objPropertiesPaths{i},'eigvec_pert_max_hub','sol_x_hub'));
    sol_x = mydata.var; clear mydata;
    mydata = load(obj.enginePaths{i});eng = mydata.obj;clear mydata;
    ss = eng.steady_state;
    dv = eng.degree_vector_weighted;
    [~,j] = max(dv);
    dist = mydatadist.var(j,:);
    pert_x = (ss' - sol_x')';clear sol_x;
    [means{i}, vars{i}, ~] = EngineClass.compute_wmeans_wvars(abs(pert_x'), dist');
end

name = ['set3a'];
figdesc = 'Perturbation Weighted Std, ';
figdesc2 = ['$\hat{\sigma_d} = \left[(\sum_i \hat{pert_i} (d_i - \hat{\mu_d})^2)\right]^{1/2}$, $\hat{pert_i} = \mid pert_i \mid/S$, $S = \sum_i \mid pert_i \mid$,'];
figdesc3 = [' $\hat{\mu_d} = \sum_i \hat{pert_i} d_i $,'];
figdesc4 = ['$pert_{t=0}$ is all mass at biggest hub'];
figdesc5 = ' $d_i$ is distance from biggest hub to node $i$';
f = figure('Name',name,'NumberTitle','off');
for i = 1:obj.numEngines
    plot(sol_ts{i},vars{i});
    hold on;
end
xlabel('t');ylabel('$\hat{\sigma}_d$','interpreter','latex');
title({[name ' ' obj.name];[figdesc figdesc4 ];[figdesc2 figdesc3 figdesc5]},'interpreter','latex');
legend(obj.legendNames);
obj.save_fig(f,name);
end
