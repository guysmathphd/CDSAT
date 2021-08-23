function obj = plot_Pk(obj)
% fig net02*

name = 'net02a';legendStr = {};
figdesc = 'Degree distribution';
f = figure('Name',name,'NumberTitle','off');
Pk = General.load_var(fullfile(obj.path,'Pk'));
U = General.load_var(fullfile(obj.path,'k_sorted_unique'));
plot(U,Pk,'-o');
legendStr{end+1} = 'Emprical';
hold on;
if obj.ER_p
    lambda = obj.ER_p * obj.N;
    a = lambda.^U;
    b = factorial(U);
    c = exp(-lambda);
    Pk_poisson = c*a./b;
    % Pk_poisson = (exp(-lambda)*lambda.^U')/factorial(U');
    plot(U,Pk_poisson,'-');
    legendStr{end+1} = ['Poisson $\lambda = ' num2str(lambda) '$'];
end
xlabel('k');
ylabel('P(k)');
legend(legendStr,'Interpreter','latex');
title({[name ' ' obj.name];figdesc});
General.save_fig(f,name,fullfile(obj.path,'netfigs'));

end