function plotTauVsQ(obj)
%set4*
foldernames = {'eigvec_pert_max_hub','eigvec_pert_min_hub',...
    'solution_random_perts','single_node_pert_sol'};
colors = {'b','r','k','m','g','c','y'};
markers = {'o','s','^','*','+','x','v','<','>','.','d','p','h','_','|'};
suf = {'','-loglog'};plotfun = {@plot,@loglog};
for j = 1:length(suf)
    name = ['set4a' suf{j}]; fname{1} = name;desc = '$\tau_A$ vs Q4';%\tau_A vs Q4
    f{1} = figure('Name',name,'NumberTitle','off');hax{1} = axes; %hold on;
    tit1 = {[name ' ' obj.name];desc};
    
    name = ['set4b' suf{j}]; fname{1} = name;desc = '$\tau_A$ vs Q6';%\tau_A vs Q6
    f{2} = figure('Name',name,'NumberTitle','off');hax{2} = axes; %hold on;
    tit2 = {[name ' ' obj.name];desc};
%     title({[name ' ' obj.name];desc},'interpreter','latex');
    
    name = ['set4c' suf{j}]; fname{3} = name;desc = '$\tau_A$ vs Q7';%\tau_A vs Q7
    f{3} = figure('Name',name,'NumberTitle','off');hax{3} = axes;% hold on;
    tit3 = {[name ' ' obj.name];desc};
        
    name = ['set4d' suf{j}]; fname{4} = name;desc = '$\tau_A$ vs Q8';%\tau_A vs Q8
    f{4} = figure('Name',name,'NumberTitle','off');hax{4} = axes;% hold on;
    tit4 = {[name ' ' obj.name];desc};
%     title({[name ' ' obj.name];desc},'interpreter','latex');
    
    name = ['set4e' suf{j}]; fname{5} = name;desc = '$\tau_A$ vs Q9';%\tau_A vs Q9
    f{5} = figure('Name',name,'NumberTitle','off');hax{5} = axes;% hold on;
    tit5 = {[name ' ' obj.name];desc};

    legendObj = [];legendStr = {};
    
    for i1=1:obj.numEngines
        disp(['plotTauVsQ i = ' num2str(i1)]);
        resultspath = obj.resultsPaths{i1};
        ind = 1;size=6;
        for i2 = 1:length(foldernames)
            folder = foldernames{i2};%folder = foldernames
            folderpath = fullfile(resultspath,'obj_properties/',folder);
            files = dir(folderpath);
            
            for i3 = 1:length(files)
                file = files(i3);
                name = file.name;
                ind1 = strfind(name,'_');
                if length(ind1)>1
                    ind1 = ind1(2);
                else
                    continue
                end
                ind2 = strfind(name,'_norms');
                ind3 = strfind(name,'_norms_normed');
                ind4 = strfind(name,'_thetas');
                if ~isempty(ind2) && isempty(ind3)
                    indsuf = ind2;
                    
                    yvalues = General.load_var(fullfile(folderpath,name));
                    tvalues = General.load_var(fullfile(folderpath,[name(1:ind1-2) 't' name(ind1:indsuf-1)]));
                    tau_A = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_tau_A']));
                    Q4 = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_Q4']));
                    p = plotfun{j}(hax{1},Q4,tau_A,colors{i1},'Marker',markers{i1},'MarkerSize',size);hold(hax{1},'on');
                    title(hax{1},tit1,'interpreter','latex');
                    xlabel(hax{1},'Q4','Interpreter','latex');ylabel(hax{1},['$\tau_A$' suf{j}],'Interpreter','latex');
                    
                    Q6 = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_Q6']));
                    p = plotfun{j}(hax{2},Q6,tau_A,colors{i1},'Marker',markers{i1},'MarkerSize',size);hold(hax{2},'on');
                    title(hax{2},tit2,'interpreter','latex');
                    xlabel(hax{2},'Q6','Interpreter','latex');ylabel(hax{2},['$\tau_A$' suf{j}],'Interpreter','latex');
                    
                    Q7 = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_Q7']));
                    plotfun{j}(hax{3},Q7,tau_A,colors{i1},'Marker',markers{i1},'MarkerSize',size);hold(hax{3},'on');
                    title(hax{3},tit3,'interpreter','latex');
                    xlabel(hax{3},'$Q7 = \frac{\mathbf{\delta x(0)}\cdot k^{\alpha}}{\| \mathbf{\delta x(0)} \|_1 \max{(k)}}$','Interpreter','latex');ylabel(hax{3},['$\tau_A$' suf{j}],'Interpreter','latex');

                    Q8 = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_Q8']));
                    plotfun{j}(hax{4},Q8,tau_A,colors{i1},'Marker',markers{i1},'MarkerSize',size);hold(hax{4},'on');
                    title(hax{4},tit4,'interpreter','latex');
                    xlabel(hax{4},'$Q8 = \frac{\mathbf{\delta x(0)}\cdot k^{\alpha}}{\| \mathbf{\delta x(0)} \|_1 }$','Interpreter','latex');ylabel(hax{4},['$\tau_A$' suf{j}],'Interpreter','latex');
                    
                    Q9 = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_Q9']));
                    plotfun{j}(hax{5},Q9,tau_A,colors{i1},'Marker',markers{i1},'MarkerSize',size);hold(hax{5},'on');
                    title(hax{5},tit5,'interpreter','latex');
                    xlabel(hax{5},'$Q8 = \frac{\mathbf{\delta x(0)}\cdot k^{-\mu} \cdot \mathbf{ss}}{\| \mathbf{\delta x(0)} \|_1 }$','Interpreter','latex');ylabel(hax{5},['$\tau_A$' suf{j}],'Interpreter','latex');

                    if ind == 1
                        legendObj(end+1) = p;
                        legendStr{end+1} = [obj.legendNames{i1}];
                        %                     size = size + 4;
                    end
                    ind = ind+1;
                end
            end
        end
        
    end
    for i=1:length(f)
        legend(hax{i},legendObj,legendStr);
        General.save_fig(f{i},fname{i},obj.figs_path);
    end
end
end