function plotTauVsQ(obj)
%set4*
foldernames = {'eigvec_pert_max_hub','eigvec_pert_min_hub',...
    'solution_random_perts'};
colors = {'b','r','k','m','g','c'};
markers = {'o','s','^','*','+','x','v','<','>','.','d','p','h','_','|'};
name = 'set4a'; fname{1} = name;desc = '$\tau_A$ vs Q';%\tau_A vs Q
f{1} = figure('Name',name,'NumberTitle','off');hax{1} = axes; hold on;
title({[name ' ' obj.name];desc},'interpreter','latex');
xlabel('Q','Interpreter','latex');ylabel('$\tau_A$','Interpreter','latex');
legendStr = {};

for i1=1:obj.numEngines
    disp(['plotTauVsQ i = ' num2str(i1)]);
    resultspath = obj.resultsPaths{i1};
    ind = 1;
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
                legendStr{end+1} = [obj.legendNames{i1} ' ' name(ind1+1:indsuf-1)];
                yvalues = General.load_var(fullfile(folderpath,name));
                tvalues = General.load_var(fullfile(folderpath,[name(1:ind1-2) 't' name(ind1:indsuf-1)]));
                Q = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_Q']));
                tau_A = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_tau_A']));
                plot(hax{1},Q,tau_A,colors{i1},'Marker',markers{ind});
                ind = ind+1;
            end
        end
    end
    
end
legend(hax{1},legendStr);
General.save_fig(f{1},fname{1},obj.figs_path);

end