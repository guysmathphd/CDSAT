function obj = plot_pert_props_vs_times(obj)
namepre = 'set5';
% different figs
Qs = {1,'Q2','Q3','Q4','Q5','Q6','Q7','Q8','Q9','Q10'};
prop_sufs = {'_norms','_norms_normed','_thetas'};
fig_suf_1 = 97;
colors  = rand(1000,3);
ylabs = {'$\mid \delta x(t)\mid$','$A(t)$','$\theta(t)$'};
xlabs = {};
% different data in each fig
solution_folders = {'solution_random_perts','single_node_pert_sol','solution_eigvec_perts'};

for i1 = 1:length(prop_sufs)
    suf1 = char(fig_suf_1);fig_suf_1 = fig_suf_1 + 1;
    for i2 = 1:length(Qs)
        Q_suf = Qs{i2};
        xlabs{i2} = ['$t/' Q_suf '$'];
        if ~isstr(Q_suf)
            Q_suf = num2str(Q_suf);
            xlabs{i2} = 't';
        end
        names{i1,i2} = [namepre suf1 '-' Q_suf];
        f{i1,i2} = figure('Name',names{i1,i2},'NumberTitle','off');hax{i1,i2} = axes;% hold on;
    end
end
legendStr = {};colorind = 0;
for j1 = 1:length(obj.objPropertiesPaths)
    proppath = obj.objPropertiesPaths{j1};
    for j2 = 1:length(solution_folders)
        solfold = solution_folders{j2};
        folderpath = fullfile(proppath,solfold);
        t_files = dir(fullfile(folderpath,'sol_t_*.mat'));
        for j3 = 1:length(t_files)
            t_file = t_files(j3).name;
            t = General.load_var(fullfile(folderpath,t_file));
            x_root = ['sol_x_' t_file(7:end-4)];
            colorind=colorind+1;
            for i1 = 1:length(prop_sufs)
                prop_suf = prop_sufs{i1};
                x_file = [x_root prop_suf '.mat'];
                x = General.load_var(fullfile(folderpath,x_file));
                for i2 = 1:length(Qs)
                    if i2 == 1
                        q = 1;
                    else
                        Q = Qs{i2};
                        Q_file = [x_root '_' Q];
                        q = General.load_var(fullfile(folderpath,Q_file));
                    end
                    disp(fullfile(folderpath,t_file));
                    disp(fullfile(folderpath,x_file));
                    loglog(hax{i1,i2},t/q,x,'Color',colors(colorind,:));
                    hold(hax{i1,i2},'on');
                end
            end
        end
    end
end
for i1 = 1:length(prop_sufs)
    for i2 = 1:length(Qs)
        xlim = hax{i1,i2}.XLim;
        ylim = hax{i1,i2}.YLim;
        loglog(hax{i1,i2},[xlim(1),1],[.5,.5],'--k','LineWidth',2);
        loglog(hax{i1,i2},[1,1],[ylim(1),.5],'--k','LineWidth',2);
        xlabel(hax{i1,i2},xlabs{i2},'interpreter','latex');
        ylabel(hax{i1,i2},ylabs{i1},'interpreter','latex');
        title(hax{i1,i2},{[names{i1,i2} ' ' obj.name]});
        General.save_fig(f{i1,i2},names{i1,i2},obj.figs_path);
    end
end
end