function obj = plot_pert_props_vs_times_1(obj,solution_folders,t_file_node_ids_str)


namepre = 'set5';
% different figs
Qs = {'1','Q2','Q3','Q4','Q5','Q6','Q7','Q8','Q9','Q10','tau_{A_1}','tau_{A_2}'};
prop_sufs = {'_norms','_norms_normed','_thetas'};
fig_suf_1 = 97;
colors  = General.colors1;
ylabs = {'$\mid \delta x(t)\mid$','$A(t)$','$\theta(t)$'};
xlabs = {};
% different data in each fig
% solution_folders = {'solution_random_perts','single_node_pert_sol'};%,'solution_eigvec_perts'};
tau_ind = 1;
for i1 = 1:length(prop_sufs)
    suf1 = char(fig_suf_1);fig_suf_1 = fig_suf_1 + 1;
    for i2 = 1:length(Qs)
        Q_suf = Qs{i2};
        if strcmp(Q_suf,'1')
            xlabs{i2} = 't';
        elseif strcmp(Q_suf(1:end-6),'tau')
            xlabs{i2} = ['$\' Q_suf '$'];            
        else
            Q_suf = ['tau-' num2str(tau_ind)];
            tau_ind = tau_ind + 1;
        end
        names{i1,i2} = [namepre suf1 '-' Q_suf];
        f{i1,i2} = figure('Name',names{i1,i2},'NumberTitle','off');hax{i1,i2} = axes;% hold on;
    end
end
legendStr = {};colorind = 0;legend_subset = cell(length(prop_sufs),length(Qs));
for j1 = 1:length(obj.objPropertiesPaths)
    proppath = obj.objPropertiesPaths{j1};
    mycolor = colors(j1,:);
    new_legend_subset = true;
    for j2 = 1:length(solution_folders)
        solfold = solution_folders{j2};
        folderpath = fullfile(proppath,solfold);
        t_files = dir(fullfile(folderpath,'sol_t_*.mat'));
        for j3 = 1:length(t_files)
            t_file = t_files(j3).name;
            if any(strcmp(t_file(7:end-4),t_file_node_ids_str))
                t = General.load_var(fullfile(folderpath,t_file));
                if t(end) < 1
                    'stop here'
                end
                x_root = ['sol_x_' t_file(7:end-4)];
                colorind=colorind+1;
                for i1 = 1:length(prop_sufs)
                    prop_suf = prop_sufs{i1};
                    x_file = [x_root prop_suf '.mat'];
                    x = General.load_var(fullfile(folderpath,x_file));
                    tau_ind = 1;
                    for i2 = 1:length(Qs)
                        
                        Q = Qs{i2};
                        if strcmp(Q,'1')
                            q = 1;
                        elseif strcmp(Q(1:end-6),'tau')
                            taus = General.load_var(fullfile(proppath,'sys_half_life_amp'));
                            tau = taus{tau_ind};tau_ind = tau_ind + 1;
                            q = tau;
                        else
                            Q_file = [x_root '_' Q];
                            q = General.load_var(fullfile(folderpath,Q_file));
                        end
                        disp(fullfile(folderpath,t_file));
                        disp(fullfile(folderpath,x_file));
                        p = loglog(hax{i1,i2},t/q,x,'Color',mycolor);
                        if new_legend_subset
                            legend_subset{i1,i2}(end+1) = p;
                        end
                        hold(hax{i1,i2},'on');
                    end
                end
                new_legend_subset = false;
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
        legend(hax{i1,i2},legend_subset{i1,i2},obj.dynamicNames,'interpreter','latex');
        General.save_fig(f{i1,i2},names{i1,i2},obj.figs_path);
    end
end
end