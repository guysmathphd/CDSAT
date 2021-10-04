function obj = plot_node_degrees_vs_time2(obj)

nodes_pu_unique_sorted = General.load_var(fullfile(obj.path,'nodes_pu_unique_sorted'));
k_outs_vs_times_out = General.load_var(fullfile(obj.path,'k_outs_vs_times_out'));
times_out = General.load_var(fullfile(obj.path,'times_out'));
node_names_out = General.load_var(fullfile(obj.path,'node_names_out'));
max_k_out = General.load_var(fullfile(obj.path,'max_k_out'));
max_k_out_inds = General.load_var(fullfile(obj.path,'max_k_out_inds'));
min_k_out = General.load_var(fullfile(obj.path,'min_k_out'));
min_k_out_inds = General.load_var(fullfile(obj.path,'min_k_out_inds'));

[max_k_out_sorted,max_k_out_sorted_inds] = sort(max_k_out,'descend');
n = 5;
top_n_inds = max_k_out_sorted_inds(1:n);

% top_n_loc_ids = nodes_pu_unique_sorted(top_n_inds);
% 
% top_n_node_names = node_names_out(top_n_inds);
name1 = 'net07a';name2 = 'net07b';name3 = 'net07c';
f1 = figure('Name',name1,'NumberTitle','off');
hax1 = axes;hold on;figdesc1 = 'k_{i,out} vs time';
f2 = figure('Name',name2,'NumberTitle','off');
hax2 = axes;hold on;figdesc2 = 'Single-Sided Amplitude Spectrum of k_{i,out}(t)';
f3 = figure('Name',name3,'NumberTitle','off');
hax3 = axes;hold on;figdesc3 = 'Single-Sided Angle (Phase shift) Spectrum of k_{i,out}(t)';

legendStr = cell(1,n);
for i1 = 1:n
    cur_ind = top_n_inds(i1);
    cur_loc_id = nodes_pu_unique_sorted(cur_ind);
    cur_node_name = node_names_out(cur_ind);
    cur_k_outs = k_outs_vs_times_out{1,cur_ind};
    cur_times_out = times_out{1,cur_ind};
    [ff,P1,P1a] = Network.myfft(cur_times_out,cur_k_outs);
    P1aa = P1a/(2*pi)./P1;
    plot(hax1,cur_times_out,cur_k_outs);hold on;
    plot(hax2,ff*24*3600,P1);hold on;
    plot(hax3,ff*24*3600,P1aa);hold on;
    legendStr{1,i1} = [cur_node_name{1} ' locID = ' num2str(cur_loc_id)];
end
legend(hax1,legendStr);
xlabel(hax1,'t');
ylabel(hax1,'k_{out}');
title(hax1,{[name1 ' ' obj.name ' ' obj.desc];[figdesc1]});
General.save_fig(f1,name1,fullfile(obj.path,'netfigs'));

legend(hax2,legendStr);
xlabel(hax2,'f (cycles/day)');
ylabel(hax2,'$\mid \mathcal{F}(k_{out}) \mid $','Interpreter','latex');
title(hax2,{[name2 ' ' obj.name ' ' obj.desc];[figdesc2]});
General.save_fig(f2,name2,fullfile(obj.path,'netfigs'));

legend(hax3,legendStr);
xlabel(hax3,'f (cycles/day)');
ylabel(hax3,'$angle[ \mathcal{F}(k_{out}) ] (days)$','Interpreter','latex');
title(hax3,{[name3 ' ' obj.name ' ' obj.desc];[figdesc3]});
General.save_fig(f3,name3,fullfile(obj.path,'netfigs'));
end