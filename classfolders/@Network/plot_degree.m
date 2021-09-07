function obj = plot_degree(obj)
%net04a net04b
data_nodes = readtable(fullfile(obj.source_data_path,'nodes'));
node_names = data_nodes.Name;
k_out = General.load_var(fullfile(obj.path,'degree_vector'));
k_in = General.load_var(fullfile(obj.path,'degree_vector_in'))';
ks = {k_out,k_in};num_bars_display_short = 100;
figdesc1 = {'k_{out} sorted','k_{in} sorted'};namepre = 'net04';
namesuf = {'a','b'};
for i1 = 1:length(ks)
    k = ks{i1};
    [k_sorted,k_sorted_inds] = sort(k,'descend');
    node_names_sorted = node_names(k_sorted_inds);

    Y = k_sorted;
    num_bars_display = [length(k_sorted),num_bars_display_short];
    for i2 = 1:length(num_bars_display)
        num_bars = num_bars_display(i2);
        X = categorical(node_names_sorted(1:num_bars));
        X = reordercats(X,node_names_sorted(1:num_bars));
        name = [namepre namesuf{i1} '-top' num2str(num_bars)];
        f = figure('Name',name,'NumberTitle','off');
        bar(X,Y(1:num_bars));
        xlabel('Location');
        ylabel('k');
        title({[name ' ' obj.name ' ' obj.desc];[figdesc1{i1} ' top ' num2str(num_bars)]});
        General.save_fig(f,name,fullfile(obj.path,'netfigs'));
    end
end

%net04c
k_tot = k_out + k_in;
[k_tot_sorted, k_tot_sorted_inds] = sort(k_tot,'descend');
node_names_sorted = node_names(k_tot_sorted_inds);
k_out_sorted = k_out(k_tot_sorted_inds);namepre = 'net04c';
figdesc = 'k_{out} and k_{in} sorted by k_{total}';
k_in_sorted = k_in(k_tot_sorted_inds);Y = [k_out_sorted k_in_sorted];
for i2 = 1:length(num_bars_display)
    num_bars = num_bars_display(i2);
    X = categorical(node_names_sorted(1:num_bars));
    X = reordercats(X,node_names_sorted(1:num_bars));
    name = [namepre '-top' num2str(num_bars)];
    f = figure('Name',name,'NumberTitle','off');
    h = bar(X,Y(1:num_bars,:),'stacked');
%     h = bar(X,Y(1:num_bars,:));
    set(h, {'DisplayName'}, {'k_{out}';'k_{in}'});
    legend();
%     xtips1 = h(1).XEndPoints;
%     ytips1 = h(1).YEndPoints;
%     labels1 = string(h(1).YData);
%     text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%         'VerticalAlignment','bottom')
%     xtips2 = h(2).XEndPoints;
%     ytips2 = h(2).YEndPoints;
%     labels2 = string(h(2).YData);
%     text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
%         'VerticalAlignment','bottom')
    xlabel('Location');
    ylabel('k');
    title({[name ' ' obj.name ' ' obj.desc];[figdesc ' top ' num2str(num_bars)]});
    General.save_fig(f,name,fullfile(obj.path,'netfigs'));
end

%net04d
k_ratio = k_out_sorted./k_tot_sorted;
namepre = 'net04d';
figdesc = 'Ratio of k_{out}/k_{total} sorted by k_{total}';
for i2 = 1:length(num_bars_display)
    num_bars = num_bars_display(i2);
    name = [namepre '-top' num2str(num_bars)];
    f = figure('Name',name,'NumberTitle','off');
    yyaxis left
    plot(k_ratio(1:num_bars));
    xlabel('Locations');
    ylabel('k_{out}/k_{total}');
    xticks(1:num_bars);
    xticklabels(node_names_sorted(1:num_bars));
    xtickangle(45);a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',8)%,'XTickLabelMode','auto');
    yyaxis right
    plot(k_tot_sorted(1:num_bars));
    ylabel('k_{total}');
    title({[name ' ' obj.name ' ' obj.desc];[figdesc ' top ' num2str(num_bars)]});
    General.save_fig(f,name,fullfile(obj.path,'netfigs'));
end



end