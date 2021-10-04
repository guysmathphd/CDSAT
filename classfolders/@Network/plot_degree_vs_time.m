function obj = plot_degree_vs_time(obj)
%net05a
total_weight = General.load_var(fullfile(obj.path,'total_weight_vs_all_times_sorted'));
all_times_sorted = General.load_var(fullfile(obj.path,'all_times_sorted'));
ys = {total_weight};
if obj.average_weight_denominator_type
    average_weight_vs_all_times_sorted_cut = General.load_var(fullfile(obj.path,'average_weight_vs_all_times_sorted'));
    ys{end+1} = average_weight_vs_all_times_sorted_cut;
end

name_suf = {'','-avg'};figdescpre = {'Total','Average'};
ylabelstr = {'$k_{total} = \sum_i k_{i}$', ...
    ['$k_{avg} = k_{total} / k_{total, ' obj.average_weight_denominator_type '}$']};
for i1 = 1:length(ys)
    name = ['net05a' name_suf{i1}];
    figdesc = [figdescpre{i1} ' network weight vs time'];
    f = figure('Name',name,'NumberTitle','off');
    plot(all_times_sorted,ys{i1});
    xlabel('Times');
    ylabel(ylabelstr{i1},'Interpreter','latex');
    title({[name ' ' obj.name ' ' obj.desc];figdesc});
    General.save_fig(f,name,fullfile(obj.path,'netfigs'));
    
    %net05b
    
    
    [ff,P1] = Network.myfft(all_times_sorted,ys{i1});
    ff = ff*24*3600; % convert to cycles per day
    
    name = ['net05b' name_suf{i1}];
    figdesc = ['Single-Sided Amplitude Spectrum of k_{' figdescpre{i1} '}(t)'];
    f = figure('Name',name,'NumberTitle','off');
    plot(ff,P1);
    xlabel('f (cycles/day)');
    ylabel('$\mid \mathcal{F}(k) \mid $','Interpreter','latex');
    title({[name ' ' obj.name ' ' obj.desc];figdesc});
    General.save_fig(f,name,fullfile(obj.path,'netfigs'));
    
end
end