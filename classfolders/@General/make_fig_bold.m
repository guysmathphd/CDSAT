function make_fig_bold(folderpath)
files = dir(fullfile(folderpath,'*.fig'));
for i1 = 1:length(files)
    file = files(i1);a = 0;
    if a
        continue
    end
    f = openfig(fullfile(folderpath,file.name));
    set(gca,'TickLength',[.025 .025]);
    set(gca,'LineWidth',5);
    set(gca,'FontSize',26);
    set(gca,'FontWeight','bold');
    xs = [];ys=[];
    h = get(gca,'Children');
    for i2 = 1:length(h)
        set(h(i2),'LineWidth',4.5);
        set(h(i2),'MarkerSize',13);
        set(h(i2),'MarkerFaceColor','none')
        xs = [xs h(i2).XData];
        ys = [ys h(i2).YData];
    end
    xmin = min(xs);xmax=max(xs);
    ymin = min(ys);ymax=max(ys);
    
    legend('off');
    %General.save_fig(f,file.name,folderpath)
end