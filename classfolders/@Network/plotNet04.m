function obj = plotNet04(obj)
netpath = obj.path;maxmarksize = 20;
% R = General.load_var(fullfile(netpath,'net09-R'));
k = General.load_var(fullfile(netpath,'degree_vector'));
rng(2);N=obj.N; marksize = log10(k./min(k));
marksize = marksize*maxmarksize/max(marksize)+1; %marksize = 1;
alpha1 = .5;alpha2 = .1;
R = rand(N,1);
thetas = General.load_var(fullfile(netpath,'net09-thetas'));
fname = 'net10a';
f = figure('Name',fname,'NumberTitle','off');
colormap(parula);
cmap = colormap;cmap_size = size(cmap,1);
[x,y] = pol2cart(thetas,R);z = zeros(size(x));
s=scatter3(x,y,z,marksize,'o',...
    'MarkerFaceColor',cmap(end,:),'MarkerEdgeColor',cmap(1,:));
%     'MarkerFaceColor','w','MarkerEdgeColor','w');
% set(gca,'Color','k');
hax = f.Children;
% hax.Color = 'k';
hax.XTick = [];
hax.YTick = [];
hax.ZTick = [];
f.InvertHardcopy = 'off';
s.MarkerFaceAlpha = alpha1;
s.MarkerEdgeAlpha = alpha1;
hax = f.Children;
hax.ZLim = [0 1];


General.save_fig(f,fname,fullfile(obj.path,'figs'));
% hold on;
% f.Name = 'net10b';
% 
% % nodes = hax.Children;
% s.MarkerFaceAlpha = alpha2;
% s.MarkerEdgeAlpha = alpha2;
% for i1 = 1:100
%     theta = thetas + 2*i1*pi/200;
%     [x,y] = pol2cart(theta,R);
%     scatter(x,y,marksize,'ow','MarkerFaceAlpha',alpha2,'MarkerEdgeAlpha',alpha2);
% 
% end

end