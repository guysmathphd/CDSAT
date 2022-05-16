function obj = plotNet02(obj)
% net09a
fname = 'net09a';
N = obj.N;
rng(1);
thetas = rand(N,1)*2*pi;General.save_var(thetas,obj.path,'net09-thetas');
phis = rand(N,1)*pi - pi/2;General.save_var(phis,obj.path,'net09-phis');
k = General.load_var(fullfile(obj.path, 'degree_vector.mat'));
maxk = max(k);
R = log(maxk./k);General.save_var(R,obj.path,'net09-R');
[X,Y,Z] = sph2cart(thetas,phis,R);
mink = min(k);
f = figure('Name',fname,'NumberTitle','off');
cmap = jet;
[ind_bins_var,bins_edges,ind_bins_var_sizes] = EngineClass.set_bins_generic(256,k,.0001,true(size(k)));
try
    for i1 = 1:N
        x = X(i1);
        y = Y(i1);
        z = Z(i1);
        ind = find(bins_edges > k(i1),1,'first')-1;
        plot3(x,y,z,'o','MarkerSize',20*k(i1)/maxk,'MarkerFaceColor',cmap(ind,:),...
            'MarkerEdgeColor',cmap(ind,:));hold on;
        %     pause(.01);
    end
catch exception
    disp(['i1 = ' num2str(i1)]);
    disp(exception.message);
end
General.save_fig(f,fname,fullfile(obj.path,'figs'));

% net09b
fname = 'net09b';
f = figure('Name',fname,'NumberTitle','off');
try
    for i1 = 1:N
        x = X(i1);
        y = Y(i1);
        z = 0;
        ind = find(bins_edges > k(i1),1,'first')-1;
        plot3(x,y,z,'o','MarkerSize',20*k(i1)/maxk,'MarkerFaceColor',cmap(ind,:),...
            'MarkerEdgeColor',cmap(ind,:));hold on;
        %     pause(.01);
    end
catch exception
    disp(['i1 = ' num2str(i1)]);
    disp(exception.message);
end
General.save_fig(f,fname,fullfile(obj.path,'figs'));

% net09c
fname = 'net09c';
f = figure('Name',fname,'NumberTitle','off');
num_circles = 7;node_inds = (1:obj.N)';
R1 = 2*floor(num_circles*node_inds/(obj.N+1)) + 2;
[X,Y,Z] = pol2cart(thetas, R1, zeros(size(R1)));S = 20;
C = [0 0.4470 0.7410];
s = scatter3(X,Y,Z,S,C,'MarkerFaceColor', C);
ax = gca;ax.ZLim = [0 1];ax.XLim = [-max(R1),max(R1)];
ax.YLim = [ax.XLim];
General.save_fig(f,fname,fullfile(obj.path,'figs'));
end