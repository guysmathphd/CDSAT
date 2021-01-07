classdef EngineSet < handle
    
    properties
        engines
        path
        name
        engine_paths
        enginesPropertiesMaps
        propertiesMapsDims
        eng_keys
        keys_sampled
        num_engines
        eigenvalues_ana
        eigenvalues_asy
        valid_ids
        index_mat
        keys_sampled_values
        figs_path
        missing_scenario
        inits_legit
    end
    
    methods
        function obj = EngineSet(path,name)
            obj.path = path;
            obj.name = name;
            if ~isfolder(obj.path)
                mkdir(obj.path)
            end
            obj.make_fig_folder();
            obj.set_properties();
            obj.set_engine_paths();
            obj.save_obj();
        end
        function obj = make_fig_folder(obj)
            obj.figs_path = fullfile(obj.path,'setfigs');
            if ~isfolder(obj.figs_path)
                mkdir(obj.figs_path);
            end
        end
        function obj = set_engine_paths(obj)
            contents = dir(obj.path);
            count = 0;
            for i = 1:size(contents,1)
                folder = contents(i).name;
                k = strfind(folder,'-');
                if ~isempty(k)
                    j = str2double(folder(k+1:end));
                    obj.engine_paths{j} = fullfile(obj.path,folder,'results',[folder 'Obj.mat']);
                    count = count+1;
                end
            end
            obj.num_engines = count;
        end
        function obj = set_properties(obj)
            obj.propertiesMapsDims = size(obj.enginesPropertiesMaps);
            inds = find(obj.propertiesMapsDims>1);
            obj.eng_keys = keys(obj.enginesPropertiesMaps{1});
            obj.keys_sampled = obj.eng_keys(inds);
            dims = obj.propertiesMapsDims(inds);
            obj.index_mat = zeros(dims);
            for i = 1:obj.num_engines
                [r,c] = ind2sub(dims,i);
                obj.index_mat(r,c) = i;
            end
%             for i = 1:36
%                 map = obj.enginesPropertiesMaps{i};
%                 str1 = func2str(map('f_M0'));
%                 str2 = func2str(map('f_M2'));
%                 disp(['i = ' num2str(i)]);
%                 disp(['f_M0 = ' str1]);
%                 disp(['f_M2 = ' str2]);
%             end
            for i = 1:length(dims)
                key = obj.keys_sampled{i};
                for j = 1:dims(i)
                    dims_j = obj.propertiesMapsDims;
                    dims_j(inds(i)) = j;
                    str = 'map=obj.enginesPropertiesMaps{';
                    for k = 1:length(dims_j)
                        str = [str num2str(dims_j(k)) ','];
                    end
                    str = [str(1:end-1) '}'];
                    disp(str);
                    eval(str); %map = obj.enginesPropertiesMaps{dims_j};
                    obj.keys_sampled_values{i,j} = map(key);
                end
            end
            
        end
        function obj = save_obj(obj)
            save(fullfile(obj.path,[obj.name 'Obj.mat']),'obj','-v7.3');
        end
        function obj = create_engines(enginesPropertiesMaps)
            obj.enginesPropertiesMaps = enginesPropertiesMaps;
            ind = 1;
            for p = obj.enginesPropertiesMaps
                disp(['EngienSet ind = ' num2str(ind)]);
                p1 = p{1};
                p1('scenarioName') = [p1('scenarioName') '-' num2str(ind)];
                p1('resultsPath') = fullfile(p1('resultsPath'),p1('scenarioName'),'results');
                EngineClass(p1);
                ind = ind+1;
            end
            obj.save_obj();
        end
        function obj = setSolve(obj,startInd,useEigVec,useSS,epsIndStart,vecIndStart,isAdjustEpsilon)
            for i = startInd:obj.num_engines
                eng_path = obj.engine_paths{i};
                disp(['eng_path = ' eng_path]);
                eng = load(eng_path);
                eng.obj.solverTimeStep = .2;
                eng.obj.maxTime = 1000;
                %%% temporary
                if ~useEigVec
                eng.obj.solution_t=[];eng.obj.solution_x=[];
                end
                eng.obj.solution_t_eigvecana=[];eng.obj.solution_t_eigvecasy=[];
                eng.obj.solution_x_eigvecana=[];eng.obj.solution_x_eigvecasy=[];
                eng.obj.solve(useEigVec,useSS,epsIndStart,vecIndStart,isAdjustEpsilon);
                clear eng;
            end
        end
        function obj = set_eigenvalues(obj,startInd,endInd)
            if isempty(obj.eigenvalues_ana)
                obj.eigenvalues_ana = cell(1,obj.num_engines);
                obj.eigenvalues_asy = cell(1,obj.num_engines);
                obj.missing_scenario = [];
            end
            for i=startInd:min(obj.num_engines,endInd)
                try
                disp(['scenario' num2str(i)]);
                eng_path = obj.engine_paths{i};
                eng = load(eng_path);
                obj.eigenvalues_ana{i} = eng.obj.eigenvalues_ana;
                obj.eigenvalues_asy{i} = eng.obj.eigenvalues_asy;
                clear eng;
                catch exception
                    disp(['Error with scenario ' num2str(i)]);
                    disp(exception.message);
                    obj.missing_scenario(end+1) = j;
                end
            end            
            obj.save_obj();
        end
        %%% To be added to obj.set_eigenvalues
        function obj = set_isInitsLegit(obj)
            inits_legit = -1*ones(1,obj.num_engines);
            for i=1:obj.num_engines
                disp(['scenario' num2str(i)]);
                eng_path = obj.engine_paths{i};
                eng = load(eng_path);
                inits_legit(i) = eng.obj.initsLegit;
            end
            obj.save_obj();
        end
                
        
        function obj = setPlot_results(obj,startInd,isDiluted)
            for i = startInd:obj.num_engines
                eng_path = obj.engine_paths{i};
                eng = load(eng_path);
                eng.obj.plot_results(isDiluted);
                clear eng;
            end
        end
        function obj = setFunction(obj,function_handle,startInd)
            for i = startInd:obj.num_engines
                eng_path = obj.engine_paths{i};
                eng = load(eng_path);
                function_handle(eng.obj);
                clear eng;
            end
        end
        function obj = plot_index(obj)
            name = 'set1a';
            figdesc = 'Scenario Index';
            f = figure('Name',name,'NumberTitle','off');
            image(obj.index_mat,'CDatamapping','scaled');
            colorbar;
            ylabel(obj.keys_sampled{1}(3:end));
            xlabel(obj.keys_sampled{2}(3:end));
            labels1 = {};
            labels2 = {};
            for j = 1:size(obj.keys_sampled_values,2)
                str1 = func2str(obj.keys_sampled_values{1,j});
                str2 = func2str(obj.keys_sampled_values{2,j});
                labels1{j} = str1(6:end-1);
                labels2{j} = str2(6:end-1);
            end
            yticklabels(labels1);
            xticklabels(labels2);
            title({[name ' ' obj.name];figdesc},'interpreter','latex');
            obj.save_fig(f,name);
        end
        function obj = plot_setEigenvalues(obj)
            %%% set2a
            n = 1:obj.num_engines;
            n(obj.missing_scenario) = [];
            
            numeigen = [size(obj.eigenvalues_ana{1},1), 5];
            suf = {['-first' num2str(numeigen(1))],'-first5'};
            for i = 1:2    
                n1 = n;
                to_delete = [];
                name = ['set2a' suf{i}];
                figdesc = ['Eigenvalues Comparison' suf{i}];
                f = figure('Name',name,'NumberTitle','off');
                ind = 1;max_val=[];min_val=[];range_val=[];
                for j = 1:length(n)
                    disp(['j=' num2str(j)]);
                    eigvals_ana = obj.eigenvalues_ana{n(j)};
                    if ~isempty(eigvals_ana)
                        eigvals_ana = eigvals_ana(1:numeigen(i));
                        max_val(ind) = max(eigvals_ana);
                        min_val(ind) = min(eigvals_ana);
                        ind = ind+1;
                    else
                        to_delete(end+1) = j;
                    end
                end
                n1(to_delete) = [];
                range_val = max_val - min_val;
                yyaxis left;
                plot(n1,max_val,'*');
                title({[name ' ' obj.name];figdesc},'interpreter','latex');
                xlabel('Scenario Index');
                ylabel('Max Eigenvalue');
                yyaxis right;
                plot(n1,range_val,'o');
                ylabel('$Range (\lambda_{max} - \lambda_{min})$','interpreter','latex');
                obj.save_fig(f,name);
            end
        end
        function obj = save_fig(obj,f,name)
            try
                saveas(f,fullfile(obj.path,'setfigs',name),'fig');
            catch exception
                disp(exception.message);
            end
            try
                saveas(f,fullfile(obj.path,'setfigs',name),'png');
            catch exception
                disp(exception.message);
            end
        end
    end
end

