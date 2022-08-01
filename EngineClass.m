classdef EngineClass <  handle
    %ENGINECLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        scenarioName %String
        desc %String
        adjacencyMatrix %2D array
        networkName
        initialValues %1D array same length as adjacency matrix
        maxTime %positive integer
        maxDerivative % positive double
        solverTimeStep % positive integer
        randSeed %positive integer
        f_M0 %function handle
        f_M1 %function handle
        f_M2 %function handle
        f_R %function handlef
        f_dM0 %function handle, derivate of M0
        f_dM1 %derivative of M1
        f_dM2 %derivative of M2
        Dii_ana % vector representing diagonal of Jacobian
        Wij_ana % matrix representing Jacobian without diagonal
        Dii_asy
        Wij_asy
        Dii_num
        Wij_num
        M2_i_bigodot
        M2_k_bigodot
        perturbations
        random_perturbations
        pert_eigvec_ana_1
        pert_eigvec_asy_1
        num_perturbations = 1;
        half_life
        eigenvalues_ana
        eigenvectors_ana
        eigenvectors_ana_binned_self
        eigenvectors_ana_binned_k
        eigenvectors_ana_binned_kinn
        eigenvalues_asy
        eigenvalues_asy_permuted
        eigenvectors_asy
        eigenvectors_asy_permuted
        eigenvectors_asy_permuted_binned_self
        eigenvectors_asy_permuted_binned_k
        eigenvectors_asy_permuted_binned_kinn
        eigenvalues_num
        eigenvectors_num
        node_relax_time
        sys_relax_time
        node_half_life
        node_half_lives_struct
        sys_half_life_amp = 1;
        
        N
        k_nn
        ki_nn
        ki_nnbinned
        numbins
        bins
        binsij
        binskinn
        binseigvecana
        binseigvecasy_permuted
        kbinned
        kbinned_mins
        kbinned_maxs
        kijbinned
        Dii_anabinned
        Dii_asybinned
        Wij_anabinned
        Wij_asybinned
        C_W_set
        C_D_set
        C_D_v3
        guesses1_v3
        guesses2_v3
        errors1_v3
        errors2_v3
        eigvals_v3
        eigvecs_v3
        numeigen;
        numeigenplots = 5;
        numeigenplots2 = 200;
        eigvec_dist_comparison_mat_ana2asy
        eigvec_angle_comparison_mat_ana2asy
        eigvec_dist_comparison_mat_ana2asy_permuted
        eigvec_angle_comparison_mat_ana2asy_permuted
        eigvec_dot_comparison_mat_ana2asy
        eigvec_dot_comparison_mat_ana2asy_permuted
        eigvec_dot_comparison_mat_ana2ana
        eigvec_dot_comparison_mat_asy_permuted2asy_permuted
        eigvec_dot_comparison_mat_asy2asy
        permutation_eigvec_ana2asy
        eps = [1];
        eps_adjusted;
        isInitsLegit;
        epsThreshold = .01;
        perturbation_factor = .1;
        pert_map_struct
        pert_map_struct_v2asy
        init_condition_str
        stop_condition_2
        opts
        epsFactor = 0.9;
        
        difEqSolver %function handle
        absTol % positive double
        relTol % positive double
        solution_t = 0;
        solution_x ;
        solution_t_eigvec;
        solution_x_eigvec;
        solution_t_eigvecana;
        solution_t_eigvecasy;
        solution_x_eigvecana;
        solution_x_eigvecasy;
        solution_t_perturbations;
        solution_x_perturbations;
        resultsPath
        header
        degree_vector
        degree_vector_weighted
        degree_weighted_average
        steady_state
        steady_state_calculated
        mu = 0
        nu = 0
        rho = 0
        eta = 0
        xi = 0
        isEngineSet = false
        init_condition
        % J_ana
        % J_asy
    end
    methods
        function obj = EngineClass(propertiesMap)
            %ENGINECLASS Construct an instance of this class
            %   Loop over keys of propertiesMap, assign each value to
            %   property of EngineClass
            for key = keys(propertiesMap)
                eval(['obj.' key{1} '=propertiesMap(key{1})']);
            end
            obj.init_object();
            obj.save_obj();
        end
        function set_num_random_perturbations(obj)
            obj.num_perturbations = 5;
            obj.save_obj();
        end
        function obj = clean_up_obj(obj)
            disp(obj.scenarioName);
            x = General.getSize(obj);
            disp(x);
            propspath = fullfile(obj.resultsPath,'obj_properties');
            obj.adjacencyMatrix = [];
            General.save_var(obj.Dii_ana,propspath,'Dii_ana');obj.Dii_ana=[];
            General.save_var(obj.Wij_ana,propspath,'Wij_ana');obj.Wij_ana=[];
            General.save_var(obj.Dii_asy,propspath,'Dii_asy');obj.Dii_asy=[];
            General.save_var(obj.Wij_asy,propspath,'Wij_asy');obj.Wij_asy=[];
            General.save_var(obj.eigenvalues_ana,propspath,'eigenvalues_ana');obj.eigenvalues_ana=[];
            General.save_var(obj.eigenvectors_ana,propspath,'eigenvectors_ana');obj.eigenvectors_ana=[];
            General.save_var(obj.eigenvalues_asy,propspath,'eigenvalues_asy');obj.eigenvalues_asy=[];
            General.save_var(obj.eigenvectors_asy,propspath,'eigenvectors_asy');obj.eigenvectors_asy=[];
            General.save_var(obj.eigenvalues_asy_permuted,propspath,'eigenvalues_asy_permuted');obj.eigenvalues_asy_permuted=[];
            General.save_var(obj.eigenvectors_asy_permuted,propspath,'eigenvectors_asy_permuted');obj.eigenvectors_asy_permuted=[];
            obj.ki_nn = [];
            General.save_var(obj.eigvec_angle_comparison_mat_ana2asy,propspath,'eigvec_angle_comparison_mat_ana2asy');obj.eigvec_angle_comparison_mat_ana2asy=[];
            General.save_var(obj.eigvec_angle_comparison_mat_ana2asy_permuted,propspath,'eigvec_angle_comparison_mat_ana2asy_permuted');obj.eigvec_angle_comparison_mat_ana2asy_permuted=[];
            General.save_var(obj.eigvec_dist_comparison_mat_ana2asy,propspath,'eigvec_dist_comparison_mat_ana2asy');obj.eigvec_dist_comparison_mat_ana2asy=[];
            General.save_var(obj.eigvec_dist_comparison_mat_ana2asy_permuted,propspath,'eigvec_dist_comparison_mat_ana2asy_permuted');obj.eigvec_dist_comparison_mat_ana2asy_permuted=[];
            General.save_var(obj.eigvec_dot_comparison_mat_ana2asy,propspath,'eigvec_dot_comparison_mat_ana2asy');obj.eigvec_dot_comparison_mat_ana2asy=[];
            General.save_var(obj.eigvec_dot_comparison_mat_ana2asy_permuted,propspath,'eigvec_dot_comparison_mat_ana2asy_permuted');obj.eigvec_dot_comparison_mat_ana2asy_permuted=[];
            General.save_var(obj.permutation_eigvec_ana2asy,propspath,'permutation_eigvec_ana2asy');obj.permutation_eigvec_ana2asy=[];
            obj.eigvec_dot_comparison_mat_ana2ana = [];
            obj.eigvec_dot_comparison_mat_asy_permuted2asy_permuted=[];
            obj.eigvec_dot_comparison_mat_asy2asy=[];
            
            x = General.getSize(obj);
            disp(x);
            obj.save_obj();
        end
        function obj = test_1(obj)
            global obj;
            addpath(fullfile(obj.resultsPath,'obj_properties'));
            y = testfunction(2);
            disp(num2str(y));
        end
        function create_odefun(obj)
            fileID = fopen(fullfile(obj.resultsPath,'obj_properties','odefun.m'),'w+');
            fprintf(fileID,'%0s\n','function dy=odefun(t,x)');
            fprintf(fileID,'%0s\n','global obj;');
            fprintf(fileID,'%0s\n',['dy = zeros(' num2str(obj.N) ',1);']);
            for i = 1:obj.N
                fprintf(fileID,'%0s',['dy(' num2str(i) ') = obj.f_M0(x(' num2str(i) ')) + obj.f_M1(x(' num2str(i) ')) * (']);
                first = true;
                for j = 1:obj.N
                    if ~first
                        if obj.adjacencyMatrix(i,j) ~=0
                            fprintf(fileID,'%0s','+');
                        end
                    end
                    if obj.adjacencyMatrix(i,j) ~=0
                        fprintf(fileID,'%0s',[num2str(obj.adjacencyMatrix(i,j)) '* obj.f_M2(x(' num2str(j) '))']);
                        first = false;
                    end
                    
                end
                fprintf(fileID,'%0s\n',');');
            end
            fprintf(fileID,'%0s\n','end');
            fclose(fileID);
            %             odefun = @(tt,x) (obj.f_M0(x) + (obj.adjacencyMatrix*obj.f_M2(x)).*obj.f_M1(x));
        end
        function obj = init_object(obj)
            obj.solution_x = obj.initialValues';
            obj.numeigen = size(obj.adjacencyMatrix,1);
            % create header
            obj.header = {'t'};
            for i = 1:length(obj.initialValues)
                obj.header{i+1} = ['x_{' num2str(i) '}'];
            end
            if ~isfolder(obj.resultsPath)
                mkdir(obj.resultsPath)
            end
            if ~isfolder(fullfile(obj.resultsPath,'figs'))
                mkdir(fullfile(obj.resultsPath,'figs'));
            end
            disp('calculate_degree');obj.calculate_degree();
            disp('calculate_degree_weighted');obj.calculate_degree_weighted();
            disp('set_dM0');obj.set_dM0();
            disp('set_dM1');obj.set_dM1();
            disp('set_dM2');obj.set_dM2();
            disp('set_R');obj.set_R();
            disp('set_const_functions');obj.set_const_functions();
            disp('set_N');obj.set_N();
            disp('set_knn');obj.set_knn();
            disp('set_Dii_asy');obj.set_Dii_asy();
            disp('set_Wij_asy');obj.set_Wij_asy();
            obj.set_J_asy();
            obj.set_J_asy_degree_vector_weighted();
            obj.set_J_asy_k_inn();
            disp('set_eig_asy');obj.set_eig_asy();
            tol = 1e-13;
            disp('obj.bins = obj.set_bins_generic');
            obj.bins = obj.set_bins_generic(obj.numbins,obj.degree_vector_weighted,tol,true(obj.N,1));
            disp('obj.kbinned = obj.set_binned_vals'); [obj.kbinned,obj.kbinned_mins,obj.kbinned_maxs] = obj.set_binned_vals(obj.degree_vector_weighted,obj.bins);
            disp('obj.Dii_asybinned = obj.set_binned_vals'); obj.Dii_asybinned = obj.set_binned_vals(obj.Dii_asy,obj.bins);
            x = obj.degree_vector_weighted;x2 = (x.^obj.nu) * (x.^obj.rho)';x3 = x2(:);
            disp('obj.binsij = obj.set_bins_generic'); obj.binsij = obj.set_bins_generic(obj.numbins,x3,tol,obj.adjacencyMatrix>0);
            disp('obj.kijbinned = obj.set_binned_vals'); obj.kijbinned = obj.set_binned_vals(x2(obj.adjacencyMatrix>0),obj.binsij);
            Wij_asy = General.load_var(fullfile(obj.resultsPath,'obj_properties','Wij_asy'));
            disp('obj.Wij_asybinned = obj.set_binned_vals'); obj.Wij_asybinned = obj.set_binned_vals(Wij_asy(obj.adjacencyMatrix>0),obj.binsij);
            disp('set_kinn'); obj.set_kinn();
            disp('obj.binskinn = obj.set_bins_generic'); obj.binskinn = obj.set_bins_generic(obj.numbins,obj.ki_nn,tol,true(obj.N,1));
            disp('obj.ki_nnbinned = obj.set_binned_vals'); obj.ki_nnbinned = obj.set_binned_vals(obj.ki_nn,obj.binskinn);
            disp('obj.set_degree_weighted_average'); obj.set_degree_weighted_average();
        end
        function obj = write_gephi_nodes_tables_for_article_fig1(obj)
            sol_path = fullfile(obj.resultsPath,'obj_properties/','eigvec_pert_max_hub/');
            sol_x_filename = 'sol_x_hub';sol_t_filename = 'sol_t_hub';
            ss = obj.steady_state;
            obj.write_gephi_nodes_table(sol_path,sol_x_filename,sol_t_filename,ss);
        end
        function obj = set_graph_object(obj)
            g = graph(obj.adjacencyMatrix);
            d = distances(g);
            obj.save_var(g,obj.resultsPath,'obj_properties','graph_object');
            obj.save_var(d,obj.resultsPath,'obj_properties','shortest_distances');
        end
        function obj = set_degree_weighted_average(obj)
            obj.degree_weighted_average = sum(obj.degree_vector_weighted)/obj.N;
        end
        function set_sys_half_life_amp(obj)
            obj.sys_half_life_amp={};
            sufs = {'','_2'};
            for i2 = 1:length(sufs)
                mypath = fullfile(obj.resultsPath,'obj_properties',['solution_random_perts' sufs{i2}]);
                half_lives = [];
                for i1 = 1:5
                    At = General.load_var(fullfile(mypath,['sol_x_random_perturbations_' num2str(i1) '_norms_normed']));
                    sol_t = General.load_var(fullfile(mypath,['sol_t_random_perturbations_' num2str(i1)]));
                    At0 = At(1);
                    %                 ind = find(At < (At0 + At(end))/2,1,'first');
                    ind = find(At<.5,1,'first');
                    half_lives(end+1) = sol_t(ind);
                end
                obj.sys_half_life_amp{i2} = mean(half_lives);
            end
            obj.save_obj();
            General.save_var(obj.sys_half_life_amp,fullfile(obj.resultsPath,'obj_properties'),'sys_half_life_amp');
        end
        function set_sys_half_life_amp_2(obj)
            mypath = fullfile(obj.resultsPath,'obj_properties','solution_random_perts_2');
            half_lives = [];
            for i1 = 1:5
                At = General.load_var(fullfile(mypath,['solution_x_random_perturbations_' num2str(i1) '_norms_normed']));
                sol_t = General.load_var(fullfile(mypath,['solution_t_random_perturbations_' num2str(i1)]));
                At0 = At(1);
                ind = find(At < (At0 + At(end))/2,1,'first');
                half_lives(end+1) = sol_t(ind);
            end
            obj.sys_half_life_amp_2 = mean(half_lives);
            obj.save_obj();
        end
        function [vals_binned,vals_binned_mins,vals_binned_maxs] = set_binned_vals(~,values_vec,bins)
            n = size(bins,1);m=size(values_vec,2);
            vals_binned = zeros(n,m);vals_binned_mins = zeros(n,m);vals_binned_maxs = zeros(n,m);
            for j = 1:m
                values_vec_j = values_vec(:,j);
                for i=1:n
                    ind = bins{i,1};
                    if ~isempty(ind)
                        vals_binned(i,j) = mean(values_vec_j(ind));
                        [vals_binned_mins(i,j),~] = min(values_vec_j(ind));
                        [vals_binned_maxs(i,j),~] = max(values_vec_j(ind));
                    end
                end
            end
        end
        function constants = find_constants_binned_sets(~,binned_vals_1,binned_vals_2)
            if size(binned_vals_1,1) ~= size(binned_vals_2,1)
                disp(['find_constants_binned_sets: size of binned_vals_1 ~= size of binned_vals_2']);
            else
                constants = binned_vals_1./binned_vals_2;
            end
        end
        function obj = set_R(obj)
            obj.f_R = @(x) (-1*(obj.f_M1(x)./obj.f_M0(x)));
            disp('set_R(obj): obj.f_R = ');
            disp(obj.f_R);
        end
        function obj = set_const_functions(obj)
            syms x;
            if nargin(obj.f_M0) == 0
                c = obj.f_M0();
                obj.f_M0 = @(x) c*ones(size(x));
            end
            if nargin(obj.f_M1) == 0
                c = obj.f_M1();
                obj.f_M1 = @(x) c*ones(size(x));
            end
            if nargin(obj.f_M2) == 0
                c = obj.f_M2();
                obj.f_M2 = @(x) c*ones(size(x));
            end
            if nargin(obj.f_dM0) == 0
                c = obj.f_dM0();
                obj.f_dM0 = @(x) c*ones(size(x));
            end
            if nargin(obj.f_dM1) == 0
                c = obj.f_dM1();
                obj.f_dM1 = @(x) c*ones(size(x));
            end
            if nargin(obj.f_dM2) == 0
                c = obj.f_dM2();
                obj.f_dM2 = @(x) c*ones(size(x));
            end
        end
        function obj = set_dM0(obj)
            syms x;
            obj.f_dM0 = matlabFunction(diff(obj.f_M0,x));
            disp('set_dM0: obj.f_dM0 = ');
            disp(obj.f_dM0);
        end
        function obj = set_dM1(obj)
            syms x;
            obj.f_dM1 = matlabFunction(diff(obj.f_M1,x));
            disp('set_dM1: obj.f_dM1 = ');
            disp(obj.f_dM1);
        end
        function obj = set_dM2(obj)
            syms x;
            obj.f_dM2 = matlabFunction(diff(obj.f_M2,x));
            disp('set_dM2: obj.f_dM2 = ');
            disp(obj.f_dM2);
        end
        function obj = write_gephi_nodes_table_sparse_perturbation(obj)
            n = [1 10 100 1000];
            mypath = fullfile(obj.resultsPath,'obj_properties','random_sample_perts');
            
            for i = n
                if i<=obj.N
                    sol_x_name = ['sol_x_n_' num2str(i)];
                    x = General.load_var(fullfile(mypath,sol_x_name));
                    p = (x' - obj.steady_state')';
                    t = General.load_var(fullfile(mypath,['sol_t_n_' num2str(i)]));
                    [norms, norms_normed,thetas,H_A,tau_A,Q] = EngineClass.calc_norms_thetas(p,t);
                    General.save_var(norms,mypath,[sol_x_name '_norms']);
                    General.save_var(thetas,mypath,[sol_x_name '_thetas']);
                    General.save_var(norms_normed,mypath,[sol_x_name '_norms_normed']);
                    [m1,i1] = max(thetas);
                    p_normed = (p'./norms)';p_prop = (p'./sum(p'))';
                    [m1,i1] = min(max(p_prop'));
                    %                     i1 = find(max(p')<.01,1);
                    %stopped here
                    ts = [0 t(floor(i1/2)) t(i1)];
                    for cur_t = ts
                        filename = ['gephi_nodes_n_' num2str(i) 't_' num2str(cur_t)];
                        EngineClass.create_gephi_nodes_table(p_prop,t,cur_t,mypath,filename);
                    end
                end
            end
        end
        function set_random_perturbations_1(obj)
            if obj.num_perturbations < 5
                obj.set_num_random_perturbations();
            end
            if norm(obj.steady_state) < obj.absTol
                isOnlyPositive = true;
            else
                isOnlyPositive = false;
            end
            obj.set_random_perturbations(isOnlyPositive,true);
        end
        function set_random_perturbations_2(obj)
            if obj.num_perturbations < 5
                obj.set_num_random_perturbations();
            end
            if norm(obj.steady_state) < obj.absTol
                isOnlyPositive = true;
            else
                isOnlyPositive = false;
            end
            obj.set_random_perturbations(isOnlyPositive,false);
        end
        function solve_random_perturbations(obj)
            obj.solve(3,1,1,0,true,false);
            %             obj.solve(4,1,1,0,true,false);
        end
        function solve_random_perturbations_2(obj)
            obj.solve(4,1,1,0,true,false);
        end
        function write_norms_thetas_multi_folders(obj)
            foldernames = {'eigvec_pert_max_hub','eigvec_pert_min_hub',...
                'eigvec_pert_min_hub_1','eigvec_power_law',...
                'random_sample_perts','single_node_pert_sol',...
                'sol_pert_k_power','solution_random_perts','solution_random_perts_2','solution_eigvec_perts'};
            for foldername = foldernames
                path = fullfile(obj.resultsPath,'obj_properties',foldername{1});
                EngineClass.write_norms_thetas_multi_sol(path,obj.steady_state,obj.sys_half_life_amp,obj.mu,obj.degree_vector_weighted,'',obj.xi);
            end
        end
        function write_norms_thetas_multi_folders_2(obj)
            foldernames = {'eigvec_pert_max_hub','eigvec_pert_min_hub',...
                'eigvec_pert_min_hub_1','eigvec_power_law',...
                'random_sample_perts','single_node_pert_sol',...
                'sol_pert_k_power','solution_random_perts','solution_random_perts_2','solution_eigvec_perts'};
            for foldername = foldernames
                path = fullfile(obj.resultsPath,'obj_properties',foldername{1});
                EngineClass.write_norms_thetas_multi_sol(path,obj.steady_state,obj.sys_half_life_amp_2,obj.mu,obj.degree_vector_weighted,'_2');
            end
        end
        function obj = set_random_perturbations(obj,isOnlyPositive,isPercentOfSs)
            if isPercentOfSs
                ss = obj.steady_state;
                pert_factor = obj.perturbation_factor;
                filenamestr = 'random_perturbations';
            else
                ss = ones(size(obj.steady_state));
                %                 pert_factor = 1;
                pert_factor = min(obj.steady_state);
                filenamestr = 'random_perturbations_2';
            end
            k = obj.num_perturbations;
            n = size(ss,2);
            seed = obj.randSeed;
            rng(seed);
            %             obj.perturbations = {};
            %             obj.perturbations{1} = zeros(n,k);
            obj.random_perturbations = zeros(n,k);
            isSSzero = norm(ss)<obj.absTol;
            ss = ones(size(ss)).*max(isSSzero,ss);
            for j = 1:k
                percent = rand(1,n)*pert_factor;
                sign = rand(1,n);
                sign = (sign > (~isOnlyPositive/2))*2 - 1;
                %                 obj.perturbations{1}(:,j) = ss.*percent.*sign;
                obj.random_perturbations(:,j) = ss.*percent.*sign;
            end
            General.save_var(obj.random_perturbations,fullfile(obj.resultsPath,'obj_properties/'),filenamestr);
            %             numvec = 5;
            %             obj.perturbations{2}(:,1) = sum(obj.eigenvectors_ana(:,1:numvec),2);
            %             obj.perturbations{3}(:,1) = sum(obj.eigenvectors_asy_permuted(:,1:numvec),2);
        end
        function obj = set_random_sparse_perturbation(obj)
            num_nodes = [1 10 100 1000];
            networkPath = fullfile('networks',obj.networkName,'random_samples');
            mypath = fullfile(obj.resultsPath,'obj_properties','random_sample_perts');
            for n = num_nodes
                if n<=obj.N
                    nodes = General.load_var(fullfile(networkPath,['n_' num2str(n)]));
                    pert = zeros(obj.N,1);
                    for i = nodes
                        pert(i) = (rand*2-1)*.1*obj.steady_state(i);
                    end
                    General.save_var(pert,mypath,['n_' num2str(n)]);
                end
            end
        end
        function obj = solve_random_sparse_perts(obj)
            num_nodes = [1 10 100 1000];
            mypath = fullfile(obj.resultsPath,'obj_properties','random_sample_perts');
            stop_cond_str = 'stepEndTime >= 100*half_life';
            eps_val = 1;ss = obj.steady_state;
            for n = num_nodes
                if n<=obj.N
                    pert = General.load_var(fullfile(mypath,['n_' num2str(n)]));
                    [~, init_out, ~] = obj.check_epsilon(eps_val, pert, ss', obj.epsThreshold,...
                        obj.epsFactor, obj.init_condition_str);
                    [sol_t,sol_x] = obj.single_solve(obj.opts,obj.difEqSolver,init_out,obj.solverTimeStep,obj.maxTime,stop_cond_str);
                    General.save_var(sol_t,mypath,['sol_t_n_' num2str(n)]);
                    General.save_var(sol_x,mypath,['sol_x_n_' num2str(n)]);
                end
            end
        end
        function obj = solve_random_bin_perts(obj)
            disp('Started solve_random_bin_perts');
           pert = General.load_var(fullfile('networks',obj.networkName,'random_bin_perts','pert1'));
           mypath = fullfile(obj.resultsPath,'obj_properties','random_bin_perts');
           eps_val = 1
           stop_cond_str = 'stepEndTime >= 100*half_life';stop_after_append_all = false;
           [~, init_out, ~] = obj.check_epsilon(eps_val, pert, obj.steady_state', obj.epsThreshold,...
               obj.epsFactor, obj.init_condition_str);
           [sol_t,sol_x] = obj.single_solve(obj.opts,obj.difEqSolver,init_out,obj.solverTimeStep,obj.maxTime,stop_cond_str,stop_after_append_all);
           General.save_var(sol_t,mypath,['sol_t_pert1']);
           General.save_var(sol_x,mypath,['sol_x_pert1']);
           disp('Finished solving_random_bin_perts');
        end
        function X = set_pert_eigvec_1(obj,eigenvectors)
            n = obj.numeigenplots;
            V = eigenvectors(:,1:n)';
            B = ones(n,1)/sqrt(n);
            X = linsolve(V,B);
        end
        function obj = set_single_node_pert_struct(obj,eps_vector,folder_name)
            perts = diag(eps_vector);
            pert_inds = [];
            for i = 1:size(obj.bins,1)
                bin = obj.bins{i};
                sbin = size(bin,1);
                randvec = rand(sbin,1);
                j = 0;
                while true
                    if ~isempty(find(randvec > j*.07 & randvec < (j+1)*.07,1))
                        inds = bin(randvec >j*.07 & randvec < (j+1)*.07);
                        break
                    else
                        j=j+1;
                    end
                end
                
                pert_inds = [pert_inds;inds];
            end
            pert_inds_sorted = sort(pert_inds,'ascend');
            struct.perts_inds = pert_inds_sorted;
            perts_reduced = perts(:,pert_inds_sorted);
            struct.perts = perts_reduced;
            n = size(perts_reduced,2);
            i = 1;
            isEnd = false;
            path = fullfile('D:\',obj.resultsPath);
            while ~isEnd
                disp(['i = ' num2str(i)]);
                if i+1 >= n
                    isEnd = true;
                    j = n;
                else
                    j = i+1;
                end
                perts_curr = perts_reduced(:,i:j);
                [new_eps, inits, is_inits_legit] = check_epsilons(obj,obj.eps,perts_curr,obj.steady_state',obj.epsThreshold,...
                    obj.epsFactor,obj.init_condition);
                [solutions_t, solutions_x] = obj.multi_solve(obj.opts,obj.difEqSolver,inits,is_inits_legit,...
                    obj.solverTimeStep,obj.maxTime,obj.stop_condition_2);
                struct.inits = inits;struct.new_epsilons = new_eps;
                struct.is_inits_legit = is_inits_legit;struct.solutions_t = solutions_t;
                struct.solutions_x = solutions_x;
                
                obj.save_var(struct,path,folder_name,['pert_' num2str(i) '-' num2str(j)]);
                i = i+2;
            end
        end
        function struct = set_pert_map(obj,a,b,v1,v2,v1str,v2str)
            struct.perts = [];
            struct.factors1 = a;
            struct.factors2 = b;
            if ~isempty(a)
                struct.perts = [struct.perts,v1+a.*v2];
            end
            if ~isempty(b)
                struct.perts = [struct.perts,v2+b.*v1];
            end
            perts = struct.perts;
            [new_eps, inits, is_inits_legit] = check_epsilons(obj,obj.eps,perts,obj.steady_state',obj.epsThreshold,...
                obj.epsFactor,obj.init_condition);
            [solutions_t, solutions_x] = obj.multi_solve(obj.opts,obj.difEqSolver,inits,is_inits_legit,...
                obj.solverTimeStep,obj.maxTime,obj.stop_condition_2);
            struct.inits = inits;struct.new_epsilons = new_eps;
            struct.is_inits_legit = is_inits_legit;struct.solutions_t = solutions_t;
            struct.solutions_x = solutions_x;struct.v1str = v1str;struct.v2str = v2str;
        end
        function obj = set_opts(obj)
            obj.opts = odeset('RelTol',obj.relTol,'AbsTol',obj.absTol);
        end
        function obj = set_pert_map_struct(obj)
            a = [0 .2 .4 .6 .8 1];
            v1 = obj.eigenvectors_ana(:,1);
            v2 = obj.eigenvectors_ana(:,2);
            obj.pert_map_struct.perts = [v1+a.*v2,v2+a(1:end).*v1];
            perts = obj.pert_map_struct.perts;
            %%%%%% Adjust Epsilons
            epsVals = obj.eps;
            ss = obj.steady_state;
            init_condition_str = '~(any(init < 0 | init > 1,''all''))';
            [obj.pert_map_struct.new_epsilons, obj.pert_map_struct.inits, obj.pert_map_struct.is_inits_legit] = check_epsilons(obj,epsVals,perts,ss',obj.epsThreshold,...
                obj.epsFactor,init_condition_str);
            %%%%%%%%%%%%%%
            stop_cond = @(t,x) (size(t,1)>100 && max(abs(x(end,:) - x(end-100,:)))<obj.absTol);
            inits = obj.pert_map_struct.inits;is_inits_legit=obj.pert_map_struct.is_inits_legit;
            [obj.pert_map_struct.solutions_t, obj.pert_map_struct.solutions_x] = obj.multi_solve(obj.opts,obj.difEqSolver,inits,is_inits_legit,...
                obj.solverTimeStep,obj.maxTime,stop_cond);
        end
        function [new_epsilons, inits_out, is_inits_legit] = check_epsilons(obj,eps_vals,perts,ss,epsThreshold,...
                epsFactor,init_condition_str)
            n1 = size(eps_vals,2);
            n2 = size(perts,2);
            new_epsilons = zeros(n1,n2);
            inits_out = cell(n1,1);
            is_inits_legit = zeros(n1,n2);
            for eps_ind = 1:n1
                eps_val = eps_vals(1,eps_ind);
                for pert_ind = 1:n2
                    pert = perts(:,pert_ind);
                    [new_epsilon, init_out, is_init_legit] = obj.check_epsilon(eps_val, pert, ss, epsThreshold,...
                        epsFactor, init_condition_str);
                    new_epsilons(eps_ind,pert_ind) = new_epsilon;
                    inits_out{eps_ind}(:,pert_ind) = init_out;
                    is_inits_legit(eps_ind,pert_ind) = is_init_legit;
                end
            end
        end
        function [new_epsilon, init_out, is_init_legit] = check_epsilon(~,eps_val, pert, ss, epsThreshold,...
                epsFactor, init_condition_str)
            is_init_legit = false;
            while ~is_init_legit
                disp(['eps = ' num2str(eps_val)]);
                if abs(eps_val) < epsThreshold
                    disp('eps too small, breaking');
                    break
                end
                eps_pert = eps_val * pert;
                init = ss + eps_pert;
                init_condition = eval(init_condition_str);
                if ~init_condition
                    disp('Adjusting Epsilon');
                    if eps_val > 0
                        eps_val = -1*eps_val;
                    else
                        eps_val = -1*eps_val*epsFactor;
                    end
                    disp(['eps = ' num2str(eps_val)]);
                else
                    is_init_legit = true;
                end
            end
            new_epsilon = eps_val;
            disp(['Final eps= ' num2str(eps_val)]);
            init_out = init;
        end
        function obj = solve_degree_weighted_perts(obj)
            powers = [1,3,10];
            path = fullfile(obj.resultsPath,'obj_properties','sol_pert_k_power');
            for i = 1:length(powers)
                power = powers(i);
                pertname = ['pert_k_power_' num2str(power)];
                pert = General.load_var(fullfile('networks',obj.networkName,'perts',pertname));
                eps_val = .1;ss = obj.steady_state;epsFactor = 0.9;
                [new_epsilon, init_out, is_init_legit] = obj.check_epsilon(eps_val, pert, ss', obj.epsThreshold,...
                    epsFactor, obj.init_condition_str);
                init = init_out;% obj.steady_state' + .1*pert;
                stop_cond_str = 'stepEndTime >= 100*half_life';
                [sol_t,sol_x] = obj.single_solve(obj.opts,obj.difEqSolver,init,obj.solverTimeStep,obj.maxTime,stop_cond_str);
                General.save_var(sol_t,path,['sol_t_' pertname]);
                General.save_var(sol_x,path,['sol_x_' pertname]);
            end
        end
        
        function obj = solve_eigvec_pert_max_hub(obj,hubtype)
            mydata = load(fullfile(obj.resultsPath,'obj_properties','eigenvectors_ana_ordered_nodes.mat'),'var');eigenvectors_ana_ordered_nodes=mydata.var;
            s = size(eigenvectors_ana_ordered_nodes,1);
            switch hubtype
                case 1 % hubs
                    f = @max;row = s;foldernamestr = 'eigvec_pert_max_hub';
                    hubs = eigenvectors_ana_ordered_nodes(row,:);
                case 2 % antihubs - actually this is the eigenvector with the least mass at any particular node
                    f = @min;row = 1;foldernamestr = 'eigvec_pert_min_hub';
                    hubs = max(abs(eigenvectors_ana_ordered_nodes));
                otherwise
            end
            [~,i] = f(abs(hubs));
            eigenvectors_ana = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvectors_ana'));  eigvec_max_hub_ana = eigenvectors_ana(:,i);
            eigenvectors_asy_permuted = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvectors_asy_permuted')); eigvec_max_hub_asy= eigenvectors_asy_permuted(:,i);
            stop_cond_str = 'stepEndTime >= 100*half_life';
            init = obj.steady_state' + .1*eigvec_max_hub_ana;
            if i>obj.numeigenplots
                [sol_t_ana,sol_x_ana] = obj.single_solve(obj.opts,obj.difEqSolver,init,obj.solverTimeStep,obj.maxTime,stop_cond_str);
                obj.save_var(sol_t_ana,fullfile(obj.resultsPath,'obj_properties'),foldernamestr,'sol_t_ana');
                obj.save_var(sol_x_ana,fullfile(obj.resultsPath,'obj_properties'),foldernamestr,'sol_x_ana');
            end
            %     solution_max_hub_eigvec_ana_1.sol_t = sol_t;
            %     solution_max_hub_eigvec_ana_1.sol_x = sol_x;
            %     obj.save_var(solution_max_hub_eigvec_ana_1,obj.resultsPath,'obj_properties','solution_max_hub_eigvec_ana_1');
            init = obj.steady_state' + .1*eigvec_max_hub_asy;
            if i>obj.numeigenplots
                [sol_t_asy,sol_x_asy] = obj.single_solve(obj.opts,obj.difEqSolver,init,obj.solverTimeStep,obj.maxTime,stop_cond_str);
                obj.save_var(sol_t_asy,fullfile(obj.resultsPath,'obj_properties'),foldernamestr,'sol_t_asy');
                obj.save_var(sol_x_asy,fullfile(obj.resultsPath,'obj_properties'),foldernamestr,'sol_x_asy');
            end
            %     solution_max_hub_eigvec_asy_1.sol_t = sol_t;
            %     solution_max_hub_eigvec_asy_1.sol_x = sol_x;
            %     obj.save_var(solution_max_hub_eigvec_asy_1,obj.resultsPath,'obj_properties','solution_max_hub_eigvec_asy_1');
            [~,i] = f(obj.degree_vector_weighted);
            pert_max_hub = zeros(obj.N,1);pert_max_hub(i) = 1;
            init = obj.steady_state' + .1*pert_max_hub;
            if hubtype == 1
                [sol_t_hub,sol_x_hub] = obj.single_solve(obj.opts,obj.difEqSolver,init,obj.solverTimeStep,obj.maxTime,stop_cond_str);
                obj.save_var(sol_t_hub,fullfile(obj.resultsPath,'obj_properties'),foldernamestr,'sol_t_hub');
                obj.save_var(sol_x_hub,fullfile(obj.resultsPath,'obj_properties'),foldernamestr,'sol_x_hub');
            end
        end
        function obj = solve_single_node_pert(obj,node_id)
            pert = zeros(obj.N,1);pert(node_id) = 1/sqrt(length(node_id));
            stop_cond_str =  'a(end) < .01'; %'stepEndTime >= 100*half_life';
            eps_val = 1;ss = obj.steady_state;epsFactor = 0.9;
            [new_epsilon, init_out, is_init_legit] = obj.check_epsilon(eps_val, pert, ss', obj.epsThreshold,...
                epsFactor, obj.init_condition_str);
            foldernamestr = 'single_node_pert_sol';
            [sol_t,sol_x] = obj.single_solve(obj.opts,obj.difEqSolver,init_out,obj.solverTimeStep,obj.maxTime,stop_cond_str,false);
            obj.save_var(sol_t,fullfile(obj.resultsPath,'obj_properties'),foldernamestr,['sol_t_' EngineClass.array2str(node_id)]);
            obj.save_var(sol_x,fullfile(obj.resultsPath,'obj_properties'),foldernamestr,['sol_x_' EngineClass.array2str(node_id)]);
        end
        function obj = solve_single_node_sum_perts_batch(obj)
            [~,I] = sort(obj.degree_vector_weighted,'descend');
            node_ids = I(1);
            for i = 2:10
                node_ids(1,end+1) = I(i);
                obj.solve_single_node_pert(node_ids);
            end
        end
        function obj = solve_single_node_combs_perts_batch(obj)
            [~,I] = sort(obj.degree_vector_weighted,'descend');
            node_ids = I(1:5);
            for j = 1:5
                C = nchoosek(node_ids,j);
                for i = 1:size(C,1)
                    node_id = C(i,:);
                    obj.solve_single_node_pert(node_id);
                end
            end
        end
        function obj = solve_single_node_perts_batch(obj)
            [~,i] = sort(obj.degree_vector_weighted,'descend');
            node_ids = i(1:10);node_ids(end+1) = floor(obj.N/2);
            a=0;
            for node_id = node_ids'
                if a
                    continue
                end
                obj.solve_single_node_pert(node_id);
            end
        end
        function obj = rename_files(obj)
            p = fullfile(obj.resultsPath, 'obj_properties','single_node_pert_sol');
            contents = dir(p);
            for j = 1:length(contents)
                content = contents(j);
                fname = content.name;
                if length(fname) > 2
                    i = find(fname == 'k');
                    newname = [fname(1:i-2) '.mat'];
                    fp = fullfile(p,fname);
                    fpnew = fullfile(p,newname);
                    status = movefile(fp,fpnew);
                    if status~=1
                        disp(['status = ' num2str(status) 'for fp = ' fp]);
                    end
                end
            end
        end
        function obj = solve_eigvec_pert_max_hub_1(obj)
            obj.solve_eigvec_pert_max_hub(1);
        end
        function [sol_t, sol_x] = single_solve(obj,opts,difEqSolver,init,timeStep,maxTime,stop_cond_str,stopAfterAppendAll)
            t = 0;sol_t = [0];sol_x = [init'];
            pert0 = obj.steady_state - init';
            setBreak1 = false;setBreak2 = false; % stop while loop?
            obj.adjacencyMatrix = General.load_var(fullfile('networks',obj.networkName,'adjacency_matrix'));
            first_appendall = true;
            half_life = maxTime;a=[];first_step = true;
            while ~setBreak1 && ~setBreak2 && (~stopAfterAppendAll || first_appendall)
                % check if this step passes maxTime and set stepEndTime
                if t + timeStep >= maxTime
                    stepEndTime = maxTime;setBreak1 = true;disp('setBreak1=true');
                else
                    stepEndTime = t + timeStep;
                end
                % run solver step
                display(['stepEndTime = ' num2str(stepEndTime)]);
                odefun = @(tt,x) obj.f_M0(x)+(obj.adjacencyMatrix*obj.f_M2(x)).*obj.f_M1(x);
                [step_sol_t,step_sol_x] = difEqSolver(odefun,[t stepEndTime],init,opts);
                t = stepEndTime; init = step_sol_x(end,:);
                % append results to solution_t and solution_x
                disp(['pertnorm = ' num2str(norm(obj.steady_state - step_sol_x(end,:)))]);
                cur_sol_x = step_sol_x(end,:);
                cur_pert = obj.steady_state - cur_sol_x;
                cur_pert_prop = cur_pert./obj.steady_state;
                a(end+1) = norm(cur_pert_prop)/norm(pert0./obj.steady_state);
                disp(['a = ' num2str(a(end))]);
                appendall = (first_appendall && (a(end) < .5*norm(pert0))) || first_step;
                if ~appendall
                    sol_t(end+1,:)=step_sol_t(end,:);
                    sol_x(end+1,:)=step_sol_x(end,:);
                else
                    sol_t = [sol_t;step_sol_t];
                    sol_x = [sol_x;step_sol_x];
                    half_life = step_sol_t(end);
                    first_appendall = false;
                end
                first_step = false;
                eval(['setBreak2 = ' stop_cond_str]);
                if setBreak2
                    disp('setBreak2 = true');
                end
            end
        end
        function [solutions_t, solutions_x] = multi_solve(obj,opts,difEqSolver,inits,is_Inits_Legit,...
                timeStep,maxTime,stop_cond)
            s1 = size(inits);s2 = size(inits{1});
            solutions_t = cell(s1);solutions_x = cell(s1);
            for epsInd = 1:s1(1)
                inits_current = inits{epsInd};n1 = size(inits_current,2);
                for initInd = 1:n1
                    if is_Inits_Legit(epsInd,initInd)
                        disp(['epsInd = ' num2str(epsInd) ', initInd = ' num2str(initInd)]);
                        init = inits_current(:,initInd);
                        [sol_t, sol_x] = obj.single_solve(opts,difEqSolver,init,timeStep,maxTime,stop_cond);
                        solutions_t{epsInd,initInd} = sol_t;
                        solutions_x{epsInd,initInd} = sol_x;
                    end
                end
            end
        end
        function obj = solve(obj,pertType,epsIndStart,pertIndStart,isAdjustEpsilonType,isBreakAfterHalfLife,isResetXSol)
            obj.adjacencyMatrix = General.load_var(fullfile('networks',obj.networkName,'adjacency_matrix'));
            %Solve the system
            tic;
            %nn=10;
            nn=5; % number of eigenvectors
            ss = obj.steady_state;
            perts = 0;
            %ssInd = useSS + 1;
            eps_vals = obj.eps;
            switch pertType
                case 1 % initial state is obj.inititalValues
                    ss = 0;
                    eps_vals = 1;
                    perts = obj.initialValues;
                    sol_t_var_str = 'sol_t_1';
                    sol_x_var_str = 'sol_x_1';
                    stop_cond_str = 'size(sol_t_1{1,pertInd,epsInd},1)>100 && max(abs(sol_x_1{1,pertInd,epsInd}(end,:)-sol_x_1{1,pertInd,epsInd}(end-100,:)))<obj.absTol';
                    solution_folder_name = 'solution';
                case 2 % perturbations are eigenvectors
                    %                     perts = [obj.eigenvectors_ana(:,1:nn), obj.eigenvectors_asy_permuted(:,1:nn),...
                    %                         sum(obj.eigenvectors_ana(:,1:nn),2),sum(obj.eigenvectors_asy_permuted(:,1:nn),2)...
                    %                         obj.pert_eigvec_ana_1,obj.pert_eigvec_asy_1];
                    obj.eigenvectors_ana = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvectors_ana'));
                    obj.eigenvectors_asy_permuted = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvectors_asy_permuted'));
                    perts = [obj.eigenvectors_ana(:,1:nn)];%, obj.eigenvectors_asy_permuted(:,1:nn)];
                    sol_t_var_str = 'sol_t_eigvec';
                    sol_x_var_str = 'sol_x_eigvec';
                    stop_cond_str = 'max(abs(sol_x_eigvec{1,pertInd,epsInd}(end,:)-obj.steady_state))<obj.absTol*100';
                    solution_folder_name = 'solution_eigvec_perts';
                case 3 % perturbations are random and small relative to steady state
                    %                     perts = obj.perturbations{1};
                    perts = General.load_var(fullfile(obj.resultsPath,'obj_properties','random_perturbations'));
                    sol_t_var_str = 'sol_t_random_perturbations';
                    sol_x_var_str = 'sol_x_random_perturbations';
                    stop_cond_str = 'max(abs(sol_x_random_perturbations{1,pertInd,epsInd}(end,:)-obj.steady_state))<obj.absTol*100';
                    solution_folder_name = 'solution_random_perts';
                case 4 % perturbations are random and small not relative to steady state
                    perts = General.load_var(fullfile(obj.resultsPath,'obj_properties','random_perturbations_2'));
                    sol_t_var_str = 'sol_t_random_perturbations';
                    sol_x_var_str = 'sol_x_random_perturbations';
                    stop_cond_str = 'max(abs(sol_x_random_perturbations{1,pertInd,epsInd}(end,:)-obj.steady_state))<obj.absTol*100';
                    solution_folder_name = 'solution_random_perts_2';
                otherwise
            end
            numperts = size(perts,2);
            numeps = length(eps_vals);
            %             eval([sol_t_var_str ' = {};']);
            %             eval([sol_x_var_str ' = {};']);
            
            %%%%%% Adjust Epsilons
            switch isAdjustEpsilonType
                case 0 %no condition
                    init_condition_str = 'true';
                case 1 %SIS
                    init_condition_str = '~(any(init < 0 | init > 1,''all''))';
                case 2 %REG
                    init_condition_str = '~(any(init < 0))';
                otherwise
            end
            
            epsFactor = .9;
            [new_epsilons, inits_out, is_inits_legit] = check_epsilons(obj,eps_vals,perts,ss',obj.epsThreshold,...
                epsFactor,init_condition_str);
            obj.eps_adjusted = new_epsilons;
            inits = inits_out;
            obj.isInitsLegit = is_inits_legit;
            
            %%%%%%%%%%%%%%
            
            %%%%%%
            eval([sol_t_var_str ' = {};']);
            eval([sol_x_var_str ' = {};']);
            
            opts = odeset('RelTol',obj.relTol,'AbsTol',obj.absTol);
            clear odefun;
            odefun = @(tt,x) (obj.f_M0(x) + (obj.adjacencyMatrix*obj.f_M2(x)).*obj.f_M1(x));
            %             addpath(fullfile(obj.resultsPath,'obj_properties'));
            for epsInd = epsIndStart:numeps
                disp(['epsInd = ' num2str(epsInd)]);
                disp(['eps_var = ' num2str(obj.eps_adjusted(epsInd))]);
                for pertInd = pertIndStart:numperts
                    if obj.isInitsLegit(epsInd,pertInd)
                        disp(['pertInd = ' num2str(pertInd)]);
                        init = inits{1,epsInd}(:,pertInd);
                        init = init';
                        t = 0;
                        %%%%%%%%%%%%%%
                        eval([sol_t_var_str '{1,pertInd,epsInd} = 0;']);
                        eval([sol_x_var_str '{1,pertInd,epsInd} = init;']);
                        %%%%%%%%%%%%%%
                        setBreak1 = false;setBreak2 = false; % stop while loop?
                        first_appendall = true;
                        if pertType~=1
                            pert0 = init-obj.steady_state;a=[];
                        end
                        aa=[];first_step = true;
                        while ~setBreak1 && ~setBreak2
                            % check if this step passes maxTime and set stepEndTime
                            %                             if t + obj.solverTimeStep >= obj.maxTime
                            %                                 stepEndTime = obj.maxTime;
                            %                                 setBreak1 = true;
                            %                                 disp('setBreakk1 = true');
                            %                             else
                            stepEndTime = t + obj.solverTimeStep;
                            %                             end
                            % run solver step
                            display(['stepEndTime = ' num2str(stepEndTime)]);
                            [sol_t,sol_x] = obj.difEqSolver(odefun,[t stepEndTime],init,opts);
                            t = stepEndTime;
                            init = sol_x(end,:);
                            % append results to solution_t and solution_x
                            %                             eval([sol_t_var_str '{1,pertInd,epsInd}(end+1:end+length(sol_t)-1,1)=sol_t(2:end,:);']);
                            %                             eval([sol_x_var_str '{1,pertInd,epsInd}(end+1:end+size(sol_x,1)-1,:)=sol_x(2:end,:);']);
                            
                            if pertType ~= 1
                                cur_sol_x = sol_x(end,:);
                                cur_pert = obj.steady_state - cur_sol_x;
                                cur_pert_prop = cur_pert./obj.steady_state;
                                a(end+1) = norm(cur_pert_prop)/norm(pert0./obj.steady_state);
                                disp(['a = ' num2str(a(end))]);
                                appendall = (first_appendall && (a(end) < .5)) || first_step; %(norm(obj.steady_state - sol_x(end,:)) < .5*norm(pert0));
                                if isBreakAfterHalfLife && a(end) < .01
                                    setBreak1 = true;
                                    disp('isBreakAfterHalfLife true, setBreak1 = true');
                                end
                            else
                                appendall = false;
                            end
                            if ~appendall
                                eval(['sol_x_num_rows = size(' sol_x_var_str '{1,pertInd,epsInd},1)']);
                                if isResetXSol && sol_x_num_rows > 1000
                                    eval([sol_t_var_str '{1,pertInd,epsInd}(1:end-101)=[];']);
                                    eval([sol_x_var_str '{1,pertInd,epsInd}(1:end-101,:)=[];']);
                                end
                                eval([sol_t_var_str '{1,pertInd,epsInd}(end+1,1)=sol_t(end,:);']);
                                eval([sol_x_var_str '{1,pertInd,epsInd}(end+1,:)=sol_x(end,:);']);
                            else
                                eval([sol_t_var_str '{1,pertInd,epsInd}(end+1:end+length(sol_t)-1,1)=sol_t(2:end,:);']);
                                eval([sol_x_var_str '{1,pertInd,epsInd}(end+1:end+size(sol_x,1)-1,:)=sol_x(2:end,:);']);
                                first_appendall = false;

                            end
                            first_step = false;
                            if pertType == 1
                                aa(end+1)=max(abs(sol_x_1{1,pertInd,epsInd}(end,:)-sol_x_1{1,pertInd,epsInd}(end-1,:)));
                                plot(aa,'.');
                            end
                            eval(['setBreak2 = ' stop_cond_str ';']);
                            if setBreak2
                                disp('setBreak2 = true');
                            end
                        end
                    else
                        disp(['epsInd = ' num2str(epsInd) ', pertInd = ' num2str(pertInd)]);
                        disp(['Did not solve, no valid inits' ]);
                    end
                    eval(['sol_t_cur = ' sol_t_var_str '{1,pertInd,epsInd};']);
                    eval(['sol_x_cur = ' sol_x_var_str '{1,pertInd,epsInd};']);
                    General.save_var(sol_t_cur, fullfile(obj.resultsPath,'obj_properties',solution_folder_name),[sol_t_var_str '_' num2str(pertInd)]);
                    General.save_var(sol_x_cur, fullfile(obj.resultsPath,'obj_properties',solution_folder_name),[sol_x_var_str '_' num2str(pertInd)]);
                    
                end
            end
            toc;
            %             str = ['EngineClass.save_var(' sol_x_var_str ',obj.resultsPath,'obj_properties','J_ana')];
            %             eval();
            if pertType == 1
                %                 obj.solution_t = obj.solution_t{1};
                %                 obj.solution_x = obj.solution_x{1};
                obj.set_steady_state();
                obj.set_M2_i_bigodot();
                obj.set_steady_state_calculated();
                obj.set_Dii_ana();
                obj.set_Wij_ana();
                obj.set_J_ana();
                obj.set_J_ana_degree_vector_weighted();
                obj.set_J_ana_k_inn();
                obj.set_eig_ana();
                obj.set_eig_ana_ordered_nodes();
                obj.Dii_anabinned = obj.set_binned_vals(obj.Dii_ana,obj.bins);
                Wij_ana = General.load_var(fullfile(obj.resultsPath,'obj_properties','Wij_ana'));
                obj.Wij_anabinned = obj.set_binned_vals(Wij_ana(obj.adjacencyMatrix>0),obj.binsij);
                %                 obj.C_D_set = obj.find_constants_binned_sets(obj.Dii_anabinned,obj.Dii_asybinned);
                %                 obj.C_W_set = obj.find_constants_binned_sets(obj.Wij_anabinned,obj.Wij_asybinned);
                %                 [obj.C_D_v3, obj.guesses1_v3, obj.guesses2_v3, obj.errors1_v3, obj.errors2_v3, obj.eigvals_v3,obj.eigvecs_v3] = obj.find_best_C_D(obj.eigenvalues_ana(1,1),1,5,0.1);
                %                 [obj.C_D_v4, obj.guesses1_v4, obj.guesses2_v4, obj.errors1_v4, obj.errors2_v4, obj.eigvals_v4,obj.eigvecs_v4] = obj.find_best_C_D(obj.eigenvalues_ana(1,1),1,5,0.5,obj.permutation_eigvec_ana2asy);
                obj.set_eigvec_comparison_mats(true,false);
                obj.set_permutation_eigvec_ana2asy(0);
                obj.set_eig_asy_permuted();
                eig_ana = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvectors_ana'));
                eig_asy_perm = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvectors_asy_permuted'));
                %                 [eigvals_asy_v2_set, eigvecs_asy_v2_set, Dii_asy_set, Wij_asy_set] = obj.set_eig_v2_sets(obj.Dii_asy, obj.Wij_asy, obj.C_D_set, obj.C_W_set, obj.numeigen,1);
                for i = 1:obj.numeigenplots
                    obj.binseigvecana{i,1} = obj.set_bins_generic(obj.numbins,abs(eig_ana(:,i)),1e-13,true(obj.N,1));
                    obj.binseigvecasy_permuted{i,1} = obj.set_bins_generic(obj.numbins,abs(eig_asy_perm(:,i)),1e-13,true(obj.N,1));
                    obj.eigenvectors_ana_binned_self{i,1} = obj.set_binned_vals(abs(eig_ana(:,i)),obj.binseigvecana{i,1});
                    obj.eigenvectors_asy_permuted_binned_self{i,1} = obj.set_binned_vals(abs(eig_asy_perm(:,i)),obj.binseigvecasy_permuted{i,1});
                end
                obj.eigenvectors_ana_binned_k = obj.set_binned_vals(eig_ana,obj.bins);
                obj.eigenvectors_ana_binned_kinn = obj.set_binned_vals(eig_ana,obj.binskinn);
                obj.eigenvectors_asy_permuted_binned_k = obj.set_binned_vals(eig_asy_perm,obj.bins);
                obj.eigenvectors_asy_permuted_binned_kinn = obj.set_binned_vals(eig_asy_perm,obj.binskinn);
                obj.set_eigvec_comparison_mats(false,true);
                obj.set_weighted_dot_products();
                obj.pert_eigvec_ana_1 = obj.set_pert_eigvec_1(eig_ana);
                obj.pert_eigvec_asy_1 = obj.set_pert_eigvec_1(eig_asy_perm);
                if norm(obj.steady_state) < obj.absTol
                    isOnlyPositive = true;
                else
                    isOnlyPositive = false;
                end
                obj.set_random_perturbations(isOnlyPositive,1);
                obj.set_random_perturbations(isOnlyPositive,0);
                %                 obj.set_eigvec_comparison_mats2();
            elseif pertType==2
                obj.split_solution_eigvec();
            end
            obj.save_obj();
        end
        function [node_half_life_1, node_half_life_2, prec1, prec2] = find_node_half_life(obj,sol_t,sol_x,node_num)
            node_sol_x = sol_x(:,node_num);
            node_pert_x = node_sol_x - obj.steady_state(node_num);
            node_half_val_1 = (node_pert_x(1) + node_pert_x(end))/2;
            node_half_val_index_1 = find(abs(node_pert_x) < abs(node_half_val_1),1,'first');
            node_half_life_1 = sol_t(node_half_val_index_1);
            prec1 = node_half_val_1 - node_pert_x(node_half_val_index_1);
            
            [~,I] = max(abs(node_pert_x));
            M = node_pert_x(I);
            node_half_val_2 = (M + node_pert_x(end))/2;
            node_pert_x_adj = node_pert_x;
            node_pert_x_adj(1:I-1) = 999999;
            node_half_val_index_2 = find(abs(node_pert_x_adj) < abs(node_half_val_2),1,'first');
            node_half_life_2 = sol_t(node_half_val_index_2) - sol_t(I);
            prec2 = node_half_val_2 - node_pert_x(node_half_val_index_2);
        end
        function [node_half_lives_1, node_half_lives_2, prec1, prec2] = find_node_half_lives(obj,sol_t,sol_x)
            n = size(sol_x,2);
            node_half_lives_1 = zeros(n,1);
            node_half_lives_2 = zeros(n,1);
            prec1 = zeros(n,1);
            prec2 = zeros(n,1);
            for node_num = 1:n
                [node_half_lives_1(node_num,1), node_half_lives_2(node_num,1), prec1(node_num,1), prec2(node_num,1)] = find_node_half_life(obj,sol_t,sol_x,node_num);
            end
        end
        function obj = plot_node_half_lives(obj)
            %% fig12-*
            name = 'fig12-1';
            figdesc = 'Node Half Lives vs Degree';
            f = figure('Name',name,'NumberTitle','off');
            hax = axes;
            S = obj.node_half_lives_struct;
            yvals = [S.nhl2{2}{1} S.nhl2{2}{2} S.nhl2{2}{3} S.nhl2{2}{4}];
            plot(obj.degree_vector_weighted,yvals,'o');
            hold on;
            l = obj.eigenvalues_ana(1:4);
            l = -1./l * log(2);
            plot(hax.XLim',[l,l]','-');
            legend('$pert = v_{1,ana}$','$pert = v_{2,ana}$','$pert = v_{3,ana}$',...
                '$pert = v_{4,ana}$','$ln2/\lambda_1$','$ln2/\lambda_2$',...
                '$ln2/\lambda_3$','$ln2/\lambda_4$','interpreter','latex');
            xlabel('k_i (weighted)');
            ylabel('node half life');
            title({[name ' ' obj.scenarioName];obj.desc;figdesc},'interpreter','latex');
            obj.save_fig(f,name);
        end
        function obj = set_node_half_lives_struct(obj)
            S.names = {'x_0 = rand','pert_0 = \varepsilon * v_{i,ana}',...
                'pert_0 = \varepsilon * v_{i,asy}', 'pert_0 = v_{1,ana} + v_{2,ana}',...
                };
            sols_t = {{obj.solution_t},obj.solution_t_eigvecana(1,1:4),...
                obj.solution_t_eigvecasy(1,1:5),obj.pert_map_struct.solutions_t(1,6)};
            sols_x = {{obj.solution_x},obj.solution_x_eigvecana(1,1:4),...
                obj.solution_x_eigvecasy(1,1:5),obj.pert_map_struct.solutions_x(1,6)};
            n1 = size(sols_t,2);
            for i=1:n1
                sols_ti = sols_t{1,i}; sols_xi = sols_x{1,i};
                n2 = size(sols_ti,2);
                for j = 1:n2
                    disp(['i= ' num2str(i)]);
                    disp(['j= ' num2str(j)]);
                    sol_t = sols_ti{1,j};sol_x = sols_xi{1,j};
                    [nhl1,nhl2,p1,p2]=obj.find_node_half_lives(sol_t,sol_x);
                    S.nhl1{i}{j} = nhl1; S.nhl2{i}{j} = nhl2; S.prec1{i}{j} = p1; S.prec2{i}{j} = p2;
                end
            end
            obj.node_half_lives_struct = S;
        end
        function obj = split_solution_eigvec(obj)
            sol_t = obj.solution_t_eigvec;
            sol_x = obj.solution_x_eigvec;
            [numeps,numvecs] = size(sol_x);
            obj.solution_t_eigvec = {};
            obj.solution_x_eigvec = {};
            for e = 1:numeps
                for v = 1:numvecs
                    if v<=obj.numeigenplots
                        if ~isempty(sol_t{1,v,e})
                            obj.solution_t_eigvecana{1,v,e} = sol_t{1,v,e};
                            obj.solution_x_eigvecana{1,v,e} = sol_x{1,v,e};
                        end
                    elseif v<=obj.numeigenplots*2
                        if ~isempty(sol_t{1,v,e})
                            obj.solution_t_eigvecasy{1,v-obj.numeigenplots,e} = sol_t{1,v,e};
                            obj.solution_x_eigvecasy{1,v-obj.numeigenplots,e} = sol_x{1,v,e};
                        end
                    else
                        obj.solution_t_eigvec{1,v-obj.numeigenplots*2,e} = sol_t{1,v,e};
                        obj.solution_x_eigvec{1,v-obj.numeigenplots*2,e} = sol_x{1,v,e};
                    end
                end
            end
        end
        function obj = set_steady_state(obj)
            sol_x = General.load_var(fullfile(obj.resultsPath,'obj_properties','solution','sol_x_1_1'));
            obj.steady_state = sol_x(end,:);
            disp('set_steady_state(obj): obj.steady_state = ');
            General.save_var(obj.steady_state,fullfile(obj.resultsPath,'obj_properties'),'steady_state');
            %disp(obj.steady_state);
        end
        function obj = set_M2_i_bigodot(obj)
            obj.M2_i_bigodot = zeros(size(obj.initialValues));
            M2_x = obj.f_M2(obj.steady_state);
            neighbor_sum_M2_x = obj.adjacencyMatrix * M2_x';
            obj.M2_i_bigodot = neighbor_sum_M2_x ./ obj.degree_vector_weighted;
            disp('set_M2_i_bigdot(obj): obj.M2_i_bigdot = ');
            %disp(obj.M2_i_bigodot);
        end
        function obj = set_steady_state_calculated(obj)
            syms x;
            R_inv = finverse(obj.f_R(x));
            x = 1./(obj.degree_vector_weighted .* obj.M2_i_bigodot);
            obj.steady_state_calculated = double(subs(R_inv));
            disp('set_steady_state_calculated(obj): obj.steady_state_calculated = ');
            %disp(obj.steady_state_calculated);
        end
        function obj = set_Dii_ana(obj)
            x = obj.steady_state';
            obj.Dii_ana = double(obj.f_dM0(x)) + ...
                double(obj.f_dM1(x)).*double(obj.adjacencyMatrix*obj.f_M2(x));
            disp('set_Dii_ana(obj): obj.Dii_ana = ');
            %disp(obj.Dii_ana);
        end
        function obj = set_Wij_ana(obj)
            x = obj.steady_state;
            Wij_ana = obj.adjacencyMatrix .* (double(obj.f_M1(x))'*...
                double(obj.f_dM2(x)));
            disp('set_Wij_ana(obj): obj.Wij_ana = ');
            General.save_var(Wij_ana,fullfile(obj.resultsPath,'obj_properties'),'Wij_ana');
            %disp(obj.Wij_ana);
        end
        function obj = set_J_ana(obj)
            Wij_ana = General.load_var(fullfile(obj.resultsPath,'obj_properties','Wij_ana'));
            J_ana = EngineClass.compute_J(obj.Dii_ana,Wij_ana,1,1);
            obj.save_var(J_ana,obj.resultsPath,'obj_properties','J_ana');
        end
        function obj = set_J_asy(obj)
            Wij_asy = General.load_var(fullfile(obj.resultsPath,'obj_properties','Wij_asy'));
            J_asy = EngineClass.compute_J(obj.Dii_asy,Wij_asy,1,1);
            obj.save_var(J_asy,obj.resultsPath,'obj_properties','J_asy');
        end
        function obj = set_J_ana_degree_vector_weighted(obj)
            mydata=load(fullfile(obj.resultsPath,'obj_properties','J_ana'));J_ana=mydata.var;
            J_ana_degree_vector_weighted = EngineClass.compute_degree_vector_weighted(J_ana);
            obj.save_var(J_ana_degree_vector_weighted,obj.resultsPath,'obj_properties','J_ana_degree_vector_weighted');
        end
        function obj = set_J_asy_degree_vector_weighted(obj)
            mydata=load(fullfile(obj.resultsPath,'obj_properties','J_asy'));J_asy=mydata.var;
            J_asy_degree_vector_weighted = EngineClass.compute_degree_vector_weighted(J_asy);
            obj.save_var(J_asy_degree_vector_weighted,obj.resultsPath,'obj_properties','J_asy_degree_vector_weighted');
        end
        function obj = set_J_asy_k_inn(obj)
            mydata=load(fullfile(obj.resultsPath,'obj_properties','J_asy'));J_asy=mydata.var;
            J_asy_k_inn = EngineClass.compute_k_inn(J_asy);
            obj.save_var(J_asy_k_inn,obj.resultsPath,'obj_properties','J_asy_k_inn');
        end
        function obj = set_J_ana_k_inn(obj)
            mydata=load(fullfile(obj.resultsPath,'obj_properties','J_ana'));J_ana=mydata.var;
            J_ana_k_inn = EngineClass.compute_k_inn(J_ana);
            obj.save_var(J_ana_k_inn,obj.resultsPath,'obj_properties','J_ana_k_inn');
        end
        function obj = set_N(obj)
            obj.N = size(obj.adjacencyMatrix,1);
        end
        function obj = set_knn(obj)
            tmp = obj.adjacencyMatrix*obj.degree_vector_weighted;
            tmp = tmp./obj.degree_vector_weighted;
            obj.k_nn = sum(tmp)/obj.N;
        end
        function obj = set_kinn(obj)
            mask = ones(obj.N) - eye(obj.N);
            tmp = obj.adjacencyMatrix .* mask;
            tmp = tmp * obj.degree_vector_weighted;
            obj.ki_nn = tmp ./ obj.degree_vector_weighted;
        end
        function obj = set_Dii_asy(obj)
            obj.Dii_asy = -obj.k_nn^obj.eta*obj.degree_vector_weighted.^obj.mu;
        end
        function obj = set_Wij_asy(obj)
            Wij_asy = obj.adjacencyMatrix .* (obj.degree_vector_weighted.^obj.nu * obj.degree_vector_weighted'.^obj.rho);
            General.save_var(Wij_asy,fullfile(obj.resultsPath,'obj_properties'),'Wij_asy');
        end
        function [C_D, guesses1, guesses2, errors1, errors2, eigvals, eigvecs] = find_best_C_D(obj,real_eigval_1,guess_1_init,guess_2_init,req_acc)
            Wij_ana = General.load_var(obj.resultsPath,'obj_properties','Wij_ana');
            Wij_asy = General.load_var(obj.resultsPath,'obj_properties','Wij_asy');
            guess1 = guess_1_init;
            guess2 = guess_2_init;
            recalc1 = true;recalc2=true;
            numeigen = 200;
            guesses1 = [];guesses2 = [];
            errors1 = []; errors2 = [];
            while true
                guesses1(end+1) = guess1;
                guesses2(end+1) = guess2;
                disp(['guess1 = ' num2str(guess1) ', guess2 = ' num2str(guess2)]);
                if recalc1
                    [eigvals_1, eigvecs_1, ~, ~] = obj.set_eig(obj.Dii_asy, Wij_asy, guess1, 1, numeigen, 1);
                end
                if recalc2
                    [eigvals_2, eigvecs_2, ~, ~] = obj.set_eig(obj.Dii_asy, Wij_asy, guess2, 1, numeigen, 1);
                end
                error1 = eigvals_1(1,1) - real_eigval_1;
                error2 = eigvals_2(1,1) - real_eigval_1;
                disp(['error1 = ' num2str(error1) ', error2 = ' num2str(error2)]);
                errors1(end+1) = error1;errors2(end+1) = error2;
                if abs(error1) < req_acc
                    disp('breaking');
                    C_D = guess1;
                    eigvals = eigvals_1;
                    eigvecs = eigvecs_1;
                    break
                elseif abs(error2) < req_acc
                    disp('breaking');
                    C_D = guess2;
                    eigvals = eigvals_2;
                    eigvecs = eigvecs_2;
                    break
                end
                if (error1 > 0 && error2 > 0) || (error1 < 0 && error2 < 0)
                    if abs(error2) > abs(error1)
                        if guess2 > guess1
                            guess2 = guess1 - (guess2 - guess1);
                        elseif guess2 < guess1
                            guess2 = guess1 + (guess1 - guess2);
                        end
                        recalc2 = true;recalc1 = false;
                    elseif abs(error1) > abs(error2)
                        if guess1 > guess2
                            guess1 = guess2 - (guess1 - guess2);
                        elseif guess1 < guess2
                            guess1 = guess2 + (guess2 - guess1);
                        end
                        recalc1 = true;recalc2 = false;
                    end
                elseif (error1 > 0 && error2 < 0) || (error1 < 0 && error2 > 0)
                    if abs(error1) > abs(error2)
                        guess1 = (guess1 + guess2)/2;
                        recalc1 = true;recalc2=false;
                    elseif abs(error2) > abs(error1)
                        guess2 = (guess1 + guess2)/2;
                        recalc2 = true;recalc1=false;
                    end
                end
            end
        end
        function [eigvals, eigvecs, Dii_v2, Wij_v2] = set_eig(~, Dii, Wij, C_D, C_W, numeigen, iscalceigen)
            v = []; d = [];
            J = Wij*C_W;
            J = J - diag(diag(J));
            D = Dii*C_D;
            J = J + diag(D);
            if iscalceigen
                [v,d] = eigs(J,min(size(J,1),numeigen),'largestreal');
            end
            eigvals = diag(d);
            eigvecs = v;
            Dii_v2 = Dii*C_D;
            Wij_v2 = Wij*C_W;
            disp('Done');
        end
        function obj = set_eig_ana(obj)
            Wij_ana = General.load_var(fullfile(obj.resultsPath,'obj_properties','Wij_ana'));
            [obj.eigenvalues_ana, eigenvectors_ana, ~, ~] = obj.set_eig(obj.Dii_ana,Wij_ana,1,1,obj.numeigen,1);
            obj.save_var(obj.eigenvalues_ana,obj.resultsPath,'obj_properties','eigenvalues_ana');
            obj.save_var(eigenvectors_ana,obj.resultsPath,'obj_properties','eigenvectors_ana');
        end
        function obj = set_eig_ana_ordered_nodes(obj)
            eigvec_ana = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvectors_ana'));
            [ind_bins_var,binned_vals] = EngineClass.set_bins_percentiles(1,obj.degree_vector_weighted);
            [ordered_eigs] = EngineClass.reorder_eigvecs_nodes(eigvec_ana, ind_bins_var);
            obj.save_var(ordered_eigs,obj.resultsPath,'obj_properties','eigenvectors_ana_ordered_nodes');
        end
        function obj = set_eig_asy(obj)
            Wij_asy = General.load_var(fullfile(obj.resultsPath,'obj_properties','Wij_asy'));
            [obj.eigenvalues_asy, eigenvectors_asy, ~, ~] = obj.set_eig(obj.Dii_asy,Wij_asy,1,1,obj.numeigen,1);
            obj.save_var(obj.eigenvalues_asy,obj.resultsPath,'obj_properties','eigenvalues_asy');
            obj.save_var(eigenvectors_asy,obj.resultsPath,'obj_properties','eigenvectors_asy');
        end
        function [eigvals_asy_v2_set, eigvecs_asy_v2_set, Dii_asy_v2_set, Wij_asy_v2_set, Dii_asy_v2_binned_set, Wij_asy_v2_binned_set] = set_eig_v2_sets(obj, Dii_asy, Wij_asy, C_D_set, C_W_set, numeigen, iscalceigen)
            num_constants = size(C_D_set,1);
            Wij_asy = General.load_var(obj.resultsPath,'obj_properties','Wij_asy');
            Dii_asy_v2_set = cell(num_constants,1);Wij_asy_v2_set = cell(num_constants,1);
            eigvals_asy_v2_set = cell(num_constants,1); eigvecs_asy_v2_set = cell(num_constants,1);
            Dii_asy_v2_binned_set = cell(num_constants,1);Wij_asy_v2_binned_set = cell(num_constants,1);
            for i = 1:num_constants
                disp(['set_eig_v2: i = ' num2str(i)]);
                C_D = C_D_set(num_constants - i + 1);
                C_W = C_W_set(i);
                [a,b,c,d] = obj.set_eig(obj.Dii_asy,Wij_asy,C_D,C_W,numeigen, iscalceigen);
                if iscalceigen
                    eigvals_asy_v2_set{i,1}=a; eigvecs_asy_v2_set{i,1}=b;
                end
                Dii_asy_v2_set{i,1}=c; Wij_asy_v2_set{i,1}=d;
                Dii_asy_v2_binned_set{i,1} = obj.set_binned_vals(c,obj.bins);
                Wij_asy_v2_binned_set{i,1} = obj.set_binned_vals(d(obj.adjacencyMatrix>0),obj.binsij);
            end
            if iscalceigen
                obj.save_var(eigvals_asy_v2_set,obj.resultsPath,'obj_properties','eigvals_asy_v2_set');
                obj.save_var(eigvecs_asy_v2_set,obj.resultsPath,'obj_properties','eigvecs_asy_v2_set');
            end
            obj.save_var(Dii_asy_v2_set, obj.resultsPath,'obj_properties','Dii_asy_v2_set');
            obj.save_var(Wij_asy_v2_set, obj.resultsPath,'obj_properties','Wij_asy_v2_set');
            obj.save_var(Dii_asy_v2_binned_set,obj.resultsPath,'obj_properties','Dii_asy_v2_binned_set');
            obj.save_var(Wij_asy_v2_binned_set,obj.resultsPath,'obj_properties','Wij_asy_v2_binned_set');
        end
        %         function obj = set_eig_ana(obj)
        %             J = obj.Wij_ana;
        %             J = J - diag(diag(J));
        %             J = J + diag(obj.Dii_ana);
        %             %[v,d] = eig(J);
        %             [v,d] = eigs(J,min(size(J,1),obj.numeigen),'largestreal');
        %             obj.eigenvalues_ana = diag(d);
        %             obj.eigenvectors_ana = v;
        %             disp('Done');
        %         end
        %         function obj = set_eig_asy(obj)
        %             J = obj.Wij_asy;
        %             J = J - diag(diag(J));
        %             J = J + diag(obj.Dii_asy);
        %             %[v,d] = eig(J);
        %             [v,d] = eigs(J,min(size(J,1),obj.numeigen),'largestreal');
        %             obj.eigenvalues_asy = diag(d);
        %             obj.eigenvectors_asy = v;
        %             disp('Done');
        %         end
        function [M1, M2, M3] = compute_comparison_matrix(obj,vectors_1,vectors_2,isSaveM2,varNameStr)
            tic;
            n1 = size(vectors_1,2);
            n2 = size(vectors_2,2);
            M1 = zeros(n1,n2);M2 = M1;M3 = M1;
            parfor i = 1:n1
                if mod(i,100) == 0
                    disp(['i = ' num2str(i)]);
                end
                ui = vectors_1(:,i);
                uimat = ui + zeros(size(vectors_2));
                m2 = abs(dot(uimat,vectors_2));M2(i,:) = m2;
                m3 = rad2deg(acos(m2)); M3(i,:) = m3;
                %                 parfor j = 1:n2
                %                     vj = vectors_2(:,j);
                %                     % L2 norm of difference between vectors
                %                     %m1 = norm(ui - vj);
                %                     % dot product of vectors
                %                     m2 = abs(dot(ui,vj));
                %                     % angle between vectors
                %                     m3 = rad2deg(acos(m2));
                %                     %M1(i,j) = m1;
                %                     M2(i,j) = m2; M3(i,j) = m3;
                %                 end
            end
            toc;
            if isSaveM2
                obj.save_var(M2,obj.resultsPath,'obj_properties',varNameStr);
            end
        end
        function obj = set_eigvec_comparison_mats3(obj)
            mydata = load(fullfile(obj.resultsPath,'obj_properties','eigvecs_asy_v2_set.mat'),'var');eigvecs_asy_v2_set=mydata.var;
            [~, ~, ~] = obj.compute_comparison_matrix(obj.eigenvectors_asy,eigvecs_asy_v2_set{14,1},true,'eigvec_dot_comparison_mat_asy2asy-v2-14');
            [~, ~, ~] = obj.compute_comparison_matrix(obj.eigenvectors_asy,eigvecs_asy_v2_set{15,1},true,'eigvec_dot_comparison_mat_asy2asy-v2-15');
            [~, ~, ~] = obj.compute_comparison_matrix(obj.eigenvectors_asy,obj.eigvecs_v3,true,'eigvec_dot_comparison_mat_asy2asy-v3');
        end
        function obj = set_eigvec_comparison_mats(obj,computeRegular,computePermuted)
            if computeRegular
                vasy = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvectors_asy'));
            elseif computePermuted
                vasy = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvectors_asy_permuted'));
            end
            vana = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvectors_ana'));
            [dist, absdot, angles] = obj.compute_comparison_matrix(vana,vasy,false,'');
            
            if computeRegular
                %                 obj.eigvec_dist_comparison_mat_ana2asy = dist;
                %                 obj.eigvec_angle_comparison_mat_ana2asy = angles;
                %                 obj.eigvec_dot_comparison_mat_ana2asy = absdot;
                EngineClass.save_var(absdot,obj.resultsPath,'obj_properties','eigvec_dot_comparison_mat_ana2asy');
            elseif computePermuted
                %                 obj.eigvec_dist_comparison_mat_ana2asy_permuted = dist;
                %                 obj.eigvec_angle_comparison_mat_ana2asy_permuted = angles;
                %                 obj.eigvec_dot_comparison_mat_ana2asy_permuted = absdot;
                EngineClass.save_var(absdot,obj.resultsPath,'obj_properties','eigvec_dot_comparison_mat_ana2asy_permuted');
            end
            if ~isequal(real(angles),angles)
                disp('Warning: eigvec_angle_comparison_mat has complex values');
            end
        end
        function obj = set_eigvec_comparison_mats2(obj)
            vana = obj.eigenvectors_ana;
            vasy = obj.eigenvectors_asy_permuted;
            [~, dotsana, ~] = obj.compute_comparison_matrix(vana,vana,false,'');
            [~, dotsasy, ~] = obj.compute_comparison_matrix(vasy,vasy,false,'');
            obj.eigvec_dot_comparison_mat_ana2ana = dotsana;
            obj.eigvec_dot_comparison_mat_asy_permuted2asy_permuted = dotsasy;
            vasy1 = obj.eigenvectors_asy;
            [~, dotsasy1, ~] = obj.compute_comparison_matrix(vasy1,vasy1,false,'');
            obj.eigvec_dot_comparison_mat_asy2asy = dotsasy1;
        end
        function obj = set_permutation_eigvec_ana2asy(obj,thresh)
            %             NN = obj.eigvec_angle_comparison_mat_ana2asy;
            NN = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigvec_dot_comparison_mat_ana2asy'));
            %             NN = abs(90 - NN);
            NN = abs(0 - NN);
            [M,I] = max(NN,[],2);
            m = M < thresh;
            disp('set_permutation_eigvec_ana2asy: ');
            disp(['sum(m) = ' num2str(sum(m))]);
            ii = 1:size(I,1);
            I(m) = ii(m);
            [val,ind] = unique(I);
            if isequal(val,ii)
                disp('I is a permutation')
            end
            obj.permutation_eigvec_ana2asy = I';
        end
        function obj = set_eig_asy_permuted(obj)
            obj.eigenvalues_asy_permuted = obj.eigenvalues_asy(obj.permutation_eigvec_ana2asy);
            eig_asy = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvectors_asy'));
            eigenvectors_asy_permuted = eig_asy(:,obj.permutation_eigvec_ana2asy);
            General.save_var(eigenvectors_asy_permuted,fullfile(obj.resultsPath,'obj_properties'),'eigenvectors_asy_permuted');
        end
        function find_eigvec_power_law(obj)
            percentiles = [100, 50, 10, 5, 1];
            propspath = fullfile(obj.resultsPath,'obj_properties');
            eigvecs_ana_ordered_nodes = General.load_var(fullfile(propspath,'eigenvectors_ana_ordered_nodes'));
            k = General.load_var(fullfile('networks',obj.networkName,'degree_vector'));
            [k_sorted,~] = sort(k);
            log_k = log10(k_sorted);
            num_vecs = size(eigvecs_ana_ordered_nodes,2);
            for percentile = percentiles
                disp(['percentile = ' num2str(percentile)]);
                num_nodes = round(num_vecs*percentile/100);
                ps = zeros(num_vecs,2);yfits = zeros(num_nodes,num_vecs);
                rsqs = zeros(num_vecs,1);
                for i1=1:num_vecs
                    log_cur_vec = log10(abs(eigvecs_ana_ordered_nodes(end-num_nodes+1:end,i1)));
                    [p_0, p_1, yfit, rsq] = General.lin_reg(log_k(1:num_nodes),log_cur_vec);
                    ps(i1,:) = [p_0,p_1];
                    yfits(:,i1) = yfit;
                    rsqs(i1,1) = rsq;
                end
                propspath = fullfile(obj.resultsPath,'obj_properties','eigvec_power_law',[num2str(percentile) 'percent']);
                General.save_var(ps,propspath,'ps');
                %             General.save_var(yfits,propspath,'yfits');
                General.save_var(rsqs,propspath,'rsqs');
            end
        end
        function obj = set_weighted_dot_products(obj)
            eigvec_ana = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvectors_ana'));
            eigvec_asy_permuted = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvectors_asy_permuted'));
            s = size(eigvec_ana,2);
            step = s/100; disp('1');
            [wdp,wdpm] = EngineClass.compute_weighted_dot_products(eigvec_ana,eigvec_asy_permuted,obj.degree_vector_weighted,obj.kbinned_mins,true,[0:step:s]',true);
            EngineClass.save_var(wdp,obj.resultsPath,'obj_properties','eigvec_ana_asy_perm_weighted_dot_products_1');
            EngineClass.save_var(wdpm,obj.resultsPath,'obj_properties','eigvec_ana_asy_perm_weighted_dot_products_mean_1');
            disp('2');
            [wdp,wdpm] = EngineClass.compute_weighted_dot_products(eigvec_ana,eigvec_asy_permuted,obj.degree_vector_weighted,obj.kbinned_maxs,false,[0:step:s]',true);
            EngineClass.save_var(wdp,obj.resultsPath,'obj_properties','eigvec_ana_asy_perm_weighted_dot_products_2');
            EngineClass.save_var(wdpm,obj.resultsPath,'obj_properties','eigvec_ana_asy_perm_weighted_dot_products_mean_2');
            disp('3');
            [ind_bins_var,binned_vals] = EngineClass.set_bins_percentiles(20,obj.degree_vector_weighted);
            [wdp,wdpm] = EngineClass.compute_weighted_dot_products(eigvec_ana,eigvec_asy_permuted,obj.degree_vector_weighted,binned_vals(1,:)',true,[0:step:s]',true);
            EngineClass.save_var(wdp,obj.resultsPath,'obj_properties','eigvec_ana_asy_perm_weighted_dot_products_3');
            EngineClass.save_var(wdpm,obj.resultsPath,'obj_properties','eigvec_ana_asy_perm_weighted_dot_products_mean_3');
            disp('4');
            [wdp,wdpm] = EngineClass.compute_weighted_dot_products(eigvec_ana,eigvec_asy_permuted,obj.degree_vector_weighted,binned_vals(1,:)',true,[0:step:s]',false);
            EngineClass.save_var(wdp,obj.resultsPath,'obj_properties','eigvec_ana_asy_perm_weighted_dot_products_4');
            EngineClass.save_var(wdpm,obj.resultsPath,'obj_properties','eigvec_ana_asy_perm_weighted_dot_products_mean_4');
            
        end
        function obj = print_output(obj) %create output file
            T = array2table([obj.solution_t,obj.solution_x],'VariableNames',obj.header);
            if ~isfolder(obj.resultsPath)
                mkdir(obj.resultsPath)
            end
            %writetable(T,fullfile(obj.resultsPath,'output.csv'));
            save(fullfile(obj.resultsPath,'output.mat'),'T');
        end
        
        function obj = calculate_degree(obj)
            B = (obj.adjacencyMatrix>0);
            obj.degree_vector = sum(B,2);
        end
        function obj = calculate_degree_weighted(obj)
            obj.degree_vector_weighted = sum(obj.adjacencyMatrix,2);
            disp('calculate_degree_weighted(obj): obj.degree_vector_weighted = ');
            disp(obj.degree_vector_weighted);
        end
        function obj = load_solution(obj, folderNameStr, solPrefixStr, solSufStr, isIndexed, tTargetVarStr,xTargetVarStr,isDilute)
            if isDilute
                step = 100;
            else
                step = 1;
            end
            tFileNameStr = [solPrefixStr 't' solSufStr];
            xFileNameStr = [solPrefixStr 'x' solSufStr];
            if ~isIndexed
                mypath_t = fullfile(obj.resultsPath,'obj_properties',folderNameStr,tFileNameStr);
                mypath_x = fullfile(obj.resultsPath,'obj_properties',folderNameStr,xFileNameStr);
                mydata = load(mypath_t); sol_t = mydata.var;
                eval([tTargetVarStr ' = sol_t;']);
                mydata = load(mypath_x); sol_x = mydata.var;
                eval([xTargetVarStr ' = sol_x;']);
            else
                break1 = false;break2 = false;
                ind = 1;
                eval([tTargetVarStr ' = []']);eval([xTargetVarStr ' = []']);
                while ~break2 || ~break1
                    if ~break1
                        try
                            mypath_t = fullfile(obj.resultsPath,'obj_properties',folderNameStr,[tFileNameStr '_' num2str(ind)]);
                            mydata = load(mypath_t); sol_t = mydata.var;
                            eval([tTargetVarStr '=  [' tTargetVarStr '; sol_t(1:' num2str(step) ':end,:)];']);clear mydata ;
                        catch exception
                            break1 = true;
                        end
                    end
                    if ~break2
                        try
                            mypath_x = fullfile(obj.resultsPath,'obj_properties',folderNameStr,[xFileNameStr '_' num2str(ind)]);disp(mypath_x);
                            mydata = load(mypath_x); sol_x = mydata.var;
                            eval([xTargetVarStr '=  [' xTargetVarStr '; sol_x(1:' num2str(step) ':end,:)];']);
                        catch exception
                            break2 = true;
                        end
                    end
                    ind = ind + 1;
                end
            end
        end
        function write_gephi_nodes_tables_2(obj)
            sol_x_rand_pert = General.load_var(fullfile(obj.resultsPath,'obj_properties','solution_random_perts_2','sol_x_random_perturbations_1'));
            pert_x_rand_pert = abs((sol_x_rand_pert' - obj.steady_state')');
            sol_t_rand_pert = General.load_var(fullfile(obj.resultsPath,'obj_properties','solution_random_perts_2','sol_t_random_perturbations_1'));
            sol_x_hub_pert = General.load_var(fullfile(obj.resultsPath,'obj_properties','single_node_pert_sol','sol_x_1'));
            pert_x_hub_pert = abs((sol_x_hub_pert' - obj.steady_state')');
            sol_t_hub_pert = General.load_var(fullfile(obj.resultsPath,'obj_properties','single_node_pert_sol','sol_t_1'));
            pert_x_rand_pert_max_t0 = max(pert_x_rand_pert(1,:));
            pert_x_rand_pert(:,obj.N+1) = pert_x_rand_pert_max_t0;
            pert_x_rand_pert_min_t0 = min(pert_x_rand_pert(1,:));
            pert_x_rand_pert(:,obj.N+2) = pert_x_rand_pert_min_t0;
            pert_x_hub_pert_max_t0 = max(pert_x_hub_pert(1,:));
            pert_x_hub_pert(:,obj.N+1) = pert_x_hub_pert_max_t0;
            pert_x_hub_pert_min_t0 = min(pert_x_hub_pert(1,:));
            pert_x_hub_pert(:,obj.N+2) = pert_x_hub_pert_min_t0;
            sys_half_life_2 = obj.sys_half_life_amp{2};
            t0 = 0; t1 = sys_half_life_2; t2 = 2*sys_half_life_2;
            ts = [t0 t1 t2];
            tau = 0;
            for t = ts
                EngineClass.create_gephi_nodes_table(pert_x_rand_pert,sol_t_rand_pert,t,...
                    fullfile(obj.resultsPath,'obj_properties','solution_random_perts_2'),['gephi_nodes_table_sol_x_random_perturbations_1_tau2_' num2str(tau)]);
                EngineClass.create_gephi_nodes_table(pert_x_hub_pert,sol_t_hub_pert,t,...
                    fullfile(obj.resultsPath,'obj_properties','single_node_pert_sol'),['gephi_nodes_table_sol_x_1_tau2_' num2str(tau)]);
                tau = tau+1;
            end
        end
        function obj = write_gephi_nodes_table_jacobian_PEV(obj)
            J = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvectors_ana'));
            nodes = 1:obj.N; nodes=nodes';
            k = obj.degree_vector_weighted;
            eig1 = J(:,1);
            values = abs(eig1)./sum(abs(eig1));
            filepath = fullfile(obj.resultsPath,'obj_properties');filename = 'gephi_nodes_jacobian_PEV';
            General.degree_radius_gephi_nodes_table(nodes,values,k,filepath,filename)
            %             General.general_gephi_nodes_table(nodes,values,x,y,filepath,filename);
        end
        function obj = calc_Q_distribution(obj)
            disp('Begin calc Q distributions');
%            rand_binary_perts = General.load_var(fullfile('networks',obj.networkName,'rand_binary_perts'));
%            cur_perts = rand_binary_perts;
%            filename = 'Q10s_rand_binary_perts';
           filenames = {};
           [~,single_node_combs] = sort(obj.degree_vector_weighted,'descend');
           filenames{end+1} = 'Q10s_single_node_perts';
           
           double_node_combs = General.load_var(fullfile('networks/',obj.networkName,'node_combs','doubles'));
%            cur_perts = zeros(obj.N,size(double_node_combs,1));
%            for i=1:size(cur_perts,2)
%                cur_perts(double_node_combs(i,:),i) = 1;
%            end
           filenames{end+1} = 'Q10s_double_node_perts';
%            cur_perts = double_node_combs';
           
           triple_node_combs = General.load_var(fullfile('networks/',obj.networkName,'node_combs/','triples.mat'));
           filenames{end+1} = 'Q10s_triple_node_perts';
           
           node_combs = {single_node_combs,double_node_combs,triple_node_combs};
%            num_perts = size(cur_perts,2);
           
           t = 0;sys_half_life_amp = obj.sys_half_life_amp;
           mu = obj.mu; k = obj.degree_vector_weighted; xi = obj.xi;ss = obj.steady_state;
           for i2 = 1:length(filenames)
               cur_perts = node_combs{i2}';
               num_perts = size(cur_perts,2);
               Qs = zeros(num_perts,1);
               filename = filenames{i2};
               parfor i1 = 1:num_perts
                   if(mod(i1,10000) ==0)
                       disp(['i1 = ' num2str(i1)])
                   end
                   pert = zeros(1,obj.N);
                   pert(cur_perts(:,i1)) = 1;
                   %                pert = cur_perts(:,i1)';
                   x = ss + pert;
                   [~, ~, ~,~,~,Q,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10] = obj.calc_norms_thetas(x,ss,t,sys_half_life_amp,mu,k,xi);
                   Qs(i1) = Q10;
               end
               disp('Finished calculating Qs. Saving...');
               General.save_var(Qs,fullfile(obj.resultsPath,'obj_properties/'),filename);
               disp('Finished saving.');
           end
        end
        % figures / plots
        function obj = plot_results(obj, isDilute)
            if isDilute
                step = 100;
                suffix = '-diluted';
            else
                step = 1;
                suffix = [];
            end
            if size(obj.solution_x,1) == 0
                if isfile(fullfile(obj.resultsPath,'obj_properties','solution','solution_x_1.mat'))
                    obj.load_solution('solution', 'solution_', '', true, 'obj.solution_t','obj.solution_x',isDilute);
                else
                    obj.load_solution('solution','sol_', '', false, 'obj.solution_t','obj.solution_x',isDilute);
                end
                step = 1;
            end
            if size(obj.steady_state,1) == 0
                obj.steady_state = General.load_var(fullfile(obj.resultsPath,'obj_properites','steady_state'));
            end
            ss = obj.steady_state';
            
            t1 = obj.solution_t;ind_t1 = 1:step:length(t1);t1 = t1(ind_t1);
            x1 = obj.solution_x(ind_t1,:); p1 = (x1'-ss)';
            t2 = [];x2=[];p2=[];
            if ~isempty(obj.solution_t_perturbations)
                t2 = obj.solution_t_perturbations{1};ind_t2 = 1:step:length(t2);
                t2 = t2(ind_t2);x2 = obj.solution_x_perturbations{1};x2 = x2(ind_t2,:);
                p2 = (x2'-ss)';
            end
            t = {t1,t1,t2,t2};x = {x1,p1,x2,p2};
            
            name1 = ['fig1a-1' suffix];name2 = ['fig1a-2' suffix];
            name3 = ['fig1a-3' suffix];name4 = ['fig1a-4' suffix];
            name = {name1,name2,name3,name4};
            ylab1 = 'x'; ylab2 = 'x-ss';ylabs = {ylab1,ylab2,ylab1,ylab2};
            figdesc1 = '$x_0 = rand$';figdesc2 = '$x_0 = SS + \varepsilon*rand$';
            figdesc = {figdesc1,figdesc1,figdesc2,figdesc2};
            for figInd = 1:size(name,2)
                f = figure('Name',name{figInd},'NumberTitle','off');
                plot(t{figInd},x{figInd},'.-','MarkerSize',12);
                legend(obj.header(2:end));
                title({[name{figInd} ' ' obj.scenarioName];obj.desc;figdesc{figInd}},...
                    'Interpreter','latex');
                xlabel(obj.header{1});
                ylabel(ylabs{figInd});
                obj.save_fig(f,name{figInd});
            end
            
            %% fig1b*
            if ~isempty(obj.solution_t_eigvecana)
                epsset = obj.eps;
                numeps = size(epsset,2);
                epsStr = cell(1,numeps);
                ana_asy_t = {obj.solution_t_eigvecana,obj.solution_t_eigvecasy};
                ana_asy_x = {obj.solution_x_eigvecana,obj.solution_x_eigvecasy};
                ana_asyStr = {'Analytical','Asymptotic'};
                n = 5;
                for epsInd = 1:numeps
                    if obj.isInitsLegit(epsInd)
                        epsStr = num2str(epsset(epsInd));
                        for i = 1:n
                            for anaasyInd = 1:2
                                ind_x = 1:step:size(ana_asy_x{anaasyInd}{2,i,epsInd},2);
                                ind_t = 1:step:length(ana_asy_t{anaasyInd}{2,i,epsInd});
                                t = ana_asy_t{anaasyInd}{2,i,epsInd}(ind_t);
                                x = ana_asy_x{anaasyInd}{2,i,epsInd}(ind_t,ind_x);
                                str = ana_asyStr{anaasyInd}(1:3);
                                name = ['fig1b' suffix ' eps' epsStr 'vec_{' num2str(i) ',' str '}'];
                                fname = ['fig1b eps ' strrep(epsStr,'.','p') ' vec' num2str(i) ' ' str];
                                f = figure('Name',name,'NumberTitle','off');
                                plot(t,x,'.-','MarkerSize',12);
                                legend(obj.header(2:step:end));
                                title({[name ' ' obj.scenarioName];obj.desc});
                                xlabel(obj.header{1});
                                ylabel('x');
                                obj.save_fig(f,fname);
                            end
                        end
                    end
                end
            end
        end
        function obj = plot_results2(obj)
            n=5;
            CM = jet(n);
            step = 1;
            ss = obj.steady_state';
            ssn = ss/norm(ss);
            figdesc{1,1} = 'State Angle from Steady State';
            figdesc{2,1} = 'Normalized Dot Product of State and Steady State';
            figdesc{3,1} = 'Absolute Value of Dot Product of State Normal and Steady State';
            figdesc{4,1} = 'State Norm';
            figdesc{5,1} = 'Absolute Value of Dot Product of State and Steady State';
            figdesc{1,2} = 'Perturbation Angle from Steady State';
            figdesc{2,2} = 'Normalized Dot Product of Perturbation and Steady State';
            figdesc{3,2} = 'Absolute Value of Dot Product of Perturbation Normal and Steady State';
            figdesc{4,2} = 'Perturbation Norm';
            figdesc{5,2} = 'Absolute Value of Dot Product of Perturbation and Steady State';
            figdesc{1,3} = 'Perturbation Angle from Eigenvector';
            figdesc{2,3} = 'Normalized Dot Product of Perturbation and Eigenvector';
            figdesc{3,3} = 'Absolute Value of Dot Product of Perturbation Normal and Eigenvector';
            figdesc{4,3} = '';
            figdesc{5,3} = 'Absolute Value of Dot Product of Perturbation and Eigenvector';
            ylabelStr{1,1} = '$\theta_{x,ss}$';
            ylabelStr{2,1} = '$\cos{\theta_{x,ss}}/\cos{\theta_{x_0,ss}}$';
            ylabelStr{3,1} = '$\left|\widehat{x}\cdot \widehat{ss}\right|$';
            ylabelStr{4,1} = '$\left|x\right|$';
            ylabelStr{5,1} = '$\left|x\cdot \widehat{ss}\right|$';
            ylabelStr{1,2} = '$\theta_{x-ss,ss}$';
            ylabelStr{2,2} = '$\cos{\theta_{x-ss,ss}}/\cos{\theta_{x_0-ss,ss}}$';
            ylabelStr{3,2} = '$\left|\widehat{x-ss}\cdot\widehat{ss}\right|$';
            ylabelStr{4,2} = '$\left|x-ss\right|$';
            ylabelStr{5,2} = '$\left|x-ss\cdot\widehat{ss}\right|$';
            ylabelStr{1,3} = '$\theta_{x-ss,v_i}$';
            ylabelStr{2,3} = '$\cos{\theta_{x-ss,v_i}}/\cos{\theta_{x_0-ss,v_i}}$';
            ylabelStr{3,3} = '$\left|\widehat{x-ss}\cdot v_i\right|$';
            ylabelStr{4,3} = '';
            ylabelStr{5,3} = '$\left|x-ss \cdot v_i\right|$';
            for epsInd = 1:length(obj.eps)
                epsStr = num2str(obj.eps(epsInd));
                angles_ana = cell(5,n,3);
                angles_asy = cell(5,n,3);
                for i = 1:n
                    disp(['vec ' num2str(i)])
                    if size(obj.isInitsLegit,2)>=i && obj.isInitsLegit(epsInd,i) && size(obj.solution_x_eigvecana,2)>=i
                        x1 = obj.solution_x_eigvecana{1,i,epsInd};
                        n1 = size(x1,1);
                        for j = 1:n1
                            if mod(j,1000) == 0
                                disp(['row ' num2str(j)]);
                            end
                            v1 = x1(j,:); % current state
                            normv1 = norm(v1);
                            u1 = v1' - ss; % current deviation from steady state
                            normu1 = norm(u1);
                            v1n = v1/normv1;
                            u1n = u1/normu1;
                            vec = obj.eigenvectors_ana(:,i);
                            dot_u1n_ssn = dot(ssn,u1n);
                            dot_v1n_ssn = dot(ssn,v1n);
                            dot_u1n_vec = dot(vec,u1n);
                            dot_u1_ssn = dot(ssn,u1);
                            dot_v1_ssn = dot(ssn,v1);
                            dot_u1_vec = dot(u1,vec);
                            if j==1
                                dot_1 = {dot_v1n_ssn,dot_u1n_ssn,dot_u1n_vec};
                            end
                            angles_ana{1,i,1}(j,1) = rad2deg(real(acos(dot_v1n_ssn)));
                            angles_ana{2,i,1}(j,1) = dot_v1n_ssn/dot_1{1};
                            angles_ana{3,i,1}(j,1) = abs(dot_v1n_ssn);
                            angles_ana{4,i,1}(j,1) = normv1;
                            angles_ana{5,i,1}(j,1) = abs(dot_v1_ssn);
                            angles_ana{1,i,2}(j,1) = rad2deg(real(acos(dot_u1n_ssn)));
                            angles_ana{2,i,2}(j,1) = dot_u1n_ssn/dot_1{2};
                            angles_ana{3,i,2}(j,1) = abs(dot_u1n_ssn);
                            angles_ana{4,i,2}(j,1) = normu1;
                            angles_ana{5,i,2}(j,1) = abs(dot_u1_ssn);
                            angles_ana{1,i,3}(j,1) = rad2deg(real(acos(dot_u1n_vec)));
                            angles_ana{2,i,3}(j,1) = dot_u1n_vec/dot_1{3};
                            angles_ana{3,i,3}(j,1) = abs(dot_u1n_vec);
                            angles_ana{4,i,3}(j,1) = 0;
                            angles_ana{5,i,3}(j,1) = abs(dot_u1_vec);
                        end
                    end
                    if size(obj.isInitsLegit,2)>=n+i && obj.isInitsLegit(epsInd,n+i)
                        x2 = obj.solution_x_eigvecasy{1,i,epsInd};
                        n2 = size(x2,1);
                        for j = 1:n2
                            v2 = x2(j,:);
                            normv2 = norm(v2);
                            u2 = v2' - ss;
                            normu2 = norm(u2);
                            v2n = v2/normv2;
                            u2n = u2/normu2;
                            vec = obj.eigenvectors_asy_permuted(:,i);
                            dot_u2n_ssn = dot(ssn,u2n);
                            dot_v2n_ssn = dot(ssn,v2n);
                            dot_u2n_vec = dot(vec,u2n);
                            dot_u2_ssn = dot(ssn,u2);
                            dot_v2_ssn = dot(ssn,v2);
                            dot_u2_vec = dot(vec,u2);
                            if j==1
                                dot_1 = {dot_v2n_ssn,dot_u2n_ssn,dot_u2n_vec};
                            end
                            angles_asy{1,i,1}(j,1) = rad2deg(real(acos(dot_v2_ssn)));
                            angles_asy{2,i,1}(j,1) = dot_v2n_ssn/dot_1{1};
                            angles_asy{3,i,1}(j,1) = abs(dot_v2n_ssn);
                            angles_asy{4,i,1}(j,1) = normv2;
                            angles_asy{5,i,1}(j,1) = abs(dot_v2_ssn);
                            angles_asy{1,i,2}(j,1) = rad2deg(real(acos(dot_u2n_ssn)));
                            angles_asy{2,i,2}(j,1) = dot_u2n_ssn/dot_1{2};
                            angles_asy{3,i,2}(j,1) = abs(dot_u2n_ssn);
                            angles_asy{4,i,2}(j,1) = normu2;
                            angles_asy{5,i,2}(j,1) = abs(dot_u2_ssn);
                            angles_asy{1,i,3}(j,1) = rad2deg(real(acos(dot_u2n_vec)));
                            angles_asy{2,i,3}(j,1) = dot_u2n_vec/dot_1{3};
                            angles_asy{3,i,3}(j,1) = abs(dot_u2n_vec);
                            angles_asy{4,i,3}(j,1) = 0;
                            angles_asy{5,i,3}(j,1) = abs(dot_u2_vec);
                        end
                    end
                end
                
                ssstr = 'ss + ';
                letters = ['a','b','c'];
                name = cell(5,3);
                fname = cell(5,3);
                for c = 1:5
                    for d=1:3
                        name{c,d} = ['fig2' letters(d) '-' num2str(c) ' eps ' epsStr];
                        fname{c,d} = ['fig2' letters(d) '-' num2str(c) ' eps ' strrep(epsStr,'.','p')];
                    end
                end
                
                numfigs = size(name,2);
                myPlot = {@plot,@semilogy};
                for figInd1 = 1:size(name,1)
                    for figInd = 1:numfigs
                        if ~isempty(figdesc{figInd1,figInd}) %&& ~isempty(obj.solution_t_eigvecana{1,i,epsInd})
                            myPlotInd = 1;
                            while true
                                namestr = name{figInd1,figInd};
                                fnamestr = fname{figInd1,figInd};
                                if myPlotInd == 2
                                    namestr = [namestr '-log'];
                                    fnamestr = [fnamestr '-log'];
                                end
                                f = figure('Name',namestr,'NumberTitle','off');
                                legendStr = {};
                                for i = 1:n
                                    if size(obj.isInitsLegit,2) >=i && obj.isInitsLegit(epsInd,i) && size(obj.solution_t_eigvecana,2)>=i && ~isempty(obj.solution_t_eigvecana{1,i,epsInd})
                                        epsStr = num2str(obj.eps_adjusted(epsInd,i));
                                        myPlot{myPlotInd}(obj.solution_t_eigvecana{1,i,epsInd}(1:step:end,1),angles_ana{figInd1,i,figInd}(1:step:end,1),'.','Color',CM(i,:),'MarkerSize',12);
                                        hold on;
                                        legendStr{end+1} = ['$x_0 = ' ssstr epsStr ' * v_{' num2str(i) ',ana}$'];
                                    end
                                end
                                for i = 1:n
                                    if size(obj.isInitsLegit,2)>=n+i && obj.isInitsLegit(n+i) && ~isempty(obj.solution_t_eigvecasy{1,i,epsInd})
                                        epsStr = num2str(obj.eps_adjusted(epsInd,n+i));
                                        myPlot{myPlotInd}(obj.solution_t_eigvecasy{1,i,epsInd}(1:step:end,1),angles_asy{figInd1,i,figInd}(1:step:end,1),'^','Color',CM(i,:));
                                        hold on;
                                        legendStr{end+1} = ['$x_0 = ' ssstr epsStr ' * v_{' num2str(i) ',asy}$'];
                                    end
                                end
                                title({[namestr ' ' obj.scenarioName];figdesc{figInd1,figInd}});
                                xlabel('Time');
                                ylabel(ylabelStr{figInd1,figInd},'Interpreter','latex','FontSize',14);
                                legend(legendStr,'Interpreter','latex','FontSize',14);
                                obj.save_fig(f,fnamestr);
                                if (figInd ~= 1 && (figInd1 == 4 || figInd1 ==5)) && (myPlotInd == 1)
                                    myPlotInd = 2;
                                else
                                    break
                                end
                            end
                        end
                    end
                end
                
            end
            %% fig2e fig2f, fig2e-polar, fig2f-polar
            t = obj.solution_t;
            x = obj.solution_x;
            xnorm = vecnorm(x');
            m = size(t,1);
            angles = cell(3,1);
            angles{1,1} = zeros(m,2);
            angles{2,1} = zeros(m,2);
            angles{3,1} = zeros(m,2);
            pert = (x' - ss)';
            pertnorm = vecnorm(pert');
            pertnorm_no0 = pertnorm(pertnorm>0);
            min_pertnorm = min(pertnorm_no0);
            pertnorm_normalized = pertnorm_no0/min_pertnorm;
            logpertnorm = log10(pertnorm_normalized);
            for i=1:m
                v = x(i,:)';
                vn = v/norm(v);
                u = pert(i,:);
                un = u/norm(u);
                dot_vn_ssn = dot(ssn,vn);
                dot_un_ssn = dot(ssn,un);
                if i==1
                    dot_1 = {dot_vn_ssn,dot_un_ssn};
                end
                angles{1,1}(i,1) = rad2deg(real(acos(dot_vn_ssn)));
                angles{1,1}(i,2) = rad2deg(real(acos(dot_un_ssn)));
                angles{2,1}(i,1) = dot_vn_ssn/dot_1{1};
                angles{2,1}(i,2) = dot_un_ssn/dot_1{2};
                angles{3,1}(i,1) = abs(dot_vn_ssn);
                angles{3,1}(i,2) = abs(dot_un_ssn);
            end
            name = {'fig2e-1','fig2f-1';'fig2e-2','fig2f-2';'fig2e-3','fig2f-3'};
            for figInd1 = 1:3
                for figInd = 1:2
                    f = figure('Name',name{figInd1,figInd},'NumberTitle','off');
                    hold on;
                    plot(t,angles{figInd1,1}(:,figInd),'.b','MarkerSize',12);
                    title({[name{figInd1,figInd} ' ' obj.scenarioName];figdesc{figInd1,figInd}});
                    xlabel('Time');
                    ylabel(ylabelStr{figInd1,figInd},'Interpreter','latex','FontSize',14);
                    legendStr = 'x_0 = rand';
                    legend(legendStr,'FontSize',14);
                    obj.save_fig(f,name{figInd1,figInd});
                end
            end
            name = {'fig2e-polar','fig2f-polar'};
            for figInd = 1:2
                f = figure('Name',name{1,figInd},'NumberTitle','off');
                %hold on;
                ax = polaraxes;
                ax.RAxis.Label.String = 'Log(Pert norm)';
                hold on;
                polarplot(deg2rad(angles{1,1}(pertnorm>0,figInd)),logpertnorm,'.b','MarkerSize',12);
                title({[name{1,figInd} ' ' obj.scenarioName];figdesc{1,figInd}});
                legendStr = 'x_0 = rand';
                legend(legendStr,'FontSize',14);
                obj.save_fig(f,name{1,figInd});
            end
            %% fig2g, fig2h, fig2i, fig2j, fig2k
            perts = {};
            ts = {};
            ms = {};
            perts{1} = pert;
            perts{2} = (obj.solution_x_perturbations{1}' - ss)';
            ts{1} = t;
            ts{2} = obj.solution_t_perturbations{1};
            epsInd = 1;
            ms{1} = m;
            ms{2} = size(ts{2},1);
            
            for i = 1:size(obj.solution_x_eigvec,2)
                if ~isempty(obj.solution_x_eigvec{1,i,1})
                    x1 = (obj.solution_x_eigvec{1,i,1}'- ss)';
                    perts{end+1} = x1;
                    t1 = obj.solution_t_eigvec{1,i,1};
                    ts{end+1} = t1;
                    ms{end+1} = size(t1,1);
                end
            end
            for i=1:n
                if obj.isInitsLegit(1,i)
                    perts{end+1} = (obj.solution_x_eigvecasy{1,i,1}' - ss)';
                    ts{end+1} = obj.solution_t_eigvecasy{1,i,1};
                    ms{end+1} = size(ts{end},1);
                end
            end
            
            
            numfigs = 4;a = ['g','h','i','j','l','m'];str = num2str(obj.perturbation_factor);
            figdesc1 = {'','Normalized dot product of ','Absolute value of dot product of ',...
                'Absolute value of dot product of '};
            figdesc2 = {'Random perturbation ',['Random ($<$' str ' of SS) perturbation '],...
                'Perturbation ($\Sigma_{i=1}^5 v_{i,ana}$) ',...
                'Perturbation ($\Sigma_{i=1}^5 v_{i,asy}$) '...
                'Perturbation ($x \cdot eigvec_{i,ana} = C$) '...
                'Perturbation ($x \cdot eigvec_{i,asy} = C$) '
                };
            figdesc3 = {'angle from eigenvector','and eigenvector'...
                'normal and eigenvector'...
                'and eigenvector'};
            ylabelStr = {'$\theta_{x-ss,Eigvec}$';...
                '$\cos{\theta_{x-ss,eigvec}}/\cos{\theta_{x_0-ss,eigvec}}$';...
                '$\left|\widehat{x-ss} \cdot eigvec \right|$';...
                '$\left|x-ss \cdot eigvec \right|$'};
            for r = 1:size(perts,2)
                angles = cell(numfigs,2); % numfigs types of graphs, 2 types of eigenvectors (ana, asy)
                for a1 = 1:size(angles,1)
                    for a2 = 1:size(angles,2)
                        angles{a1,a2} = zeros(ms{r},n);
                    end
                end
                %                 angles{1,1} = zeros(ms{r},n); % angles with analytical eigenvectors
                %                 angles{2,1} = zeros(ms{r},n); % normalized dot product with analytical eigenvectors
                %                 angles{3,1} = zeros(ms{r},n); % abs(dot(pert,eig_ana))
                %                 angles{1,2} = zeros(ms{r},n); % angles with asymptotic eigenvectors
                %                 angles{2,2} = zeros(ms{r},n); % normalized dot product with asymptotic eigenvectors
                %                 angles{3,2} = zeros(ms{r},n); % abs(dot(pert,eig_asy))
                
                for v = 1:n
                    vec_ana = obj.eigenvectors_ana(:,v);
                    vec_asy = obj.eigenvectors_asy_permuted(:,v);
                    for i = 1:ms{r}
                        u = perts{r}(i,:);
                        un = u/norm(u);
                        dot_u_vec_ana = dot(u,vec_ana);
                        dot_u_vec_asy = dot(u,vec_asy);
                        dot_un_vec_ana = dot(un,vec_ana);
                        dot_un_vec_asy = dot(un,vec_asy);
                        if i==1
                            dot_1 = {dot_un_vec_ana, dot_un_vec_asy};
                        end
                        angles{1,1}(i,v) = rad2deg(real(acos(dot_un_vec_ana)));
                        angles{2,1}(i,v) = dot_un_vec_ana/dot_1{1};
                        angles{3,1}(i,v) = abs(dot_un_vec_ana);
                        angles{4,1}(i,v) = abs(dot_u_vec_ana);
                        angles{1,2}(i,v) = rad2deg(real(acos(dot_un_vec_asy)));
                        angles{2,2}(i,v) = dot_un_vec_asy/dot_1{2};
                        angles{3,2}(i,v) = abs(dot_un_vec_asy);
                        angles{4,2}(i,v) = abs(dot_u_vec_asy);
                    end
                end
                
                for figInd = 1:numfigs
                    myPlotInd = 1;
                    if r>size(a,2)
                        str1 = ['k-vec' num2str(r-size(a,2))];
                        vecind = r-size(a,2);
                        str2 = ['Perturbation ' num2str(obj.eps_adjusted(5 + vecind))...
                            '*$eigvec_{asy,' num2str(vecind) '}$ '];
                    else
                        str1 = a(r);
                        str2 = figdesc2{r};
                    end
                    name = ['fig2' str1 '-' num2str(figInd)];
                    figdesc = [figdesc1{figInd} str2 figdesc3{figInd}];
                    while true
                        if myPlotInd == 2
                            name = [name '-log'];
                        end
                        f = figure('Name',name,'NumberTitle','off');
                        legendStr = {};
                        for v=1:n
                            myPlot{myPlotInd}(ts{r},angles{figInd,1}(:,v),'.','Color',CM(v,:),'MarkerSize',12);
                            legendStr{1,v} = ['Eigvec_{ana,' num2str(v) '}'];
                            hold on;
                        end
                        if ~(r>size(a,2))
                            for v=1:n
                                myPlot{myPlotInd}(ts{r},angles{figInd,2}(:,v),'^','Color',CM(v,:),'MarkerSize',10);
                                legendStr{1,n+v} = ['Eigvec_{asy,' num2str(v) '}'];
                                hold on;
                            end
                        end
                        
                        title({[name ' ' obj.scenarioName];figdesc},'Interpreter','latex');
                        xlabel('Time');
                        ylabel(ylabelStr{figInd,1},'Interpreter','latex','FontSize',14);
                        legend(legendStr,'FontSize',14);
                        obj.save_fig(f,name);
                        if (figInd == 4 || figInd ==3) && myPlotInd == 1
                            myPlotInd = 2;
                        else
                            break
                        end
                    end
                end
            end
        end
        function obj = plot_steady_vs_degree(obj)
            %% fig3-*
            names = '123';
            xvals = {obj.degree_vector,obj.degree_vector_weighted,obj.ki_nn};
            xlabs = {'k','k-weighted','k_{i,nn}'};
            for i1 = 1:length(names)
                name = ['fig3-' names(i1)];
                f = figure('Name',name,'NumberTitle','off');
                plot(xvals{i1},obj.steady_state,'.','MarkerSize',12);
                title({[name ' ' obj.scenarioName];obj.desc});
                xlabel(xlabs{i1});
                ylabel('Steady State');
                obj.save_fig(f,name);
            end
            name = 'fig3a';
            f = figure('Name',name,'NumberTitle','off');
            loglog(obj.degree_vector,obj.steady_state,'.','MarkerSize',12);
            title({[name ' ' obj.scenarioName];obj.desc});
            xlabel('node degree');
            ylabel('Steady State');
            obj.save_fig(f,name);
            name = 'fig3b';
            f = figure('Name',name,'NumberTitle','off');
            plot(obj.degree_vector_weighted,obj.steady_state,'.','MarkerSize',12);
            hold on;
            plot(obj.degree_vector_weighted,obj.steady_state_calculated,'^r');
            title({[name ' ' obj.scenarioName];obj.desc});
            xlabel('k - weighted node degree');
            ylabel('Steady State');
            legend('Numerical Integration','x_i = R^{-1}(q_i)');
            obj.save_fig(f,name);
            % Relevant for SIS
            name = 'fig3c';
            f = figure('Name',name,'NumberTitle','off');
            plot(obj.degree_vector_weighted, -obj.steady_state+1,'.','MarkerSize',12);
            hold on;
            plot(obj.degree_vector_weighted, 1./obj.degree_vector_weighted,'sk','MarkerSize',10);
            plot(obj.degree_vector_weighted, 1./obj.degree_vector_weighted - 1./obj.degree_vector_weighted.^2,'om','MarkerSize',10);
            legend('Numerical Integration', 'Second order approximation k^{-1}','Third order approximation k^{-1} - k^{-2}');
            xlabel('k - weighted node degree');
            ylabel('1 - Steady State');
            title({[name ' ' obj.scenarioName];obj.desc});
            obj.save_fig(f,name);
            % Relevant for SIS
            name = 'fig3d';
            f = figure('Name',name,'NumberTitle','off');
            loglog(obj.degree_vector_weighted, 1 - obj.steady_state,'.','MarkerSize',10);
            hold on;
            loglog(obj.degree_vector_weighted, 1./obj.degree_vector_weighted,'sk','MarkerSize',12);
            loglog(obj.degree_vector_weighted, 1./obj.degree_vector_weighted - 1./obj.degree_vector_weighted.^2,'om','MarkerSize',12);
            legend('Numerical Integration', 'Second order approximation k^{-1}','-1*Third order approximation k^{-2} - k^{-1}');
            xlabel('k - weighted node degree');
            ylabel('1 - Steady State');
            title({[name ' ' obj.scenarioName];obj.desc});
            obj.save_fig(f,name);
            %Relevant for REG
            name = 'fig3e';
            f = figure('Name',name,'NumberTitle','off');
            plot(obj.degree_vector_weighted, obj.steady_state,'.','MarkerSize',12);
            hold on;
            plot(obj.degree_vector_weighted, obj.degree_vector_weighted,'sk','MarkerSize',10);
            legend('Numerical Integration', 'First order approximation k');
            xlabel('k - weighted node degree');
            ylabel('Steady State');
            title({[name ' ' obj.scenarioName];obj.desc});
            obj.save_fig(f,name);
            % Relevant for REG
            name = 'fig3f';
            f = figure('Name',name,'NumberTitle','off');
            loglog(obj.degree_vector_weighted, obj.steady_state,'.','MarkerSize',10);
            hold on;
            loglog(obj.degree_vector_weighted, obj.degree_vector_weighted,'sk','MarkerSize',12);
            legend('Numerical Integration', 'First order approximation k');
            xlabel('log(k) - weighted node degree');
            ylabel('log(Steady State)');
            title({[name ' ' obj.scenarioName];obj.desc});
            obj.save_fig(f,name);
        end
        function obj = plot_jacobian(obj)
            Wij_ana = General.load_var(obj.resultsPath,'obj_properties','Wij_ana');
            Wij_asy = General.load_var(obj.resultsPath,'obj_properties','Wij_asy');
            mydata = load(fullfile(obj.resultsPath,'obj_properties','Wij_asy_v2_set.mat'),'var');Wij_asy_v2_set=mydata.var;
            mydata = load(fullfile(obj.resultsPath,'obj_properties','Dii_asy_v2_set.mat'),'var');Dii_asy_v2_set=mydata.var;
            mydata = load(fullfile(obj.resultsPath,'obj_properties','Wij_asy_v2_binned_set.mat'),'var');Wij_asy_v2_binned_set=mydata.var;
            mydata = load(fullfile(obj.resultsPath,'obj_properties','Dii_asy_v2_binned_set.mat'),'var');Dii_asy_v2_binned_set=mydata.var;
            name = 'fig4a';
            f = figure('Name',name,'NumberTitle','off');
            x = obj.degree_vector_weighted;
            y = obj.Dii_ana;
            loglog(x,y,'.','MarkerSize',12);
            hold on;
            y1 = obj.Dii_asy;
            loglog(x,y1,'^','MarkerSize',10);
            y2 = Dii_asy_v2_set{1,1};
            loglog(x,y2,'v','MarkerSize',10);
            xlabel('k_i - weighted');
            ylabel('J_{ii}');
            title({[name ' ' obj.scenarioName];obj.desc});
            legend('Analytic Jacobian', 'Asymptotic Jacobian', 'Asymptotic Jacobian v2');
            obj.save_fig(f,name);
            %%%%
            name = 'fig4b';
            f = figure('Name',name,'NumberTitle','off');
            x1 = x * x';
            y = Wij_ana;
            loglog(x1(:),y(:),'.','MarkerSize',12);
            hold on;
            y1 = Wij_asy;
            loglog(x1(:),y1(:),'^','MarkerSize',10);
            y2 = Wij_asy_v2_set{1,1};
            loglog(x1(:),y2(:),'v','MarkerSize',10);
            xlabel('k_ik_j - weighted');
            ylabel('W_{ij}');
            title({[name ' ' obj.scenarioName];obj.desc});
            legend('Analytic Jacobian','Asymptotic Jacobian',...
                'Asymptotic Jacobian v2');
            obj.save_fig(f,name);
            %%%%
            name = 'fig4c';
            f = figure('Name',name,'NumberTitle','off');
            x2 = (x.^obj.nu) * (x.^obj.rho)';
            loglog(x2(:),y(:),'.','MarkerSize',12);
            hold on;
            loglog(x2(:),y1(:),'^','MarkerSize',10);
            loglog(x2(:),y2(:),'v','MarkerSize',10);
            xlabel(['k_i^{\nu}k_j^{\rho} - weighted, \nu = ' num2str(obj.nu)...
                ', \rho = ' num2str(obj.rho)]);
            ylabel('W_{ij}');
            title({[name ' ' obj.scenarioName];obj.desc});
            legend('Analytic Jacobian','Asymptotic Jacobian',...
                'Asymptotic Jacobian v2');
            obj.save_fig(f,name);
            %%%%
            name = 'fig4d';
            figdesc = 'Dii Real and Theoretical binned';
            f = figure('Name',name,'NumberTitle','off');
            x = obj.kbinned;
            y1 = obj.Dii_anabinned;
            y2 = obj.Dii_asybinned;
            y3 = Dii_asy_v2_binned_set{1,1};
            loglog(x,y1,'s','MarkerSize',12);
            hold on;
            loglog(x,y2,'^','MarkerSize',10);
            loglog(x,y3,'v','MarkerSize',10);
            xlabel('k_i - weighted');
            ylabel('J_{ii}');
            title({[name ' ' obj.scenarioName];obj.desc;figdesc});
            legend('Analytic Jacobian','Asymptotic Jacobian',...
                'Asymptotic Jacobian v2');
            obj.save_fig(f,name);
            %%%%
            name = 'fig4e';
            figdesc = 'Wij Real and Theoretical binned';
            f = figure('Name',name,'NumberTitle','off');
            x = obj.kijbinned;
            y1 = obj.Wij_anabinned;
            y2 = obj.Wij_asybinned;
            y3 = Wij_asy_v2_binned_set{1,1};
            loglog(x,y1,'s','MarkerSize',12);
            hold on;
            loglog(x,y2,'^','MarkerSize',10);
            loglog(x,y3,'v','MarkerSize',10);
            xlabel(['k_i^{\nu}k_j^{\rho} - weighted, \nu = ' num2str(obj.nu)...
                ', \rho = ' num2str(obj.rho)]);
            ylabel('W_{ij}');
            title({[name ' ' obj.scenarioName];obj.desc;figdesc});
            legend('Analytic Jacobian','Asymptotic Jacobian',...
                'Asymptotic Jacobian v2');
            obj.save_fig(f,name);
        end
        function plot_eigenvalues(obj)
            try
                mydata = load(fullfile(obj.resultsPath,'obj_properties','eigvals_asy_v2_set.mat'),'var');eigvals_asy_v2_set=mydata.var;
            catch exception
                disp(['plot_eigenvalues: load eigvals_asy_v2_set.mat exception: ' exception.message]);
            end
            try
                mydata = load(fullfile(obj.resultsPath,'obj_properties','eigvals_asy_v5.mat'),'var');eigvals_asy_v5=mydata.var;
            catch exception
                disp(['plot_eigenvalues: load eigvals_asy_v5.mat exception: ' exception.message]);
            end
            obj.eigenvalues_asy = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvalues_asy'));
            obj.eigenvalues_ana = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvalues_ana'));
            suf1 = {'-Real', '-Complex','-abs'};
            desc1 = {' real part', ' complex valued', ' absolute value'};
            xlabs = {'$i$', '$Re(\lambda_i)$', '$i$'};
            ylabs = {'$Re(\lambda_i)$','$Im(\lambda_i)$' ,'$\mid \lambda_i \mid$'};
            fs = {@real, @(x) (1*x), @abs};
            for i1 = 1:length(fs)
                name = ['fig5a' suf1{i1}];
                figdesc = ['Jacobian Eigenvalues' desc1{i1}];
                legendStr = {};
                legendInd = 1;
                f = figure('Name',name,'NumberTitle','off');
                plot(fs{i1}(obj.eigenvalues_ana),'*-');
                legendStr{legendInd} = 'Analytic Jacobian'; legendInd = legendInd + 1;
                hold on;
                plot(fs{i1}(obj.eigenvalues_asy),'^-');
                legendStr{legendInd} = 'Asymptotic Jacobian'; legendInd = legendInd + 1;
                if exist('eigvals_asy_v2_set','var')
                    n = size(eigvals_asy_v2_set,1);
                    for i =1:n
                        plot(fs{i1}(eigvals_asy_v2_set{i,1}),'v-');
                        legendStr{legendInd} = ['Asymptotic Jacobian v2-' num2str(i)];legendInd=legendInd+1;
                    end
                end
                plot(fs{i1}(obj.eigvals_v3),'s-');
                legendStr{legendInd} = ['Asymptotic Jacobian v3, C_D = ' num2str(obj.C_D_v3)];legendInd=legendInd+1;
                if exist('eigvals_asy_v5','var')
                    plot(fs{i1}(eigvals_asy_v5),'o-');
                    legendStr{legendInd} = ['Asymptotic Jacobian v5'];legendInd=legendInd+1;
                end
                xlabel(xlabs{i1},'Interpreter','latex');
                ylabel(ylabs{i1},'Interpreter','latex');
                title({[name ' ' obj.scenarioName];obj.desc;figdesc});
                legend(legendStr);
                obj.save_fig(f,name);
            end
        end
        %% fig5b-*
        function plot_eigenvalues2(obj)
            obj.eigenvalues_asy_permuted = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvalues_asy_permuted'));
            obj.eigenvalues_asy = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvalues_asy'));
            obj.eigenvalues_ana = General.load_var(fullfile(obj.resultsPath,'obj_properties','eigenvalues_ana'));
            if obj.mu < 0
                dir = 'descend';
            else
                dir = 'ascend';
            end
            [deg_sorted,~] = sort(obj.degree_vector_weighted,dir);
            name = 'fig5b'; f = figure('Name',name,'NumberTitle','off');
            loglog(deg_sorted,abs(obj.eigenvalues_ana),'-*');hold on;
            loglog(deg_sorted,abs(obj.eigenvalues_asy),'-^');hold on;
            loglog(deg_sorted,abs(obj.eigenvalues_asy_permuted),'-v');hold on;
            loglog(deg_sorted,deg_sorted.^obj.mu,'-s');
            xlabel('$k_i$','Interpreter','latex');ylabel('$\mid\lambda_{k_i}\mid$','Interpreter','latex');
            figdesc = 'Comparison of eigenvalues and $k_i^\mu$';
            title({[name ' ' obj.scenarioName];obj.desc;figdesc},'interpreter','latex');
            legend('$\lambda_{i,Re}$ sorted','$\lambda_{i,Th}$ sorted','$\lambda_{i,Th Permuted}$ sorted','$k_i^\mu$','interpreter','latex');
            obj.save_fig(f,name);
        end
        function plot_eigenvectors(obj,usePermuted)
            %n = obj.numeigenplots2;
            n = obj.numeigen;
            if usePermuted
                e2 = obj.eigenvectors_asy_permuted(:,1:n);
                suf = '-Permuted';
                dist = obj.eigvec_dist_comparison_mat_ana2asy_permuted(1:n,1:n);
                angle = obj.eigvec_angle_comparison_mat_ana2asy_permuted(1:n,1:n);
                vecdot = obj.eigvec_dot_comparison_mat_ana2asy_permuted(1:n,1:n);
            else
                e2 = obj.eigenvectors_asy(:,1:n);
                suf = '';
                dist = obj.eigvec_dist_comparison_mat_ana2asy(1:n,1:n);
                angle = obj.eigvec_angle_comparison_mat_ana2asy(1:n,1:n);
                vecdot = obj.eigvec_dot_comparison_mat_ana2asy(1:n,1:n);
            end
            e1 = obj.eigenvectors_ana(:,1:n);
            % fig6a
            name = ['fig6a' suf];
            figdesc = ['Jacobian Eigenvectors Comparison' suf];
            f = figure('Name',name,'NumberTitle','off');
            vhat = e1 - e2;
            y = vecnorm(vhat);
            plot(y,'*-');
            xlabel('n');
            ylabel('$\mid v-\hat{v} \mid$','interpreter','latex');
            title({[name ' ' obj.scenarioName];obj.desc;figdesc});
            legend('Analytic Jacobian vs Asymptotic Jacobian');
            obj.save_fig(f,name);
            %%% fig6b
            name = ['fig6b' suf];
            figdesc = ['Jacobian Eigenvectors Angle Comparison' suf];
            f = figure('Name',name,'NumberTitle','off');
            d = dot(e1,e2);
            th = acos(d)*180/pi;
            plot(real(th),'*-');
            xlabel('n');
            ylabel('$\theta_{v,\hat{v}} [deg]$','interpreter','latex');
            title({[name ' ' obj.scenarioName];obj.desc;figdesc});
            legend('Analytic Jacobian vs Asymptotic Jacobian');
            obj.save_fig(f,name);
            %%% fig6b-1
            suf2 = {'','-1'};d = dot(e1,e2);
            x1 = 1:size(d,2);x2 = -1./obj.eigenvalues_ana(1:n,:);
            x = {x1,x2};xlabs = {'v_i','$-1/\lambda_i$'};
            figdesc = ['Jacobian Eigenvectors Dot Comparison' suf];
            for i1 = 1:size(x,2)
                name = ['fig6b-1' suf suf2{i1}];
                f = figure('Name',name,'NumberTitle','off');
                plot(x{i1},abs(d),'*-');
                xlabel(xlabs{i1},'Interpreter','latex');
                ylabel('$\mid v \cdot \hat{v} \mid$','interpreter','latex');
                title({[name ' ' obj.scenarioName];obj.desc;figdesc});
                legend('Analytic Jacobian vs Asymptotic Jacobian');
                obj.save_fig(f,name);
            end
            %% fig6c
            name = ['fig6c' suf];
            figdesc = ['Jacobian Eigenvectors Distance Comparison Matrix' suf];
            f = figure('Name',name,'NumberTitle','off');
            image(dist,'CDatamapping','scaled');
            colorbar;
            ylabel('Analytic Eigenvectors');
            xlabel(['Asymptotic Eigenvectors' suf]);
            title({[name ' ' obj.scenarioName];obj.desc;figdesc;'$\mid v-\hat{v} \mid$'},'interpreter','latex');
            obj.save_fig(f,name);
            %% fig6d
            name = ['fig6d' suf];
            figdesc = ['Jacobian Eigenvectors Angle Comparison Matrix' suf];
            f = figure('Name',name,'NumberTitle','off');
            image(angle,'CDatamapping','scaled');
            colorbar;
            ylabel('Analytic Eigenvectors');
            xlabel(['Asymptotic Eigenvectors' suf]);
            title({[name ' ' obj.scenarioName];obj.desc;figdesc;'$\theta_{v,\hat{v}} [deg]$'},'interpreter','latex');
            obj.save_fig(f,name);
            
            %% fig6d-1
            name = ['fig6d-1' suf];
            figdesc = ['Jacobian Eigenvectors Dot Product Comparison Matrix' suf];
            f = figure('Name',name,'NumberTitle','off');
            image(vecdot,'CDatamapping','scaled');
            colorbar;
            ylabel('Analytic Eigenvectors');
            xlabel(['Asymptotic Eigenvectors' suf]);
            title({[name ' ' obj.scenarioName];obj.desc;figdesc;'$\mid v \cdot \hat{v} \mid $'},'interpreter','latex');
            obj.save_fig(f,name);
            %% fig6e
            name = ['fig6e' suf];
            figdesc = ['Jacobian Eigenvectors vs Steady State Angle Comparison' suf];
            f = figure('Name',name,'NumberTitle','off');
            ss = obj.steady_state';
            ss = ss/norm(ss);
            ss = ss .* ones(size(e1));
            dana = dot(ss,e1);
            thetaana = acos(dana)*180/pi;
            dasy = dot(ss,e2);
            thetaasy = acos(dasy)*180/pi;
            hold on;
            plot(real(thetaana),'*-');
            plot(real(thetaasy),'o-');
            xlabel('n');
            ylabel('$\theta_{v,ss} [deg]$','interpreter','latex');
            title({[name ' ' obj.scenarioName];obj.desc;figdesc});
            legend('Analytic Jacobian','Asymptotic Jacobian');
            obj.save_fig(f,name);
            %% fig6f* fig6g* fig6h* matrix image
            n = obj.numeigenplots2;
            namepre = {'fig6f-','fig6g-','fig6h-'};
            namesuf = {'1'};
            f = @(x) (x);
            namesuf2 = {'','-log'};
            myfun = {f, @log10};
            figdescstr1 = {'Analytic', 'Asymptotic-Permuted', 'Asymptotic'};
            figdescstr2 = {'Dot Product Absolute Value'};
            figdescstr3 = {'','Log of '};
            mat = {obj.eigvec_dot_comparison_mat_ana2ana(1:n,1:n),...
                obj.eigvec_dot_comparison_mat_asy_permuted2asy_permuted(1:n,1:n),...
                obj.eigvec_dot_comparison_mat_asy2asy(1:n,1:n)};
            for i1 = 1:size(namepre,2)
                for i2 = 1:size(namesuf,2)
                    for i3 = 1:size(namesuf2,2)
                        name = [namepre{i1} namesuf{i2} namesuf2{i3}];
                        figdesc = [figdescstr1{i1} ' Jacobian Eigenvectors ' figdescstr3{i3} figdescstr2{i2} ' Comparison Matrix'];
                        f = figure('Name',name,'NumberTitle','off');
                        image(myfun{i3}(mat{i1}),'CDatamapping','scaled');
                        colorbar;
                        ylabel([figdescstr1{i1} ' Eigenvectors']);
                        xlabel([figdescstr1{i1} ' Eigenvectors']);
                        title({[name ' ' obj.scenarioName];obj.desc;figdesc;'$\mid v_i \cdot v_j \mid$'},'interpreter','latex');
                        obj.save_fig(f,name);
                    end
                end
            end
            %% fig6i* histogram
            namepre = {'fig6i'}; namesuf1 = {'-1'}; namesuf2 = {'','-log'};
            f = @(x) (x); myfun = {f, @log10};
            figdescstr1 = {'Absolute Value of Eigenvector Dot Products'};
            figdescstr2 = {'', '-log'};
            labels = {'','\log'};
            mat = {obj.eigvec_dot_comparison_mat_ana2ana,...
                obj.eigvec_dot_comparison_mat_asy2asy,...
                obj.eigvec_dot_comparison_mat_ana2asy_permuted};
            for i1 = 1:size(namepre,2)
                for i2 = 1:size(namesuf1,2)
                    for i3 = 1:size(namesuf2,2)
                        name = [namepre{i1} namesuf1{i2} namesuf2{i3}];
                        figdesc = [figdescstr1{i1} figdescstr2{i3} ' Histogram'];
                        f = figure('Name',name,'NumberTitle','off');
                        for i4 = 1:size(mat,2)
                            mat_curr = mat{i4};
                            a = ones(size(mat_curr)) - eye(size(mat_curr,1));
                            a = a > 0;
                            mat_curr = mat_curr(a);
                            histogram(myfun{i3}(mat_curr));
                            hold on;
                        end
                        legend('$\mid v_{i,ana} \cdot v_{j,ana} \mid$',...
                            '$\mid v_{i,asy} \cdot v_{j,asy} \mid$',...
                            '$\mid v_{i,asy_perm} \cdot v_{j,asy_perm} \mid$',...
                            'Interpreter','latex');
                        xlabel(['$' labels{i3} '\mid v_i \cdot v_j \mid$'],'Interpreter','latex');
                        title({[name ' ' obj.scenarioName];obj.desc;figdesc;['$' labels{i3} '\mid v_i \cdot v_j \mid$']},'interpreter','latex');
                        obj.save_fig(f,name);
                    end
                end
            end
        end
        function obj = plot_eigenvectors2(obj,isFirst,isKinn,isPermuted,isLogLog,isBinned)
            %%% fig7a,fig7b, fig7c,fig7d, -permuted -loglog -binned
            eigvecana = obj.eigenvectors_ana;
            if isKinn
                fnumab = 'fig7b';
                fnumcd = 'fig7d';
                x = obj.ki_nn;
                if isBinned
                    x = obj.ki_nnbinned;
                    eigvecana = obj.eigenvectors_ana_binned_k;
                end
                figdescab = 'Eigenvector Elements vs Node Average Nearest Neighbor Degree';
                figdesccd = 'Eigenvector Elements Mean vs Node Average Nearest Neighbor Degree';
                xlab = 'k_{nn,i} - weighted';
            else
                fnumab = 'fig7a';
                fnumcd = 'fig7c';
                x = obj.degree_vector_weighted;
                if isBinned
                    x=obj.kbinned;
                    eigvecana = obj.eigenvectors_ana_binned_kinn;
                end
                figdescab = 'Eigenvector Elements vs Node Degree';
                figdesccd = 'Eigenvector Elements Mean vs Node Degree';
                xlab = 'k_i - weighted';
            end
            if isFirst
                suffix = ['-first ' num2str(obj.numeigenplots)];
                n = obj.numeigenplots;
            else
                n = size(obj.eigenvectors_ana,2);
                suffix = '';
            end
            if isPermuted
                suffix2 = '-permuted';
                eigvecasy = obj.eigenvectors_asy_permuted;
                if isBinned && ~isKinn
                    eigvecasy = obj.eigenvectors_asy_permuted_binned_k;
                elseif isBinned && isKinn
                    eigvecasy = obj.eigenvectors_asy_permuted_binned_kinn;
                end
            else
                suffix2 = '';
                eigvecasy = obj.eigenvectors_asy;
            end
            if isLogLog
                suffix3 = '-loglog';
                myplot = @loglog;
                abs4log = @abs;
                ylab = '\mid v_i \mid';
            else
                suffix3 = '';
                myplot = @plot;
                abs4log = @(x) (x);
                ylab = 'v_i';
            end
            if isBinned
                suffix4 = '-binned';
            else
                suffix4 = '';
            end
            CM = jet(n);
            %% fig7a fig7b
            name = [fnumab suffix suffix2 suffix3 suffix4];
            f = figure('Name',name,'NumberTitle','off');
            for i=1:n
                myplot(x,abs4log(eigvecana(:,i)),'.','MarkerSize',18,'Color',CM(i,:));
                %plot(log10(x),log10(abs4log(obj.eigenvectors_ana(:,i))),'.','MarkerSize',18,'Color',CM(i,:));
                hold on;
            end
            for i=1:n
                myplot(x,abs4log(eigvecasy(:,i)),'^','MarkerSize',12,'Color',CM(i,:));
                %plot(log10(x),log10(abs4log(eigvecasy(:,i))),'^','MarkerSize',12,'Color',CM(i,:));
                hold on;
            end
            xlabel(xlab);
            ylabel(ylab);
            legendStr = cell(1, 2*n);
            for i = 1:n
                str = ['v^{(' num2str(i) ')} analytical'];
                legendStr{i} = str;
            end
            for i = n+1:2*n
                str = ['v^{(' num2str(i-n) ')} asymptotic' suffix2];
                legendStr{i} = str;
            end
            legend(legendStr);
            title({[name ' ' obj.scenarioName];obj.desc;figdescab;'v_i^{(j)} = i^{th} element of j^{th} eigenvector'});
            obj.save_fig(f,name);
            %% fig7c fig7d
            name = [fnumcd suffix suffix2 suffix3 suffix4];
            f = figure('Name',name,'NumberTitle','off');
            hold on;
            myplot(x,abs4log(mean(eigvecana(:,1:n),2)),'.','MarkerSize',18);
            myplot(x,abs4log(mean(eigvecasy(:,1:n),2)),'^','MarkerSize',12);
            xlabel(xlab);
            ylabel('$\mid \bar{v_i} \mid$','Interpreter','Latex');
            legend('Analytic Eigenvectors Mean', 'Asymptotic Eigenvectors Mean');
            title({[name ' ' obj.scenarioName];obj.desc;figdesccd;'v_i^{(j)} = i^{th} element of j^{th} eigenvector'});
            obj.save_fig(f,name);
        end
        %% fig8a fig8b random states angles to SS
        function obj = plot_random_states(obj)
            myfactor = {1,2};
            myshift = {0,-1};
            namesuf = {' 1st quad',' all quads'};
            for ind1 = 1:2
                ss = obj.steady_state';
                ssn = ss/norm(ss);
                seed = obj.randSeed;
                rng(seed);
                n = 300;
                x0 = rand(6000,n)*myfactor{ind1}+myshift{ind1};
                nrm = vecnorm(x0);
                x0n = x0./nrm;
                pert0 = x0-ss;
                pert0nrm = vecnorm(pert0);
                pert0n = pert0./pert0nrm;
                ssn = zeros(size(x0n)) + ssn;
                x0n1 = zeros(size(x0n)) + x0n(:,1);
                x2 = rand(6000,n)*myfactor{ind1}+myshift{ind1};
                nrm2 = vecnorm(x2);
                x2n = x2./nrm2;
                dims = [2:10,20:10:100,200:100:1000,2000:1000:10000];
                %dimplotind = [28:36];
                dimplotind = [1,9,24, 32,36];
                angleset = zeros(n,size(dims,2));
                ind = 1;
                for i = dims
                    x1 = rand(i,n)*myfactor{ind1}+myshift{ind1};
                    x2 = rand(i,n)*myfactor{ind1}+myshift{ind1};
                    x1 = x1./vecnorm(x1);
                    x2 = x2./vecnorm(x2);
                    angleset(:,ind) = rad2deg(real(acos(dot(x1,x2))))';
                    ind = ind+1;
                end
                anglesetmean = mean(angleset);
                anglesetstd = std(angleset);
                x3 = rand(2,n)*myfactor{ind1}+myshift{ind1};
                nrm3 = vecnorm(x3);
                x3n = x3./nrm3;
                angles = rad2deg(real(acos(dot(x0n,ssn))));
                angles2 = rad2deg(real(acos(dot(x0n1,x0n))));
                angles3 = rad2deg(real(acos(dot(x0n,x2n))));
                angles4 = rad2deg(real(acos(dot(pert0n,ssn))));
                %% fig8a
                name = {['fig8a' namesuf{ind1}],['fig8f' namesuf{ind1}]};
                figdesc = {['Angles between random states' namesuf{ind1}],...
                    ['Distribution of angles between random states' namesuf{ind1}]};
                f = figure('Name',name{1},'NumberTitle','off');
                plot(angles,'ob');
                hold on;
                plot(angles2,'*r');
                plot(angles3,'.k','MarkerSize',12);
                plot(angles4,'^m');
                legendStr = {};
                legendStr{1} = '$\theta_{x_{rand},SS}$';
                legendStr{2} = '$\theta_{x_{rand,1},x_{rand}}$';
                legendStr{3} = '$\theta_{x_{rand},x_{rand}}$';
                legendStr{4} = '$\theta_{x_{rand}-SS,SS}$';
                for ind = dimplotind
                    dim = dims(ind);
                    anglesdim = angleset(:,ind);
                    plot(anglesdim,'.','MarkerSize',12)
                    legendStr{end+1} = ['$\theta_{x_{rand},x_{rand}}$ dim = ' num2str(dim)];
                end
                xlabel('Index');
                ylabel('$\theta$','Interpreter','Latex');
                legend(legendStr,'Interpreter','Latex','FontSize',14);
                title({[name{1} ' ' obj.scenarioName];obj.desc;figdesc{1}});
                obj.save_fig(f,name{1});
                %% fig8f
                f = figure('Name',name{2},'NumberTitle','off');
                histogram(angles);
                hold on;
                histogram(angles2(2:end));
                histogram(angles3);
                histogram(angles4);
                %                 legendStr = {};
                for ind = dimplotind
                    dim = dims(ind);
                    anglesdim = angleset(:,ind);
                    histogram(anglesdim)
                    %                     legendStr{end+1} = ['$\theta_{x_{rand},x_{rand}}$ dim = ' num2str(dim)];
                end
                xlabel('$\theta$','Interpreter','Latex');
                %             ylabel('$\theta$','Interpreter','Latex');
                legend(legendStr,'Interpreter','Latex','FontSize',14);
                title({[name{2} ' ' obj.scenarioName];obj.desc;figdesc{2}});
                obj.save_fig(f,name{2});
                %% fig8b fig8c fig8d fig8e
                name = {['fig8b' namesuf{ind1}],['fig8c' namesuf{ind1}],...
                    ['fig8d' namesuf{ind1}],['fig8e' namesuf{ind1}]};
                figdesc = {['Average angle between random states vs dimension' namesuf{ind1}],...
                    ['Average angle between random states vs log(dimension)' namesuf{ind1}],...
                    ['Standard Deviation of angles between random states vs dimension' namesuf{ind1}],...
                    ['Standard Deviation of angles between random states vs log(dimension)' namesuf{ind1}]};
                x = {dims, log10(dims)};
                xlab = {'Dimension','log(Dimension)'};
                y = {anglesetmean,anglesetstd};
                ylab = {'$\overline{\theta}$','$\sigma_{\theta}$'};
                for ind = [1:4]
                    f = figure('Name',name{ind},'NumberTitle','off');
                    plot(x{mod(ind-1,2)+1},y{(ind>2) + 1},'ob');
                    xlabel(xlab{mod(ind-1,2)+1});
                    ylabel(ylab{(ind>2)+1},'Interpreter','Latex');
                    %legend(legendStr,'Interpreter','Latex','FontSize',14);
                    title({[name{ind} ' ' obj.scenarioName];obj.desc;figdesc{ind}});
                    obj.save_fig(f,name{ind});
                end
            end
        end
        %% fig9*
        function plot_network1(obj)
            sizes = [.01 .1 .5 .75 1 ...
                1.25 1.5  1.75  2 2.25 ...
                2.5  2.75  3  5 10];
            %colors = jet(obj.numbins);
            %             layouts = {'force','force3','subspace','subspace3'};
            %             linestyle = {'none','-'};
            layouts = {'force','subspace'};
            linestyle = {'-'};
            linewidthsstr = {'none','0.1'};
            coloringstr = {'Degree','Steady State','abs(Eigvec_{ana,', 'abs(Eigvec_{asy,'};
            sizestr = {'Degree'};
            tickinds = [1 10 11 12 13 14 15];
            ssbins = obj.set_bins_generic(obj.numbins,1-obj.steady_state',1e-13,true(obj.N,1));
            ssbinned = obj.set_binned_vals(1-obj.steady_state',ssbins);
            %             sstickvals = ssbinned(tickinds);sstickvals = sstickvals(end:-1:1);
            %             sstickstr = string(ssbinned(tickinds));sstickstr = sstickstr(end:-1:1);
            %             ssbins = ssbins(end:-1:1);
            %             ssbinned = 1-ssbinned(end:-1:1);
            colors_list{1} = {jet(size(obj.bins,1))};
            colorstemp = jet(size(ssbins,1));
            colors_list{2} = {colorstemp(end:-1:1,:)};
            caxislim_2 = [-log10(ssbinned(end)) -log10(ssbinned(1))];
            
            for i=1:size(obj.binseigvecana,1)
                binseigvecana{i} = obj.binseigvecana{i,1}(obj.eigenvectors_ana_binned_self{i}~=0);
                valseigvecana{i} = obj.eigenvectors_ana_binned_self{i}(obj.eigenvectors_ana_binned_self{i}~=0);
                binseigvecasy{i} = obj.binseigvecasy_permuted{i,1}(obj.eigenvectors_asy_permuted_binned_self{i}~=0);
                valseigvecasy{i} = obj.eigenvectors_asy_permuted_binned_self{i}(obj.eigenvectors_asy_permuted_binned_self{i}~=0);
                valseigvec{i} = sort([valseigvecana{i};valseigvecasy{i}],'ascend');
                colors34{i} = jet(size(valseigvec{i},1));
                caxis_lim34{i} = [log10(valseigvec{i}(1)) log10(valseigvec{i}(end))];
            end
            colors_list{3} = colors34;colors_list{4} = colors34;
            caxis_lim_list = {{[obj.kbinned(1),obj.kbinned(end)]},{caxislim_2},caxis_lim34,caxis_lim34};
            sym x; tickstrfun_list = {@(x) (x), @(x) (1-10.^(-x)), @(x) (10.^x), @(x) (10.^x)};
            coloringbins = {{obj.bins},{ssbins},binseigvecana,binseigvecasy};
            valscolorbar_list = {{obj.kbinned},{ssbinned},valseigvec,valseigvec};
            binnedcolorings = {{obj.kbinned},{ssbinned},valseigvecana,...
                valseigvecasy};
            %tickstr = compose("%.2e",tickvals);
            a = obj.adjacencyMatrix;
            a_uw = (a>0); %uw = unweighted
            g_uw = graph(a_uw);
            p = {};
            ind = 1;
            for i1 = 1:length(layouts)
                for i2 = 1:length(linestyle)
                    for i3 = 1:length(coloringbins)
                        tickstrfun = tickstrfun_list{i3};
                        for i4 = 1:length(coloringbins{i3})
                            colorbins = coloringbins{i3}{i4};
                            vals = binnedcolorings{i3}{i4};
                            colors = colors_list{i3}{i4};
                            valscolorbar = valscolorbar_list{i3}{i4};
                            name = ['fig9-' num2str(i1) '-' num2str(i2) '-' num2str(i3) '-' num2str(i4)];
                            if length(coloringbins{i3})>1
                                str1 = [num2str(i4) '})'];
                            else
                                str1 = '';
                            end
                            figdesc = ['Network Visualization ' layouts{i1} ' layout, Size: ' ...
                                sizestr{1} ', Coloring: ' ...
                                coloringstr{i3} str1];
                            f = figure('Name',name,'NumberTitle','off');
                            colormap jet;
                            p{ind} = plot(g_uw,'LineWidth',.1,'LineStyle',linestyle{i2},'Marker','o','layout',layouts{i1});
                            for i = 1:size(obj.bins,1)
                                inds1 = obj.bins{i};
                                highlight(p{ind},inds1,'MarkerSize',sizes(i))
                            end
                            for i = 1:size(colorbins,1)
                                inds2 = colorbins{i};
                                highlight(p{ind},inds2,'NodeColor',colors(vals(i)==valscolorbar,:));
                            end
                            caxis_lim =caxis_lim_list{i3}{i4};
                            tickvals = caxis_lim(1):(caxis_lim(2)-caxis_lim(1))/10:caxis_lim(2);
                            caxis(caxis_lim);
                            ticks = tickvals;
                            %ticks = [binnedcolorings{i3}(1)', binnedcolorings{i3}(end-1:end)'];
                            ticksstr = compose("%.2e",tickstrfun(ticks));
                            colorbar('Ticks',ticks,'TickLabels',ticksstr);
                            title({[name ' ' obj.scenarioName];obj.desc;figdesc});
                            obj.save_fig(f,name);
                            ind = ind+1;
                        end
                    end
                end
            end
        end
        function obj = plot_degree(obj)
            %% fig10a
            names = 'ab';
            figdescs = {'Node Degrees','Node Average Nearest Neighbor Degrees'};
            values_list = {obj.degree_vector_weighted,obj.ki_nn};
            legendStrs = {'k_i','k_{i,nn}'};
            bins_list = {obj.bins,obj.binskinn};
            binned_list = {obj.kbinned,obj.ki_nnbinned};
            ylabstr = {'k_i','k_{i,nn}'};
            myplots = {@loglog,@semilogx};
            numfigs = length(names);
            for figind1 = 1:numfigs
                name = ['fig10' names(figind1)];
                figdesc = figdescs{figind1};
                f = figure('Name',name,'NumberTitle','off');
                myplot = myplots{figind1};
                values = values_list{figind1};
                myplot(values,'.');
                legendStr1 = legendStrs{figind1};
                legendStr{1} = legendStr1;
                hold on;
                values_sorted = sort(values,'descend');
                myplot(values_sorted,'.r');
                legendStr{2} = [legendStr1 ' - Sorted'];
                binsizes = [0];
                cumbinsizes = [0];
                bins1 = bins_list{figind1};
                for i = obj.numbins:-1:1
                    binsizes(end+1) =  size(bins1{i},1);
                    cumbinsizes(end+1) = cumbinsizes(end) + size(bins1{i},1);
                end
                midbins = (cumbinsizes(1:end-1) + cumbinsizes(2:end))/2;
                %             halfsizes = (binsizes(1:end-1) + binsizes(2:end))/2;
                %             midbins = binsizes(1:end-1) + halfsizes;
                binned = binned_list{figind1}(end:-1:1);
                myplot(midbins,binned,'or','MarkerSize',14);
                legendStr{3} = [legendStr1 ' - Binned'];
                title({[name ' ' obj.scenarioName];obj.desc;figdesc});
                legend(legendStr);
                xlabel('i');
                ylabel(ylabstr{figind1});
                obj.save_fig(f,name);
            end
        end
        function obj = plot_eigvec_map(obj)
            structs{1} = obj.pert_map_struct;
            m.new_epsilons = obj.eps_adjusted(1,7);
            m.is_inits_legit = [1];
            m.solutions_t = obj.solution_t_eigvecasy(1,2,1);
            m.solutions_x = obj.solution_x_eigvecasy(1,2,1);
            structs{2} = m;
            [f1,hax] = obj.plot_eigvec_map1(structs,obj.eigenvectors_ana(:,1),...
                obj.eigenvectors_ana(:,2),obj.steady_state');
        end
        function [f, hax] = plot_eigvec_map1(obj,map_structs,v1,v2,ss)
            numfigs = 8;
            for i = 1:numfigs
                name = ['fig11-' num2str(i)];
                f{i} = figure('Name',name,'NumberTitle','off');
                hax{i} = axes;
            end
            m1 = size(map_structs,2);
            markers = ['.*'];
            legendStr = {};
            legendInd = 1;
            for i1 = 1:m1
                map_struct = map_structs{1,i1};
                n1 = size(map_struct.is_inits_legit,2);
                marker = markers(i1);
                first = true;
                for i = 1:n1
                    if map_struct.is_inits_legit(1,i)
                        if isfield(map_struct,'factors1')
                            if i<=length(map_struct.factors1)
                                str = ['x_0 = ' num2str(map_struct.factors1(i)) '*' map_struct.v1str ' + ' map_struct.v2str];
                            else
                                str = ['x_0 = ' map_struct.v2str ' + ' num2str(map_struct.factors2(i-length(map_struct.factors1))) '*' map_struct.v1str];
                            end
                            legendStr{legendInd} = str;
                        else
                            legendStr{legendInd} = map_struct.legendStr{i};
                        end
                        legendInd = legendInd + 1;
                        pert_t = map_struct.solutions_t{1,i};
                        pert_x = map_struct.solutions_x{1,i}' - ss;
                        epsilon = map_struct.new_epsilons(1,i);
                        v1mat = v1 + zeros(size(pert_x));
                        pert_x_norm = vecnorm(pert_x).*ones(size(pert_x));
                        pert_x_normalized = pert_x./pert_x_norm;
                        pert_x_normalized_1 = pert_x_normalized * epsilon;
                        dotv1 = dot(v1mat,pert_x_normalized);
                        dotv1_1 = dot(v1mat,pert_x_normalized_1);
                        dotv1_2 = dot(v1mat,pert_x);
                        v2mat = v2 + zeros(size(pert_x));
                        dotv2 = dot(v2mat,pert_x_normalized);
                        dotv2_1 = dot(v2mat,pert_x_normalized_1);
                        dotv2_2 = dot(v2mat,pert_x);
                        x = {pert_t,pert_t,dotv1,dotv1_1,pert_t,pert_t,dotv1_2,dotv1_2};
                        y = {dotv1,dotv2,dotv2,dotv2_1,dotv1,dotv2,dotv2_2,dotv2_2};
                        myplot = {@plot,@plot,@plot,@plot,@semilogy,@semilogy,@plot,@loglog};
                        for j = 1:numfigs
                            myplot{j}(hax{j},x{j},y{j},marker);
                            if first
                                hold(hax{j},'on');
                            end
                            legend(hax{j},legendStr);
                        end
                        first = false;
                    end
                end
            end
            figdesc = {'Normalized perturbation projection on v_{1,ana} vs time',...
                'Normalized perturbation projection on v_{2,ana} vs time',...
                'Normalized perturbation projection',...
                'Normalized perturbation projection scaled to x_0',...
                'Normalized perturbation projection on v_{1,ana} vs time (log)',...
                'Normalized perturbation projection on v_{2,ana} vs time (log)',...
                'Perturbation projection','Perturbation projection (log)'};
            str1 = '$\textit{proj}_{\textbf{v}_{1,ana}}\widehat{\textbf{pert}}$';
            str2 = '$\textit{proj}_{\textbf{v}_{2,ana}}\widehat{\textbf{pert}}$';
            str3 = '$\textit{proj}_{\textbf{v}_{1,ana}}\textbf{pert}$';
            str4 = '$\textit{proj}_{\textbf{v}_{2,ana}}\textbf{pert}$';
            xlabs = {'t','t',str1,str1,'t','t',str3,str3};
            ylabs = {str1,str2,str2,str2,str1,str2,str4,str4};
            for j = 1:numfigs
                name = ['fig11-' num2str(j)];
                title(hax{j},{[name ' ' obj.scenarioName];obj.desc;figdesc{j}});
                xlabel(hax{j},xlabs{j},'Interpreter','latex');
                ylabel(hax{j},ylabs{j},'Interpreter','latex');
                obj.save_fig(f{j},name);
            end
        end
        function pert_approx_x = find_nth_order_approx(~,sol_t,pert,eigvecs,eigvals,order)
            approx = 0;
            pert_0 = pert(1,:);
            pert_approx_x = cell(1,order);
            for i = 1:order
                disp(['i = ' num2str(i)]);
                eigvec_i = eigvecs(:,i);
                eigval_i = eigvals(i);
                approx = approx + dot(pert_0,eigvec_i)*eigvec_i*exp(eigval_i*sol_t');
                pert_approx_x{1,i} = approx';
            end
        end
        function obj = plot_pert_approx(obj,drawNodes,drawBinned,drawWholeSystem,drawDecayTime)
            %% fig 13*
            mydata = load(fullfile(obj.resultsPath,'obj_properties','eigvals_asy_v2_set.mat'),'var');eigvals_asy_v2_set=mydata.var;
            mydata = load(fullfile(obj.resultsPath,'obj_properties','eigvecs_asy_v2_set.mat'),'var');eigvecs_asy_v2_set=mydata.var;
            namepre = 'fig13-';
            [~,node_maxk_ind] = max(obj.degree_vector_weighted);
            [~,node_avgk_ind] = min(abs(obj.degree_weighted_average - obj.degree_vector_weighted));
            [~,node_mink_ind] = min(obj.degree_vector_weighted);
            [~,node_maxkinn_ind] = max(obj.ki_nn);
            [~,node_knn_ind] = min(abs(obj.k_nn - obj.ki_nn));
            [~,node_minkinn_ind] = min(obj.ki_nn);
            node_ids = [node_maxk_ind,node_avgk_ind,node_mink_ind,node_maxkinn_ind,...
                node_knn_ind,node_minkinn_ind];
            node_bins = {obj.bins, obj.binskinn};
            sol_t_runs = {{obj.solution_t},obj.solution_t_eigvecana,obj.solution_t_eigvecasy,...
                obj.solution_t_perturbations};
            sol_x_runs = {{obj.solution_x},obj.solution_x_eigvecana,obj.solution_x_eigvecasy,...
                obj.solution_x_perturbations};
            figdesc1 = {'pert = rand', 'pert = v_{i,ana}', 'pert = v_{i,asy}', 'pert = .1*rand*ss'};
            figdesc2 = {'max(k_j)', 'mean(k_j)', 'min(k_j)', 'max(k_{j,nn})',...
                'mean(k_{j,nn})', 'min(k_{j,nn})','max(pert)'};
            figdesc1b = {'k', 'k_{nn}'};
            % pert_0 = rand, obj.solution_t,...
            num_nodes = size(node_ids,2);num_runs = size(sol_t_runs,2);
            num_node_bins = size(node_bins,2);
            num_runs_str = 'abcdefghjijklmnopqrstuvwxyz';
            eigvecs_sets = {obj.eigenvectors_ana,obj.eigenvectors_asy_permuted,eigvecs_asy_v2_set{14,1},eigvecs_asy_v2_set{15,1},obj.eigvecs_v3};
            eigvals_sets = {obj.eigenvalues_ana,obj.eigenvalues_asy_permuted,eigvals_asy_v2_set{14,1},eigvals_asy_v2_set{15,1},obj.eigvals_v3};
            num_eigvec_sets = size(eigvecs_sets,2);
            eigvecsets_str = {'ana','asy','asy-v2-14','asy-v2-15','asy-v3'};
            order = 5;
            linestyles = {'--',':','-.','-.','-'};markers = {'o','*','+','x','s'};
            colors = {'r','k','m','b','g'};step = 100;
            myplotfuns = {@plot,@semilogy};numplotfuns = size(myplotfuns,2);plottypestr = {'','-logabs'};
            myabsfuns = {@(x) x, @(x) abs(x)};
            norms = [1 2];
            num_norms = size(norms,2);
            for i1 = 1:num_runs
                sol_ts = sol_t_runs{1,i1};sol_xs = sol_x_runs{1,i1};
                num_sols = size(sol_ts,2);
                for i2 = 1:num_sols
                    if ~isempty(sol_ts{1,i2})
                        if (i1 == 2) || (i1 == 3)
                            figdesc3 = ['i = ' num2str(i2)];
                        else
                            figdesc3 = '';
                        end
                        sol_t = sol_ts{1,i2}; sol_x = sol_xs{1,i2};
                        pert = (sol_x' - obj.steady_state')';
                        pert0 = pert(1,:);
                        [~, node_maxpert_ind] = max(pert(1,:));
                        node_ids(1,num_nodes + 1) = node_maxpert_ind;
                        if drawNodes
                            for i3 = 1:num_nodes + 1
                                node_id = node_ids(i3);
                                %                         for i3 = num_nodes + 1
                                for i6 = 1:numplotfuns
                                    name = [namepre num_runs_str(i1) '-' num2str(i2) '-node' num2str(node_id)...
                                        plottypestr{i6}];
                                    f = figure('Name',name,'NumberTitle','off');
                                    legendStr = {}; legendInd = 1;
                                    myplotfuns{i6}(sol_t,myabsfuns{i6}(pert(:,node_id)),'-','LineWidth',3);
                                    hold on;
                                    legendStr{legendInd} = 'Numerical Solution (non-linear)';
                                    legendInd = legendInd + 1;
                                    for i4 = 1:num_eigvec_sets
                                        disp(['i4 = ' num2str(i4)]);
                                        eigvecs = eigvecs_sets{1,i4};
                                        eigvals = eigvals_sets{1,i4};
                                        pert_approx_x = obj.find_nth_order_approx(sol_t,pert,eigvecs,eigvals,order);
                                        num_pert_approx_x = size(pert_approx_x,2);
                                        for i5 = 1:num_pert_approx_x
                                            pert_approx_x_curr = pert_approx_x{1,i5};
                                            myplotfuns{i6}(sol_t(1:step:end),myabsfuns{i6}(pert_approx_x_curr(1:step:end,node_id)),'LineStyle',...
                                                linestyles{i4},'Marker',markers{i5},'Color',colors{i4},'MarkerEdgeColor',colors{i4});
                                            legendStr{legendInd} = ['Order ' num2str(i5) ' approx, v_{ ' eigvecsets_str{i4} '}'];
                                            legendInd = legendInd + 1;
                                        end
                                    end
                                    xlabel('t');
                                    ylabel(['x_{' num2str(node_id) '}']);
                                    figdesc = ['Linear approximations,' figdesc1{i1} ', ' figdesc3 ', node with ' figdesc2{i3}];
                                    figdesc4 =  ['k_i = ' num2str(obj.degree_vector_weighted(node_id)) ', k_{i,nn} = ' num2str(obj.ki_nn(node_id))];
                                    title({[name ' ' obj.scenarioName];obj.desc;figdesc;figdesc4});
                                    legend(legendStr);
                                    obj.save_fig(f,name);
                                end
                            end
                        end
                        if drawBinned
                            for i7 = 1:num_node_bins
                                curr_bins = node_bins{1,i7};
                                num_curr_bins = size(curr_bins,1);
                                curr_bins_desc = figdesc1b{1,i7};
                                for i8 = 1:num_curr_bins
                                    name = [namepre num_runs_str(i1) '-' num2str(i2) '-binnedby-' curr_bins_desc...
                                        'bin-' num2str(i8) plottypestr{2}];
                                    f = figure('Name',name,'NumberTitle','off');
                                    legendStr = {}; legendInd = 1;
                                    curr_bin = curr_bins{i8,1};
                                    pert_binned = mean(pert(:,curr_bin),2);
                                    myplotfuns{2}(sol_t,myabsfuns{2}(pert_binned),'-','LineWidth',3);
                                    hold on;
                                    legendStr{legendInd} = 'Numerical Solution';
                                    legendInd = legendInd + 1;
                                    for i9 = 1:num_eigvec_sets
                                        disp(['i9 = ' num2str(i9)]);
                                        eigvecs = eigvecs_sets{1,i9};
                                        eigvals = eigvals_sets{1,i9};
                                        pert_approx_x = obj.find_nth_order_approx(sol_t,pert,eigvecs,eigvals,order);
                                        num_pert_approx_x = size(pert_approx_x,2);
                                        pert_approx_x_curr = pert_approx_x{1,num_pert_approx_x};
                                        pert_approx_x_curr_binned = mean(pert_approx_x_curr(:,curr_bin),2);
                                        myplotfuns{2}(sol_t(1:step:end),myabsfuns{2}(pert_approx_x_curr_binned(1:step:end,:)),'LineStyle',...
                                            linestyles{i9},'Marker',markers{num_pert_approx_x},'Color',colors{i9},'MarkerEdgeColor',colors{i9});
                                        legendStr{legendInd} = ['Order ' num2str(num_pert_approx_x) ' approx, v_{ ' eigvecsets_str{i9} '}'];
                                        legendInd = legendInd + 1;
                                    end
                                    xlabel('t');
                                    ylabel(['x']);
                                    figdesc = ['Linear approximations,' figdesc1{i1} ', ' figdesc3 ', binned by ' curr_bins_desc ', bin ' num2str(i8)];
                                    figdesc4 =  ['Num nodes in bin: ' num2str(size(curr_bin,1)),...
                                        ', k_{min} = ' num2str(min(obj.degree_vector_weighted(curr_bin))),...
                                        ', k_{max} = ' num2str(max(obj.degree_vector_weighted(curr_bin))),...
                                        ', k_{nn,min} = ' num2str(min(obj.ki_nn(curr_bin))),...
                                        ', k_{nn,max} = ' num2str(max(obj.ki_nn(curr_bin)))];
                                    title({[name ' ' obj.scenarioName];obj.desc;figdesc;figdesc4});
                                    legend(legendStr);
                                    obj.save_fig(f,name);
                                end
                            end
                        end
                        if drawWholeSystem
                            for i10 = 1:num_norms
                                normP = norms(i10);
                                pert_normed = vecnorm(pert',normP)';
                                name = [namepre num_runs_str(i1) '-' num2str(i2) '-WholeSystemNorm-L' num2str(normP)...
                                    '-' plottypestr{2}];
                                f = figure('Name',name,'NumberTitle','off');
                                legendStr = {}; legendInd = 1;
                                myplotfuns{2}(sol_t,myabsfuns{2}(pert_normed),'-','LineWidth',3);
                                hold on;
                                legendStr{legendInd} = 'Numerical Solution';
                                legendInd = legendInd + 1;
                                for i9 = 1:num_eigvec_sets
                                    %disp(['i9 = ' num2str(i9)]);
                                    eigvecs = eigvecs_sets{1,i9};
                                    eigvals = eigvals_sets{1,i9};
                                    pert_approx_x = obj.find_nth_order_approx(sol_t,pert,eigvecs,eigvals,order);
                                    num_pert_approx_x = size(pert_approx_x,2);
                                    pert_approx_x_curr = pert_approx_x{1,num_pert_approx_x};
                                    pert_approx_x_curr_normed = vecnorm(pert_approx_x_curr',normP)';
                                    myplotfuns{2}(sol_t(1:step:end),myabsfuns{2}(pert_approx_x_curr_normed(1:step:end,:)),'LineStyle',...
                                        linestyles{i9},'Marker',markers{num_pert_approx_x},'Color',colors{i9},'MarkerEdgeColor',colors{i9});
                                    legendStr{legendInd} = ['Order ' num2str(num_pert_approx_x) ' approx, v_{ ' eigvecsets_str{i9} '}'];
                                    legendInd = legendInd + 1;
                                end
                                xlabel('t');
                                ylabel(['|x|_{L' num2str(normP) '}']);
                                figdesc = ['Linear approximations,' figdesc1{i1} ', ' figdesc3 ', System L' num2str(normP) ' norm'];
                                title({[name ' ' obj.scenarioName];obj.desc;figdesc});
                                legend(legendStr);
                                obj.save_fig(f,name);
                            end
                        end
                        if drawDecayTime
                            xvalues = {obj.degree_vector_weighted, obj.ki_nn};
                            num_xvalues = size(xvalues,2);
                            xlabelstrsf = {'ki','kinn'};
                            xlabelstrs = {'k_i','k_{i,nn}'};
                            nodesDecayTime = obj.find_decay_times(pert,sol_t);
                            for i11 = 1:num_xvalues
                                xvalue = xvalues{i11};xlabelstr = xlabelstrs{i11};
                                xlabelstrf = xlabelstrsf{i11};
                                name = [namepre num_runs_str(i1) '-' num2str(i2) '-DecayTimeVs' xlabelstrf...
                                    ];
                                f = figure('Name',name,'NumberTitle','off');
                                legendStr = {}; legendInd = 1;
                                semilogx(xvalue,nodesDecayTime,'s');
                                hold on;
                                legendStr{legendInd} = 'Numerical Solution';
                                legendInd = legendInd + 1;
                                for i9 = 1:num_eigvec_sets
                                    %disp(['i9 = ' num2str(i9)]);
                                    eigvecs = eigvecs_sets{1,i9};
                                    eigvals = eigvals_sets{1,i9};
                                    pert_approx_x = obj.find_nth_order_approx(sol_t,pert,eigvecs,eigvals,order);
                                    num_pert_approx_x = size(pert_approx_x,2);
                                    pert_approx_x_curr = pert_approx_x{1,num_pert_approx_x};
                                    pert_approx_x_curr_decayTimes = obj.find_decay_times(pert_approx_x_curr,sol_t);
                                    semilogx(xvalue,pert_approx_x_curr_decayTimes,'^','MarkerEdgeColor',colors{i9});
                                    legendStr{legendInd} = ['Order ' num2str(num_pert_approx_x) ' approx, v_{ ' eigvecsets_str{i9} '}'];
                                    legendInd = legendInd + 1;
                                end
                                xlabel(xlabelstr);
                                ylabel(['t_{1/e}']);
                                figdesc = ['Linear approximations,' figdesc1{i1} ', ' figdesc3 ', Node Decay Times vs ' xlabelstr];
                                title({[name ' ' obj.scenarioName];obj.desc;figdesc});
                                legend(legendStr);
                                obj.save_fig(f,name);
                            end
                        end
                    end
                end
            end
        end
        %%
        function obj = plot_eigenvectors3(obj)
            mydata = load(fullfile(obj.resultsPath,'obj_properties','eigvec_dot_comparison_mat_asy2asy-v2-14.mat'),'var');
            eigvec_dot_comparison_mat_asy2asy_v2_14 = mydata.var;
            mydata = load(fullfile(obj.resultsPath,'obj_properties','eigvec_dot_comparison_mat_asy2asy-v2-15.mat'),'var');
            eigvec_dot_comparison_mat_asy2asy_v2_15 = mydata.var;
            mydata = load(fullfile(obj.resultsPath,'obj_properties','eigvec_dot_comparison_mat_asy2asy-v3.mat'),'var');
            eigvec_dot_comparison_mat_asy2asy_v3 = mydata.var;
            mats = {eigvec_dot_comparison_mat_asy2asy_v2_14,...
                eigvec_dot_comparison_mat_asy2asy_v2_15,...
                eigvec_dot_comparison_mat_asy2asy_v3};
            namepre = 'fig14-';
            figdesc1 = 'Dot product of Absolute value of eigenvectors';
            xlabelstrs = {'Asymptotic Jacobian v2 bin 14 Eigenvectors',...
                'Asymptotic Jacobian v2 bin 15 Eigenvectors',...
                'Asymptotic Jacobian v3 Eigenvectors'};
            num_mats = size(mats,2);
            for i1 = 1:num_mats
                name = [namepre num2str(i1)];
                figdesc = [figdesc1];
                f = figure('Name',name,'NumberTitle','off');
                image(mats{i1},'CDatamapping','scaled');
                colorbar;
                ylabel('Asymptotic Jacobian Eigenvectors');
                xlabel(xlabelstrs{i1});
                title({[name ' ' obj.scenarioName];obj.desc;figdesc;'$\mid v_i \cdot v_j \mid$'},'interpreter','latex');
                obj.save_fig(f,name);
            end
        end
        %% fig17*
        function obj = plot_eigenvectors4(obj)
            for i=1:4
                str = ['mydata = load(fullfile(obj.resultsPath,''obj_properties'',''eigvec_ana_asy_perm_weighted_dot_products_' num2str(i) '''),''var'');'];
                eval(str);
                str = ['wdp' num2str(i) '= mydata.var;'];
                eval(str);
                str = ['mydata = load(fullfile(obj.resultsPath,''obj_properties'',''eigvec_ana_asy_perm_weighted_dot_products_mean_' num2str(i) '''),''var'');'];
                eval(str);
                str = ['wdpm' num2str(i) '= mydata.var;'];
                eval(str);
            end
            %             mydata = load(fullfile(obj.resultsPath,'obj_properties','eigvec_ana_asy_perm_weighted_dot_products_1'),'var');
            %             wdp1 = mydata.var;
            %             mydata = load(fullfile(obj.resultsPath,'obj_properties','eigvec_ana_asy_perm_weighted_dot_products_mean_1'),'var');
            %             wdpm1 = mydata.var;
            %             mydata = load(fullfile(obj.resultsPath,'obj_properties','eigvec_ana_asy_perm_weighted_dot_products_2'),'var');
            %             wdp2 = mydata.var;
            %             mydata = load(fullfile(obj.resultsPath,'obj_properties','eigvec_ana_asy_perm_weighted_dot_products_mean_2'),'var');
            %             wdpm2 = mydata.var;
            %             mydata = load(fullfile(obj.resultsPath,'obj_properties','eigvec_ana_asy_perm_weighted_dot_products_3'),'var');
            %             wdp3 = mydata.var;
            %             mydata = load(fullfile(obj.resultsPath,'obj_properties','eigvec_ana_asy_perm_weighted_dot_products_mean_3'),'var');
            %             wdpm3 = mydata.var;
            name = 'fig17-1';
            f = figure('Name',name,'NumberTitle','off');
            figdesc = 'Weighted dot products (using only $k_i > k_0$) of corresponding eigenvectors';
            %            semilogx(obj.eigenvalues_ana,abs(wdp1),'o');
            m = size(wdp1,2);
            legendStr1 = cell(m,1);
            for i=1:m
                legendStr1{i,1} = ['k_0 = ' num2str(obj.kbinned_mins(i))];
            end
            legend(legendStr1);
            xlabel('$\lambda$','Interpreter','latex');
            ylabel('$\mid v \cdot_{>k_0} \hat{v}\mid$','Interpreter','latex');
            title({[name ' ' obj.scenarioName];obj.desc;figdesc},'interpreter','latex');
            obj.save_fig(f,name);
            %%%%
            name = 'fig17-2';
            f = figure('Name',name,'NumberTitle','off');
            figdesc = 'Weighted dot products (using only $k_i < k_0$) of corresponding eigenvectors';
            %            semilogx(obj.eigenvalues_ana,abs(wdp2),'o');
            m = size(wdp2,2);
            legendStr2 = cell(m,1);
            for i=1:m
                legendStr2{i,1} = ['k_0 = ' num2str(obj.kbinned_maxs(i))];
            end
            legend(legendStr2);
            xlabel('$\lambda$','Interpreter','latex');
            ylabel('$\mid v \cdot_{<k_0} \hat{v}\mid$','Interpreter','latex');
            title({[name ' ' obj.scenarioName];obj.desc;figdesc},'interpreter','latex');
            obj.save_fig(f,name);
            %%%%%%
            name = 'fig17-3';
            f = figure('Name',name,'NumberTitle','off');
            figdesc = 'Weighted dot products (using only $k_i > k_0$) of corresponding eigenvectors historgrams';
            for i=1:m
                histogram(abs(wdp1(:,i)));
                hold on;
            end
            legend(legendStr1);
            title({[name ' ' obj.scenarioName];obj.desc;figdesc},'interpreter','latex');
            obj.save_fig(f,name);
            %%%%%%
            name = 'fig17-4';
            f = figure('Name',name,'NumberTitle','off');
            figdesc = 'Weighted dot products (using only $k_i < k_0$) of corresponding eigenvectors historgrams';
            for i=1:m
                histogram(abs(wdp2(:,i)));
                hold on;
            end
            legend(legendStr2);
            title({[name ' ' obj.scenarioName];obj.desc;figdesc},'interpreter','latex');
            obj.save_fig(f,name);
            %%%%%
            name = 'fig17-5';
            f = figure('Name',name,'NumberTitle','off');
            figdesc = 'Average Weighted dot products of corresponding eigenvectors';
            image(wdpm1,'CDatamapping','scaled');
            colorbar;
            title({[name ' ' obj.scenarioName];obj.desc;figdesc},'interpreter','latex');
            xlabel('k_0');
            s = size(wdpm1,1);n = 100/s;ytickstep = 5;
            xticks(1:15);xticklabels(strsplit(num2str(obj.kbinned_mins')));
            yticks(ytickstep:ytickstep:s);yticklabels(strsplit(num2str(s - ytickstep:-ytickstep:ytickstep)));
            ylabel('Percentile of eigenvectors (0 = all, 99 = first 60)');
            obj.save_fig(f,name);
            %%%%%
            name = 'fig17-7';
            f = figure('Name',name,'NumberTitle','off');
            figdesc = 'Average Weighted dot products of corresponding eigenvectors';
            image(wdpm3,'CDatamapping','scaled');
            colorbar;
            title({[name ' ' obj.scenarioName];obj.desc;figdesc},'interpreter','latex');
            s = size(wdpm3,1);n = 100/s;ytickstep = 5;
            xlabel('Node degree percentile (0 = all, 90 = hubs)');
            xticks(1:20);xticklabels(strsplit(num2str(5:5:100)));
            yticks(ytickstep:ytickstep:s);yticklabels(strsplit(num2str(s - ytickstep:-ytickstep:ytickstep)));
            ylabel('Percentile of eigenvectors (0 = all, 99 = first 60)');
            obj.save_fig(f,name);
            %%%%%
            name = 'fig17-8';
            f = figure('Name',name,'NumberTitle','off');
            figdesc = 'Average Weighted dot products of corresponding eigenvectors (reverse order)';
            image(wdpm4,'CDatamapping','scaled');
            colorbar;
            title({[name ' ' obj.scenarioName];obj.desc;figdesc},'interpreter','latex');
            s = size(wdpm4,1);n = 100/s;ytickstep = 5;
            xlabel('Node degree percentile (0 = all, 90 = hubs)');
            xticks(1:20);xticklabels(strsplit(num2str(5:5:100)));
            yticks(ytickstep:ytickstep:s);yticklabels(strsplit(num2str(s - ytickstep:-ytickstep:ytickstep)));
            ylabel('Percentile of eigenvectors (0 = all, 99 = last 60)');
            obj.save_fig(f,name);
            %%%%%
        end
        %%
        function obj = plot_eigenvectors5(obj)
            %%% fig18-1
            if ~isequal(sum(sum(abs(imag(obj.eigenvalues_ana)))),0)
                disp('there are complex ana eigenvalues');
            end
            if ~isequal(sum(sum(abs(imag(obj.eigenvalues_asy)))),0)
                disp('there are complex asy eigenvalues');
            end
            if ~isequal(sum(sum(abs(imag(obj.eigenvectors_ana)))),0)
                disp('there are complex ana eigenvectors');
            end
            if ~isequal(sum(sum(abs(imag(obj.eigenvectors_asy)))),0)
                disp('there are complex asy eigenvectors');
            end
            [B,I] = sort(obj.degree_vector_weighted);
            mydata = load(fullfile(obj.resultsPath,'obj_properties','eigenvectors_ana_ordered_nodes'),'var');
            tmp = mydata.var;evoana = tmp(end:-1:1,:);
            mydata = load(fullfile(obj.resultsPath,'obj_properties','eigenvectors_asy'),'var');
            tmp = mydata.var(I,obj.permutation_eigvec_ana2asy);evoasy = tmp(end:-1:1,:);
            evos = {evoana,evoasy};suf2 = {'-ana','-asy'};
            obj.eigenvectors_ana = General.load_var(fullfile(obj.resultsPath,'obj_properties/','eigenvectors_ana.mat'));
            obj.eigenvectors_asy_permuted = General.load_var(fullfile(obj.resultsPath,'obj_properties/','eigenvectors_asy_permuted.mat'));
            obj.eigenvalues_ana = General.load_var(fullfile(obj.resultsPath,'obj_properties/','eigenvalues_ana.mat'));
            for i2 = 1:length(evos)
                evo = evos{i2};
                [ind_bins_var,binned_vals] = EngineClass.set_bins_percentiles(1,obj.degree_vector_weighted);
                w = obj.degree_vector_weighted.*ones(size(obj.eigenvectors_ana));
                fs = {@(x) (1*x),@real,@imag};suf1 = {'-abs','-real','-imag'};
                desc1 = {'\mid', '\mid Re(', '\mid Im('};desc2 = {'\mid', ')\mid', ')\mid'};
                for i1 = 1:length(fs)
                    name = ['fig18-1' suf1{i1} suf2{i2}];
                    f = figure('Name',name,'NumberTitle','off');
                    figdesc = ['Eigenvectors mass distribution $ ' desc1{i1} ' v ' desc2{i1} ' $'];
                    image(abs(fs{i1}(evo)),'CDatamapping','scaled');
                    colorbar;set(gca,'ColorScale','log')
                    title({[name ' ' obj.scenarioName];obj.desc;figdesc},'interpreter','latex');
                    yticks(1000:1000:6000);yticklabels(strsplit(num2str(binned_vals(end-1000:-1000:1000)')));
                    xlabel('$v^{(i)}$','Interpreter','latex');ylabel('$k_i$','Interpreter','latex');
                    obj.save_fig(f,name);
                end
            end
            %%% fig18-2
            %%% fig18-3
            w2 = obj.degree_vector_weighted/norm(obj.degree_vector_weighted) .* ones(size(obj.eigenvectors_ana));
            ws = {w,w2};
            x1 = 1:size(obj.eigenvectors_ana,2);
            x2 = -1./obj.eigenvalues_ana;x3 = real(obj.eigenvalues_ana);
            x2asy = -1./obj.eigenvalues_asy_permuted;
            x = {x1,x2,x3}; xasy = {x1,x2asy};suf = {'','-2','-3'};
            xlabs = {'$v^{(i)}$','$Re(-1/\lambda_i)$','$Re(\lambda_i)$'};
            midstr = {'','\mid'};
            
            suf2 = {'','-1'}; suf3 = {'','-absv'};
            descsuf = {'','-normalized'}; fs = {@(x) (1*x), @abs};
            for i3 = 1:length(fs)
                f = fs{i3};g = fs{end - i3 + 1};
                ylabs = {[ '$' midstr{end-i3+1} midstr{i3} ' v ' midstr{i3} ' \cdot k ' midstr{end-i3+1} '$' ],...
                    ['$' midstr{end-i3+1} midstr{i3} ' v ' midstr{i3} ' \cdot \hat{k}' midstr{end-i3+1} '$']};
                for i2 = 2:length(ws)
                    ws1 = g(dot(f(obj.eigenvectors_ana),ws{i2}));
                    ws2 = g(dot(f(obj.eigenvectors_asy_permuted),ws{i2}));
                    figdesc = ['Eigenvectors weighted sum' descsuf{i2}];
                    for i1 = 2:length(x)
                        name = ['fig18-3' suf{i1} suf2{i2} suf3{i3}];
                        f = figure('Name',name,'NumberTitle','off');
                        plot(x{i1},ws1,'-o');hold on;
                        plot(x{i1},ws2,'o-');
                        %                     plot(xasy{i1},ws2,'o-');
                        legend('Real Eigenvectors','Theoretical Eigenvectors');
                        xlabel(xlabs{i1},'Interpreter','latex');
                        ylabel(ylabs{i2},'Interpreter','latex');
                        title({[name ' ' obj.scenarioName];obj.desc;figdesc},'interpreter','latex');
                        obj.save_fig(f,name);
                    end
                end
            end
            %%% fig18-4
            ws1 = abs(dot(obj.eigenvectors_ana,w2));
            [maxana,imaxana] = max(ws1); [minana,iminana] = min(ws1);
            ws2 = abs(dot(obj.eigenvectors_asy_permuted,w2));
            [maxasy,imaxasy] = max(ws2); [minasy,iminasy] = min(ws2);
            [B,I] = sort(obj.degree_vector_weighted);
            parvecana = obj.eigenvectors_ana(I,imaxana);
            parvecanabinned = obj.eigenvectors_ana_binned_k(:,imaxana);
            perpvecana = obj.eigenvectors_ana(I,iminana);
            perpvecanabinned = obj.eigenvectors_ana_binned_k(:,iminana);
            parvecasy = obj.eigenvectors_asy_permuted(I,imaxasy);
            parvecasybinned = obj.eigenvectors_asy_permuted_binned_k(:,imaxasy);
            perpvecasy = obj.eigenvectors_asy_permuted(I,iminasy);
            perpvecasybinned = obj.eigenvectors_asy_permuted_binned_k(:,iminasy);
            y = [parvecana,perpvecana,parvecasy,perpvecasy];
            ybinned = [parvecanabinned,perpvecanabinned,parvecasybinned,perpvecasybinned];
            xs = {B, obj.kbinned};ys = {y,ybinned}; suf = {'','-binned'};
            for i1=1:length(xs)
                name = ['fig18-4' suf{i1}];
                figdesc = ['Eigenvector node values vs. Degree' suf{i1}];
                f = figure('Name',name,'NumberTitle','off');
                loglog(xs{i1},abs(ys{i1}),'*-');
                legend(['v_{re,par}, ind = ' num2str(imaxana) ', dot = ' num2str(maxana)],...
                    ['v_{re,perp}, ind = ' num2str(iminana) ', dot = ' num2str(minana)],...
                    ['v_{th,par}, ind = ' num2str(imaxasy) ', dot = ' num2str(maxasy)],...
                    ['v_{th,perp}, ind = ' num2str(iminasy) ', dot = ' num2str(minasy)]);
                xlabel('$k_i$','Interpreter','latex');ylabel('$ \mid v_i \mid $ ($i^{th}$ element of v, sorted by degree)','Interpreter','latex');
                title({[name ' ' obj.scenarioName];obj.desc;figdesc},'interpreter','latex');
                obj.save_fig(f,name);
            end
        end
        %%
        function obj = plot_jacobian2(obj)
            Wij_ana = General.load_var(fullfile(obj.resultsPath,'obj_properties','Wij_ana'));
            Wij_asy = General.load_var(fullfile(obj.resultsPath,'obj_properties','Wij_asy'));
            Dii_ana = General.load_var(fullfile(obj.resultsPath,'obj_properties','Dii_ana'));
            Dii_asy = General.load_var(fullfile(obj.resultsPath,'obj_properties','Dii_asy'));
            Jana = Wij_ana + diag(Dii_ana);
            Jasy = Wij_asy + diag(Dii_asy);
            jacobians = {Jana,Jasy};
            namepre = 'fig15-';
            path = obj.resultsPath;
            title = {'Analytic Jacobian','Asymptotic Jacobian'};
            xlabel = '';
            ylabel = '';
            n = size(jacobians,2);
            for i1 = 1:n
                name = [namepre num2str(i1)];
                EngineClass.plot_image_static(jacobians{i1},name,name,path,title{i1},xlabel,ylabel);
            end
        end
        function obj = plot_jacobian3(obj)
            %% fig16-1*
            mydata=load(fullfile(obj.resultsPath,'obj_properties','J_ana_degree_vector_weighted'));
            J_ana_degree_vector_weighted = mydata.var;
            mydata=load(fullfile(obj.resultsPath,'obj_properties','J_asy_degree_vector_weighted'));
            J_asy_degree_vector_weighted = mydata.var;
            name = 'fig16-1';
            figdesc = 'Jacobian node degree vs. node degree';
            f = figure('Name',name,'NumberTitle','off');
            semilogx(obj.degree_vector_weighted,J_ana_degree_vector_weighted,'.');
            legendStr{1} = 'Real Jacobian';
            hold on;
            semilogx(obj.degree_vector_weighted,J_asy_degree_vector_weighted,'o');
            legendStr{2} = 'Theoretical Jacobian';
            ylabel('Jacobian Node Degree');
            xlabel('Adjacency Matrix Node Degree');
            title({[name ' ' obj.scenarioName];obj.desc;figdesc});
            legend(legendStr);
            obj.save_fig(f,name);
            %% fig16-2*
            mydata=load(fullfile(obj.resultsPath,'obj_properties','J_ana_k_inn'));
            J_ana_k_inn = mydata.var;
            mydata=load(fullfile(obj.resultsPath,'obj_properties','J_asy_k_inn'));
            J_asy_k_inn = mydata.var;
            name = 'fig16-2';
            figdesc = 'Jacobian vs. adjacency matrix nearest neighbor degree';
            f = figure('Name',name,'NumberTitle','off');
            loglog(obj.ki_nn,J_ana_k_inn,'.');
            legendStr{1} = 'Real Jacobian';
            hold on;
            loglog(obj.ki_nn,J_asy_k_inn,'o');
            legendStr{2} = 'Theoretical Jacobian';
            ylabel('Jacobian Nearest Neighbor degree k_{i,nn}');
            xlabel('Adjacency Matrix Nearest Neighbor Degree');
            title({[name ' ' obj.scenarioName];obj.desc;figdesc});
            legend(legendStr);
            obj.save_fig(f,name);
        end
        %% fig19* **fig20* fig21*
        function obj = plot_localization(obj)
            namepre = 'fig19-';
            figdesc = {'Perturbation Weighted Std',...
                'Perturbation Weighted Mean',...
                'Perturbation node max mass'};
            
            ylabs = {'$\hat{\sigma}_k$','$\hat{\mu}_k$','$\mid pert_{i,max}\mid /\sum_i \mid pert_i \mid$';...
                '$\hat{\sigma}_d$','$\hat{\mu}_d$',''};
            numfigs = length(figdesc);namesuf1 = {'-k','-d'};
            for i2 = 1:size(ylabs,1)
                kd = {'k','d'};
                figdesc2 = {['$\hat{\sigma} = \left[(\sum_i \hat{pert_i} (' kd{i2} '_i - \hat{\mu})^2)/S\right]^{1/2}$, $\hat{pert_i} = \mid pert_i \mid/S$, $S = \sum_i \mid pert_i \mid$'],...
                    ['$\hat{\mu} = \sum_i \mid pert_i \mid ' kd{i2} '_i / \sum_i \mid pert_i \mid$'],...
                    ''};
                for i1 = 1:size(ylabs,2)
                    if ~isempty(ylabs{i2,i1})
                        name = [namepre num2str(i1) namesuf1{i2}];
                        f{i2,i1} = figure('Name',name,'NumberTitle','off');
                        ax{i2,i1} = gca;hold on;
                        xlabel('t');ylabel(ylabs{i2,i1},'interpreter','latex');
                        title({[name ' ' obj.scenarioName];obj.desc;figdesc{i1};figdesc2{i1}},'interpreter','latex');
                    end
                end
            end
            step = 1;
            foldernamestrs = {'eigvec_pert_max_hub','eigvec_pert_min_hub'};
            desc1 = {'pert is real eigvec with most mass at biggest hub ($v^{(2)}$)', ...
                'pert is theoretical eigvec with most mass at biggest hub ($v^{(2)}$)',...
                'pert is all mass at biggest hub',...
                'pert is real eigvec with least mass at any node ($v^{(1181)}$)',...
                'pert is first real eigvec'};
            filenamesufstrs = {'ana','asy','hub'};
            fig20ind = 1;fig21ind = 1;
            legendInd = 1;legendStr = {};
            mypathst = {fullfile(obj.resultsPath,'obj_properties','sol_t_ana_v2'),...
                fullfile(obj.resultsPath,'obj_properties','sol_t_asy_v2'),...
                fullfile(obj.resultsPath,'obj_properties','eigvec_pert_max_hub','sol_t_hub'),...
                fullfile(obj.resultsPath,'obj_properties','eigvec_pert_min_hub','sol_t_ana'),...
                fullfile(obj.resultsPath,'obj_properties','sol_t_ana_v1')
                };
            mypathsx = {fullfile(obj.resultsPath,'obj_properties','sol_x_ana_v2'),...
                fullfile(obj.resultsPath,'obj_properties','sol_x_asy_v2'),...
                fullfile(obj.resultsPath,'obj_properties','eigvec_pert_max_hub','sol_x_hub'),...
                fullfile(obj.resultsPath,'obj_properties','eigvec_pert_min_hub','sol_x_ana'),...
                fullfile(obj.resultsPath,'obj_properties','sol_x_ana_v1')
                };
            mydata = load(fullfile('networks',obj.networkName, 'shortest_distances'));
            shortest_distances = mydata.var;clear mydata;
            descInd = length(desc1);ylabs2 = {'i (6000 = biggest hub)','i (6000 = farthest node from initial pert mass concentration'};
            for i2 = 1:length(ylabs2)
                for i4 = 1:length(mypathst)
                    try
                        mydata = load(mypathst{i4});sol_t = mydata.var(1:step:end);
                        mydata = load(mypathsx{i4});sol_x = mydata.var(1:step:end,:);
                        clear mydata;
                        num_quants = 1;
                        pert_x = (obj.steady_state' - sol_x')';clear sol_x;
                        if i2 == 1 % weiging by degree
                            weights = obj.degree_vector_weighted;
                        elseif i2 == 2 %weighing by distance from initial mass concentration node
                            pert_x0 = pert_x(1,:);
                            [~,i] = max(abs(pert_x0));
                            weights = shortest_distances(i,:)';
                        end
                        [means, vars, pert_x_ordered_quants] = EngineClass.compute_wmeans_wvars(abs(pert_x'), weights);
                        pert_x_ordered_quants = pert_x_ordered_quants';
                        %                         [bdegree,Idegree] = sort(obj.degree_vector_weighted);
                        %                         sd = shortest_distances(i,:)';
                        %                         [bdistance,Idistance] = sort(sd);
                        %                         spread_space = {bdegree,bdistance};
                        %                         spread_space_ind = {Idegree,Idistance};
                        %                         B = spread_space{i2};I = spread_space_ind{i2};
                        %                         pert_x_ordered = pert_x(:,I);clear pert_x;
                        %                         sum_pert_x_ordered = sum(abs(pert_x_ordered),2);
                        %                         pert_x_ordered_proportions = abs(pert_x_ordered)./sum_pert_x_ordered;clear pert_x_ordered sum_pert_x_ordered;
                        %                         pert_x_ordered_quants = num_quants*pert_x_ordered_proportions;clear pert_x_ordered_proportions;
                        %
                        %                         B_mat = B.*ones(size(pert_x_ordered_quants'));
                        %                         means = dot(B_mat,pert_x_ordered_quants')/num_quants;clear B_mat;
                        %                         prevar1 = B - means; %B is the ordered degree vector, should be a column,
                        %                         % means should be a row vector
                        %                         prevar2 = prevar1.^2; clear prevar1;
                        %                         prevar3 = pert_x_ordered_quants' .* prevar2; clear prevar2
                        %                         prevar4 = sum(prevar3,1); clear prevar3;
                        %                         var = sqrt(prevar4 / num_quants); clear prevar4;
                        plot(ax{i2,1},sol_t,vars);
                        plot(ax{i2,2},sol_t,means);
                        if ~isempty(ylabs{i2,3})
                            plot(ax{i2,3},sol_t,max(pert_x_ordered_quants,[],2));
                        end
                        descInd1 = mod(descInd,length(desc1))+1;
                        if legendInd <= length(desc1)
                            legendStr{legendInd} = desc1{descInd1};legendInd = legendInd+1;
                        end
                        %%% fig20-*
                        name = ['fig20-' num2str(i4) namesuf1{i2}];
                        f20 = figure('Name',name,'NumberTitle','off');
                        xlabel(ylabs2{i2});x = 1:obj.N;
                        semilogy(x,pert_x_ordered_quants(1,:),'*-');hold on;
                        tind = find(sol_t >= 46, 1);
                        semilogy(x,pert_x_ordered_quants(tind,:),'*-');
                        semilogy(x,pert_x_ordered_quants(end,:),'*-');
                        title({[name ' ' obj.scenarioName];obj.desc;desc1{descInd1}},'interpreter','latex');
                        fig20ind = fig20ind+1;xlabel('i');
                        ylabel('$\mid pert_i\mid/\sum_i \mid pert_i \mid$','Interpreter','latex');
                        legend('t = 0',['t = ' num2str(sol_t(tind))],'t = end');obj.save_fig(f20,name);
                        %%% fig21-*
                        name = ['fig21-' num2str(i4) namesuf1{i2}];
                        f21 = figure('Name',name,'NumberTitle','off');
                        image(abs(pert_x_ordered_quants(:,end:-1:1)'),'CDatamapping','scaled');
                        colorbar;set(gca,'ColorScale','log');
                        inds = 100:100:length(sol_t);
                        xticks(inds);xticklabels(strsplit(num2str(round(sol_t(inds)'))));
                        n = size(pert_x_ordered_quants,2);
                        inds = find(mod(1:n,500)==1); yticks(inds);
                        inds = find(mod(1:n,500)==0);inds=inds(end:-1:1);yticklabels(strsplit(num2str(inds)));
                        xlabel('t');ylabel(ylabs2{i2});
                        title({[name ' ' obj.scenarioName];obj.desc;desc1{descInd1};'$\mid pert_i\mid/\sum_i \mid pert_i \mid$'},'interpreter','latex');
                        descInd=descInd+1;obj.save_fig(f21,name);fig21ind = fig21ind+1;
                        
                    catch exception
                        disp(exception.message);
                        descInd=descInd+1;
                    end
                end
            end
            for i2 = 1:size(ylabs,1)
                for i1 = 1:size(ylabs,2)
                    if ~isempty(ylabs{i2,i1})
                        legend(ax{i2,i1},legendStr,'interpreter','latex');
                        name = [namepre num2str(i1) namesuf1{i2}];obj.save_fig(f{i2,i1},name);
                    end
                end
            end
        end
        %% fig 22*
        function obj = plot_localization2(obj)
            [ids_sorted_by_deg,B] = obj.get_ids_sorted_by_degs();
            ids = [ids_sorted_by_deg(1),ids_sorted_by_deg(end)];
            networkPath = fullfile('networks',obj.networkName);
            num_nodes = 100;
            num_times = 3;
            k = General.load_var(fullfile(networkPath,'degree_vector'));
            A = General.load_var(fullfile(networkPath,'adjacency_matrix'));
            dist = General.load_var(fullfile(networkPath,'shortest_distances'));
            %             [~,start_node] = max(k);
            
            
            %             sol_path_t = fullfile(obj.resultsPath,'obj_properties','eigvec_pert_max_hub','sol_t_hub');
            %             sol_path_x = fullfile(obj.resultsPath,'obj_properties','eigvec_pert_max_hub','sol_x_hub');
            %             sol_t = General.load_var(sol_path_t);
            %             sol_x = General.load_var(sol_path_x);
            %             ss = obj.steady_state;
            %             pert_x = (sol_x' - ss')';
            %             [means, vars, pert_x_ordered_quants] = EngineClass.compute_wmeans_wvars(abs(pert_x'), dist_start_node');
            %             TF = find(islocalmax(vars) | islocalmin(vars));
            %             num_rows = size(pert_x,1);
            %             step = floor(num_rows/(num_times-1));
            % %             pert_inds = 1:step:num_rows;
            %             pert_inds = [1 floor(TF(1)/2) TF];
            %             if abs(vars(TF(end)) - vars(end)) > .1
            %                 pert_inds = [pert_inds num_rows];
            %             end
            %             l = length(pert_inds);
            %             num_tiles_rows = floor(sqrt(l));
            %             num_tiles_cols = ceil(l/num_tiles_rows);
            
            num_tiles_rows = 2;num_tiles_cols = 2;pert_norm_hat_snapshots = [1,.5,.1,.02];
            threshold = 1e-5;
            we = {'none','inverse'};inds_top_num_nodes = [];
            for i2 = 1:length(ids)
                id = ids(i2);
                [~,~,pert_x,~,pert_x_norm_hat] = obj.get_pert(id);
                start_node = id;
                for i = pert_norm_hat_snapshots
                    ind = find(pert_x_norm_hat<=i,1);
                    pert_x1 = abs(pert_x(ind,:));
                    [B,inds_cur] = sort(pert_x1,'descend');
                    stop = find(B<threshold, 1);
                    if isempty(stop)
                        stop = num_nodes+1;
                    end
                    stop = min(stop-1,num_nodes);
                    inds_top_num_nodes = union(inds_top_num_nodes,inds_cur(1:stop));
                end
            end
            for i2 = 1:length(ids)
                id = ids(i2);
                [~,~,pert_x,~,pert_x_norm_hat] = obj.get_pert(id);
                name = ['fig22a-id-' num2str(id)];
                f = figure('Name',name,'NumberTitle','off');
                t = tiledlayout(num_tiles_rows,num_tiles_cols);
                %                 dist_start_node = dist(start_node,:);
                %                 [dist_start_node_sorted,inds] = sort(dist_start_node);
                %                 inds_top_num_nodes = inds(1:num_nodes);
                A_num_nodes = A(inds_top_num_nodes,inds_top_num_nodes);
                k_num_nodes = k(inds_top_num_nodes);
                k_max = max(k_num_nodes);
                k_min = min(max(k_num_nodes)-1,0);
                k_range = k_max-k_min;
                G = graph(A_num_nodes);
                max_marker_size = 16;
                desc1 = ['Perturbation mass concentration diffusion, $k_i = ' num2str(obj.degree_vector_weighted(id)) '$'];
                
                for i1 = pert_norm_hat_snapshots
                    nexttile;
                    p=plot(G,'LineWidth',.1,'LineStyle','-','Marker','o','layout','force','WeightEffect','direct','NodeLabel',{});%'Center',find(inds_top_num_nodes==id,1), 'WeightEffect','inverse'
                    %                         p=plot(G,'LineWidth',.1,'LineStyle','-','Marker','o','layout','subspace','NodeLabel',{});
                    for i = 1:length(inds_top_num_nodes)
                        k_i = k_num_nodes(i);
                        markerSize = (k_i-k_min)/k_range*max_marker_size+4;
                        highlight(p,i,'MarkerSize',markerSize);
                    end
                    ind = find(pert_x_norm_hat<=i1,1);
                    pert_x1 = abs(pert_x(ind,:));
                    pert_x1_num_nodes = pert_x1(1,inds_top_num_nodes);
                    norm_pert_x1_num_nodes = norm(pert_x1_num_nodes);
                    pert_x1_num_nodes = pert_x1_num_nodes/norm_pert_x1_num_nodes;
                    pert_x1_max = max(pert_x1_num_nodes);
                    pert_x1_min = min(pert_x1_num_nodes);
                    pert_x1_range = pert_x1_max - pert_x1_min;
                    pert_x1_min = pert_x1_min - threshold^2*pert_x1_range;
                    pert_x1_range = pert_x1_max - pert_x1_min;
                    
                    for i = 1:length(inds_top_num_nodes)
                        pert_x1_i = pert_x1_num_nodes(i);
                        marker_color1 = (pert_x1_i - pert_x1_min)/pert_x1_range;
                        marker_color = 1/(1 - log10(marker_color1));
                        %                         marker_color = log10(pert_x1_i * pert_x1_max / pert_x1_min) / (log10(pert_x1_max) - log10(pert_x1_min));
                        highlight(p,i,'NodeColor',[1 1-marker_color 1-marker_color]);
                    end
                    
                    title(['$\mid pert \mid = ' num2str(i1) '*\mid pert_0 \mid$'],'Interpreter','latex');
                end
                title(t,{[name ' ' obj.scenarioName];obj.desc;desc1},'interpreter','latex');
                obj.save_fig(f,name);
            end
        end
        %% fig23*
        function obj = plot_localization3(obj)
            powers = [1 3 10];
            name = 'fig23-1-t'; fname{1} = name;%pert_norm vs t
            f{1} = figure('Name',name,'NumberTitle','off');
            hax{1} = axes;hold on;xlabel('t','Interpreter','latex');ylabel('$\mid pert \mid$','Interpreter','latex');
            desc1 = 'Perturbation norm vs t';
            title({[name ' ' obj.scenarioName];obj.desc;desc1},'interpreter','latex');
            name = 'fig23-1-tau'; fname{2} = name;%pert_norm vs tau
            f{2} = figure('Name',name,'NumberTitle','off');
            hax{2} = axes;hold on;xlabel('$\tau$','Interpreter','latex');ylabel('$\mid pert \mid$','Interpreter','latex');
            desc1 = 'Perturbation norm vs $\tau$ ($\tau$ = half-life)';
            title({[name ' ' obj.scenarioName];obj.desc;desc1},'interpreter','latex');
            name = 'fig23-2-t'; fname{3} = name;%pert_pert0_dot_normed vs t
            f{3} = figure('Name',name,'NumberTitle','off');
            hax{3} = axes;hold on;xlabel('t','Interpreter','latex');ylabel('$\widehat{pert} \cdot \widehat{pert_0}$','Interpreter','latex');
            desc1 = 'Perturbation shape compared to initial perturbation vs t';
            title({[name ' ' obj.scenarioName];obj.desc;desc1},'interpreter','latex');
            name = 'fig23-2-tau'; fname{4} = name; %pert_pert0_dot_normed vs tau
            f{4} = figure('Name',name,'NumberTitle','off');
            hax{4} = axes;hold on;xlabel('$\tau$','Interpreter','latex');ylabel('$\widehat{pert} \cdot \widehat{pert_0}$','Interpreter','latex');
            desc1 = 'Perturbation shape compared to initial perturbation vs $\tau$ ($\tau$ = half-life)';
            title({[name ' ' obj.scenarioName];obj.desc;desc1},'interpreter','latex');
            name = 'fig23-3'; fname{5} = name; %half-lives vs t
            f{5} = figure('Name',name,'NumberTitle','off');
            hax{5} = axes; hold on; xlabel('t','Interpreter','latex');ylabel('Half-life','Interpreter','latex');
            desc1 = 'Perturbation half-life vs t';
            title({[name '' obj.scenarioName];obj.desc;desc1},'interpreter','latex');
            legendInd = 1;legendStr = cell(1,length(powers));
            for i = powers
                fnamet = ['sol_t_pert_k_power_' num2str(i)];path = fullfile(obj.resultsPath,'obj_properties','sol_pert_k_power',fnamet);
                sol_t = General.load_var(path);
                fnamex = ['sol_x_pert_k_power_' num2str(i)];path = fullfile(obj.resultsPath,'obj_properties','sol_pert_k_power',fnamex);
                sol_x = General.load_var(path);
                [pert_norm, pert_pert0_dot_normed, half_life, sol_tau,ind1,ind2] = EngineClass.compute_pert_props(sol_t,sol_x,obj.steady_state);
                plot(hax{1},sol_t,pert_norm);legend;
                plot(hax{2},sol_tau,pert_norm);legend;
                plot(hax{3},sol_t,pert_pert0_dot_normed);legend;
                plot(hax{4},sol_tau,pert_pert0_dot_normed);legend;
                plot(hax{5},sol_t(ind1:ind2),half_life);
                legendStr{legendInd} = ['$pert_0 \sim k_i^{' num2str(i) '}$'];legendInd = legendInd+1;
            end
            
            for i = 1:length(f)
                legend(hax{i},legendStr,'interpreter','latex');
                EngineClass.save_fig_static(f{i},fname{i},obj.resultsPath);
            end
        end
        %% Helper function
        function [ids_sorted_by_deg,B] = get_ids_sorted_by_degs(obj)
            p = fullfile(obj.resultsPath,'obj_properties','single_node_pert_sol');
            files = dir(p);
            %names = cell(1,length(files)-2);
            ids = zeros(1,length(files)/2 - 1);ind = 1;
            for i = 1:length(files)
                file = files(i);
                if file.name(1) ~= '.' && file.name(5) == 't'
                    %                     names{i} = files.name;
                    ind1 = find(file.name == '.');
                    ids(ind) = str2double(file.name(7:ind1-1));
                    ind = ind+1;
                end
            end
            ids = sort(ids);
            degs = obj.degree_vector_weighted(ids);
            [B,I] = sort(degs,'descend');
            ids_sorted_by_deg = ids(I);
        end
        function [sol_t,sol_x,pert_x,pert_x_0,pert_x_norm_hat] = get_pert(obj,id)
            p = fullfile(obj.resultsPath,'obj_properties','single_node_pert_sol');
            sol_t = General.load_var(fullfile(p,['sol_t_' num2str(id)]));
            sol_x = General.load_var(fullfile(p,['sol_x_' num2str(id)]));
            pert_x = (sol_x' - obj.steady_state')';
            pert_x_0 = pert_x(1,:);
            pert_x_norm = vecnorm(pert_x');
            pert_x_norm_hat = pert_x_norm./norm(pert_x_0);
            ind = find(pert_x_norm_hat < .02,1);
            sol_t = sol_t(1:ind);
            sol_x = sol_x(1:ind,:);
            pert_x = pert_x(1:ind,:);
            pert_x_norm_hat = pert_x_norm_hat(1:ind);
        end
        %% fig24*
        function obj = plot_single_node_pert(obj)
            [ids_sorted_by_deg,B] = obj.get_ids_sorted_by_degs();
            %%%
            name = 'fig24a';fname{1} = name;%pert_norm/pert0_norm vs t
            f{1} = figure('Name',name,'NumberTitle','off');
            hax{1} = axes;hold on;xlabel('t','Interpreter','latex');ylabel('$\mid pert \mid / \mid pert_0 \mid$','Interpreter','latex');
            desc1 = 'Perturbation norm (relative to $pert_0$) vs $t$, $pert_0 = \delta(i)$';
            title({[name ' ' obj.scenarioName];obj.desc;desc1},'interpreter','latex');
            
            name = 'fig24b';fname{2} = name;%tau vs k_i
            f{2} = figure('Name',name,'NumberTitle','off');
            hax{2} = axes;hold on;xlabel('$k_i$','Interpreter','latex');ylabel('$\tau$','Interpreter','latex');
            desc1 = '$\tau = $ time when $\mid pert \mid = .02*\mid pert_0 \mid$';
            title({[name ' ' obj.scenarioName];obj.desc;desc1},'interpreter','latex');
            
            taus = zeros(length(ids_sorted_by_deg,1));legendInd = 1;legendStr = cell(1,length(ids));
            for i = 1:length(ids_sorted_by_deg)
                id = ids_sorted_by_deg(i);
                [sol_t,~,~,~,pert_x_norm_hat] = obj.get_pert(id);
                taus(i) = sol_t(ind);
                plot(hax{1},sol_t,pert_x_norm_hat);legendStr{legendInd} = ['$k_i = ' num2str(B(i)) '$'];
                legendInd = legendInd+1;
            end
            legend(hax{1},legendStr,'interpreter','latex');
            plot(hax{2},B,taus,'o');
            
            for i = 1:length(f)
                EngineClass.save_fig_static(f{i},fname{i},obj.resultsPath);
            end
        end
        %% fig25*
        function obj = plot_eigvecs_power_law(obj)
            percentiles = [100, 50, 10, 5, 1];
            %%%
            name{1} = 'fig25a';f{1} = figure('Name',name{1},'NumberTitle','off');
            desc1 = 'Linear regression $\theta$ values, $v_i = C(k)k_i^\theta$';
            hax{1} = axes; hold on; xlabel('$j$ ($j^{th}$ eigenvector)','Interpreter','latex'); ylabel('$\theta$','Interpreter','latex');
            title({[name{1} ' ' obj.scenarioName];obj.desc;desc1},'interpreter','latex');
            
            name{2} = 'fig25b';f{2} = figure('Name',name{2},'NumberTitle','off');
            desc1 = 'Linear regression $R^2$ values, $v_i = C(k)k_i^\theta$';
            hax{2} = axes; hold on; xlabel('$j$ ($j^{th}$ eigenvector)','Interpreter','latex'); ylabel('$R^2$','Interpreter','latex');
            title({[name{2} ' ' obj.scenarioName];obj.desc;desc1},'interpreter','latex');
            legendStr = {};
            for percentile = percentiles
                cur_path = fullfile(obj.resultsPath,'obj_properties','eigvec_power_law',[num2str(percentile) 'percent']);
                ps = General.load_var(fullfile(cur_path,'ps'));
                plot(hax{1},ps(:,1));
                rsqs = General.load_var(fullfile(cur_path,'rsqs'));
                plot(hax{2},rsqs);
                legendStr{end+1} = ['Top ' num2str(percentile) '% of nodes by degree'];
            end
            for i = 1:length(f)
                legend(hax{i},legendStr);
                General.save_fig(f{i},name{i},fullfile(obj.resultsPath,'figs'));
            end
        end
        %% fig26*
        function obj = plot_concentration(obj)
            num_nodes = [1 10 100 1000];
            netPath = fullfile('networks',obj.networkName,'random_samples');
            solPath = fullfile(obj.resultsPath,'obj_properties','random_sample_perts');
            name{1} = 'fig26a';f{1} = figure('Name',name{1},'NumberTitle','off');
            desc = 'Random Perturbation concentration vs time';hax{1} = axes;hold on;
            xlabel('t'); ylabel('C(t)');title({[name{1} ' ' obj.scenarioName];obj.desc;desc},'interpreter','latex');
            name{2} = 'fig26b';f{2} = figure('Name',name{2},'NumberTitle','off');
            desc = 'Random Perturbation conservation vs time';hax{2} = axes;hold on;
            xlabel('t'); ylabel('K(t)');title({[name{2} ' ' obj.scenarioName];obj.desc;desc},'interpreter','latex');
            name{3} = 'fig26c';f{3} = figure('Name',name{3},'NumberTitle','off');
            desc = 'Hub Perturbation concentration vs time';hax{3} = axes;hold on;
            xlabel('t');ylabel('C(t)');title({[name{3} ' ' obj.scenarioName];obj.desc;desc},'interpreter','latex');
            name{4} = 'fig26d';f{4} = figure('Name',name{4},'NumberTitle','off');
            desc = 'Hub Perturbation conservation vs time';hax{4} = axes;hold on;
            xlabel('t');ylabel('K(t)');title({[name{4} ' ' obj.scenarioName];obj.desc;desc},'interpreter','latex');
            name{5} = 'fig26e';f{5} = figure('Name',name{5},'NumberTitle','off');
            desc = 'Hubs Perturbation concentration vs time';hax{5} = axes;hold on;
            xlabel('t');ylabel('C(t)');title({[name{5} ' ' obj.scenarioName];obj.desc;desc},'interpreter','latex');
            name{6} = 'fig26f';f{6} = figure('Name',name{6},'NumberTitle','off');
            desc = 'Hubs Perturbation conservation vs time';hax{6} = axes;hold on;
            xlabel('t');ylabel('K(t)');title({[name{6} ' ' obj.scenarioName];obj.desc;desc},'interpreter','latex');
            legendStr = {};
            [~,J] = sort(obj.degree_vector_weighted,'descend');
            for n = num_nodes
                I = General.load_var(fullfile(netPath,['n_' num2str(n)]));
                pert_0 = General.load_var(fullfile(solPath,['n_' num2str(n)]));
                sol_t = General.load_var(fullfile(solPath,['sol_t_n_' num2str(n)]));
                sol_x = General.load_var(fullfile(solPath,['sol_x_n_' num2str(n)]));
                %                 pert_0 = sol_x(1,:);
                r = size(sol_x,1);concentration = zeros(r,1);conservation = zeros(r,1);
                for i = 1:r
                    pert_t = sol_x(i,:) - obj.steady_state;
                    concentration(i,1) = EngineClass.compute_concentration(pert_t,I);
                    conservation(i,1) = EngineClass.compute_conservation(pert_t,pert_0);
                end
                plot(hax{1},sol_t,concentration);
                plot(hax{2},sol_t,conservation);
                legendStr{end+1} = ['$\mid I \mid$ = ' num2str(n)];
            end
            for i = 1:2
                legend(hax{i},legendStr,'interpreter','latex');
                General.save_fig(f{i},name{i},fullfile(obj.resultsPath,'figs'));
            end
            [K,J] = sort(obj.degree_vector_weighted,'descend');
            solPath = fullfile(obj.resultsPath,'obj_properties','single_node_pert_sol');
            legendStr = {};ind=1;
            for j = J(1:4)'
                k = K(ind);
                sol_t = General.load_var(fullfile(solPath,['sol_t_' num2str(j)]));
                sol_x = General.load_var(fullfile(solPath,['sol_x_' num2str(j)]));
                pert_0 = sol_x(1,:) - obj.steady_state;
                r = size(sol_x,1);concentration = zeros(r,1);conservation = zeros(r,1);
                for i = 1:r
                    pert_t = sol_x(i,:) - obj.steady_state;
                    concentration(i,1) = EngineClass.compute_concentration(pert_t,j);
                    conservation(i,1) = EngineClass.compute_conservation(pert_t,pert_0);
                end
                plot(hax{3},sol_t,concentration);
                plot(hax{4},sol_t,conservation);
                legendStr{end+1} = ['Hub ' num2str(ind) ', i = ' num2str(j) ', k_i = ' num2str(k)];
                ind = ind+1;
            end
            for i = 3:4
                legend(hax{i},legendStr);
                General.save_fig(f{i},name{i},fullfile(obj.resultsPath,'figs'));
            end
            legendStr = {};ind=1;
            for j = 2:4
                nodes = J(1:j)';
                sol_t = General.load_var(fullfile(solPath,['sol_t_' EngineClass.array2str(nodes)]));
                sol_x = General.load_var(fullfile(solPath,['sol_x_' EngineClass.array2str(nodes)]));
                pert_0 = sol_x(1,:) - obj.steady_state;
                r = size(sol_x,1);concentration = zeros(r,1);conservation = zeros(r,1);
                for i = 1:r
                    pert_t = sol_x(i,:) - obj.steady_state;
                    concentration(i,1) = EngineClass.compute_concentration(pert_t,nodes);
                    conservation(i,1) = EngineClass.compute_conservation(pert_t,pert_0);
                end
                plot(hax{5},sol_t,concentration);
                plot(hax{6},sol_t,conservation);
                legendStr{end+1} = ['Nodes ' num2str(nodes)];
            end
            for i = 5:6
                legend(hax{i},legendStr);
                General.save_fig(f{i},name{i},fullfile(obj.resultsPath,'figs'));
            end
        end
        %% fig28*
        function plot_pert_amp_phase_vs_t(obj)
            %             foldernames = {'eigvec_pert_max_hub','eigvec_pert_min_hub',...
            %                 'eigvec_pert_min_hub_1','eigvec_power_law',...
            %                 'random_sample_perts','single_node_pert_sol',...
            %                 'sol_pert_k_power','solution_random_perts'};
            %             foldernames = {'solution_random_perts'};
            %             foldernames = {'eigvec_pert_max_hub','eigvec_pert_min_hub',...
            %                 'solution_random_perts','solution_random_perts_2'};
            foldernames = {'solution_random_perts','solution_random_perts_2','single_node_pert_sol','solution_eigvec_perts'};
            colors = {'b','r','k','m','g','c'};
            markers = {'o','s','^','*','+','x','v','<','>','.','d','p','h','_','|'};
            marksize = 6;
            nameroot = 'fig28';suf1s = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n'};
            funs = {@plot,@loglog};funssuf = {'','-loglog'};
            for i3 = 1:length(suf1s)
                suf1 = suf1s{i3};
                for i2 = 1:length(funs)
                    suf = funssuf{i2};
                    names{1,i3,i2} = [nameroot suf1 suf]; fname{1,i3,i2} = names{1,i3,i2};
                    f{1,i3,i2} = figure('Name',names{1,i3,i2},'NumberTitle','off');hax{1,i3,i2} = axes; %hold on;
                    if i3>3
                        for i1 = 1:length(obj.sys_half_life_amp)
                            names{i1,i3,i2} = [nameroot suf1 '-' num2str(i1) suf]; fname{i1,i3,i2} = names{i1,i3,i2};
                            f{i1,i3,i2} = figure('Name',names{i1,i3,i2},'NumberTitle','off');hax{i1,i3,i2} = axes;% hold on;
                        end
                    end
                end
            end
            legendStr = {};
            for j1 = 1:length(foldernames)
                folder = foldernames{j1};%folder = foldernames
                folderpath = fullfile(obj.resultsPath,'obj_properties/',folder);
                files = dir(folderpath);
                ind = 1;
                for i1 = 1:length(files)
                    file = files(i1);
                    name = file.name;
%                     ind1 = strfind(name,'_');
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
                        legendStr{end+1} = EngineClass.filename2legend(name(ind1+1:indsuf-1));
                        for i2 = 1:length(funs)
                            fun = funs{i2};
                            
                            
                            yvalues = General.load_var(fullfile(folderpath,name));
                            tvalues = General.load_var(fullfile(folderpath,[name(1:ind1-2) 't' name(ind1:indsuf-1)]));
                            fun(hax{1,1,i2},tvalues,yvalues,colors{j1},'LineStyle','-','Marker',markers{ind},'MarkerSize',marksize);hold(hax{1,1,i2},'on');
                            yvalues = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_norms_normed']));
                            fun(hax{1,2,i2},tvalues,yvalues,colors{j1},'LineStyle','-','Marker',markers{ind},'MarkerSize',marksize);hold(hax{1,2,i2},'on');
                            yvalues = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_thetas']));
                            fun(hax{1,3,i2},tvalues,rad2deg(abs(yvalues)),colors{j1},'LineStyle','-','Marker',markers{ind},'MarkerSize',marksize);hold(hax{1,3,i2},'on');
                            Q = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_Q']));
                            Q2 = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_Q2']));
                            Q3 = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_Q3']));
                            Q4 = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_Q4']));
                            Q5 = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_Q5']));
                            Q6 = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_Q6']));
                            Q7 = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_Q7']));
                            Q8 = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_Q8']));
                            Q9 = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_Q9']));
                            Q10 = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_Q10']));
                            tau_A = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_tau_A']));
                            for i1 = 1:length(tau_A)
                                tau = tau_A{i1};
                                fun(hax{i1,4,i2},Q,tau,colors{j1},'LineStyle','none','Marker',markers{ind},'MarkerSize',marksize);hold(hax{i1,4,i2},'on');
                                %                         tau_A_2 = General.load_var(fullfile(folderpath,[name(1:indsuf-1) '_tau_A_2']));
                                %                         plot(hax{5},Q,tau_A_2,colors{j1},'LineStyle','none','Marker',markers{ind},'MarkerSize',size);
                                
                                fun(hax{i1,6,i2},Q2,tau,colors{j1},'LineStyle','none','Marker',markers{ind},'MarkerSize',marksize);hold(hax{i1,6,i2},'on');
                                
                                fun(hax{i1,7,i2},Q3,tau,colors{j1},'LineStyle','none','Marker',markers{ind},'MarkerSize',marksize);hold(hax{i1,7,i2},'on');
                                
                                fun(hax{i1,8,i2},Q4,tau,colors{j1},'LineStyle','none','Marker',markers{ind},'MarkerSize',marksize);hold(hax{i1,8,i2},'on');
                                
                                fun(hax{i1,9,i2},Q5,tau,colors{j1},'LineStyle','none','Marker',markers{ind},'MarkerSize',marksize);hold(hax{i1,9,i2},'on');
                                
                                fun(hax{i1,10,i2},Q6,tau,colors{j1},'LineStyle','none','Marker',markers{ind},'MarkerSize',marksize);hold(hax{i1,10,i2},'on');
                                
                                fun(hax{i1,11,i2},Q7,tau,colors{j1},'LineStyle','none','Marker',markers{ind},'MarkerSize',marksize);hold(hax{i1,11,i2},'on');
                                
                                fun(hax{i1,12,i2},Q8,tau,colors{j1},'LineStyle','none','Marker',markers{ind},'MarkerSize',marksize);hold(hax{i1,12,i2},'on');
                                fun(hax{i1,13,i2},Q9,tau,colors{j1},'LineStyle','none','Marker',markers{ind},'MarkerSize',marksize);hold(hax{i1,13,i2},'on');
                                fun(hax{i1,14,i2},Q10,tau,colors{j1},'LineStyle','none','Marker',markers{ind},'MarkerSize',marksize);hold(hax{i1,14,i2},'on');
                            end
                        end
                        ind = ind+1;
                        if ind > length(markers)
                            ind = 1;
                            marksize = marksize + 4;
                        end
                        
                    end
                end
            end
            descs = {'Perturbation Amplitude vs t','Perturbation Amplitude normed vs t',...
                'Perturbation Phase [deg] vs t'};
            desc1 = ['$\tau_{A_' num2str(i1) '}$ vs Q'];
            ylabs = {'$\|\delta(t)\|$','$A(t)$','$\phi(t)$'};
            xlabs = {'$Q7 = \frac{\mathbf{\delta x(0)}\cdot k^{\alpha}}{\| \mathbf{\delta x(0)} \|_1 \max{(k)}}$',...
                '$Q8 = \frac{\mathbf{\delta x(0)}\cdot k^{\alpha}}{\| \mathbf{\delta x(0)} \|_1 }$',...
                '$Q9 = \frac{\mathbf{\delta x(0)}\cdot k^{-\mu} \cdot \mathbf{ss}}{\| \mathbf{\delta x(0)} \|_1 }$',...
                '$Q10 = \frac{\mathbf{\delta x(0)}\cdot k^{\alpha}}{\| \mathbf{\delta x(0)}\cdot k^{-\xi} \|_1 }$',};

            for i3 = 1:size(f,1)
                for i4 = 1:size(f,2)
                    if i4<4
                        xlab = 't';ylab = ylabs{i4};desc = descs{i4};
                    else
                        xlab = ['Q' num2str(i4-3)];ylab =['$\tau_{A_' num2str(i3) '}$'];
                        desc = [desc1 num2str(i4-3)];
                    end
                    if i4>10
                        xlab = xlabs{i4-7-3};
                    end
                    for i5 = 1:size(f,3)
                        if ~isempty(f{i3,i4,i5})
                            title(hax{i3,i4,i5},{[names{i3,i4,i5} ' ' obj.scenarioName];obj.desc;desc},'interpreter','latex');
                            xlabel(hax{i3,i4,i5},xlab,'Interpreter','latex');
                            ylabel(hax{i3,i4,i5},ylab,'Interpreter','latex');
                            legend(hax{i3,i4,i5},legendStr);
                            General.save_fig(f{i3,i4,i5},fname{i3,i4,i5},fullfile(obj.resultsPath,'figs'));
                        end
                    end
                end
            end            
        end
        %%
        %% fig29*
        function plot_Q_partial_sums(obj)
            x_sol_random_pert = General.load_var(fullfile(obj.resultsPath,'obj_properties','solution_random_perts','sol_x_random_perturbations_3'));
            legendStr{1} = 'Random Perturbation 3';
            x_sol_hub_low_Q = General.load_var(fullfile(obj.resultsPath,'obj_properties','single_node_pert_sol','sol_x_7'));
            legendStr{2} = 'Node 7 perturbation';
            ss = General.load_var(fullfile(obj.resultsPath,'obj_properties','steady_state'));
            k = obj.degree_vector_weighted;k_alpha = k.^(-obj.xi - obj.mu);k_alpha_normed = k_alpha/max(k_alpha);
            x_sols = {x_sol_random_pert,x_sol_hub_low_Q};
            partial_sums = cell(size(x_sols));
            for i1 = 1:length(x_sols)
                x = x_sols{i1};
                x = (x' - ss')';
                % Calc pert/ss, or dlogx (dx/x)
                x = (x'./ss')';
                x_init = x(1,:);
                x_init_normed = x_init/norm(x_init,1);
                partial_sum = zeros(length(k),1);
                prev = 0;
                for i2 = 1:length(k)
                    partial_sum(i2) = prev + (abs(x_init_normed(i2)) * k_alpha_normed(i2));
                    prev = partial_sum(i2);
                end
                partial_sums{i1} = partial_sum;
            end
            name = 'fig29a'; fname{1} = name;desc = 'Q7 Partial Sum Calculation';
            f{1} = figure('Name',name,'NumberTitle','off');hax{1} = axes; hold on;
            title({[name ' ' obj.scenarioName];obj.desc;desc},'interpreter','latex');
            xlabel('Node index $i$','Interpreter','latex');ylabel('$\sum_{j=0}^i \frac{\delta x_j(0)\cdot k_j}{\mid \mathbf{\delta x(0)} \max{(k)}}$','Interpreter','latex');
            
            for i1 = 1:length(partial_sums)
                plot(hax{1},partial_sums{i1},'-*');
            end
            legend(hax{1},legendStr);
            
            
        end
        %% fig30*
        function plot_network_dts(obj)
            fname = 'fig30a';
            sfactor = 1;
            network = obj.networkName;
            netpath = fullfile('networks/',network);
            k = General.load_var(fullfile(netpath,'degree_vector'));
            [k_sorted,I] = sort(k,'descend');
            A = General.load_var(fullfile(netpath,'adjacency_matrix'));
            f = openfig(fullfile(netpath,'figs/',"net09a.fig"));
            hax = f.Children;
            nodes = hax.Children;
            sol_x = General.load_var(fullfile(obj.resultsPath,'obj_properties',...
                'single_node_pert_sol','sol_x_1.mat'));
            tau_As = General.load_var(fullfile(obj.resultsPath,'obj_properties',...
                'single_node_pert_sol/','sol_x_1_tau_A.mat'));
            tau_A = tau_As{1};
            pert = (sol_x' - obj.steady_state')';
            pert0 = pert(1,:);
            R = General.load_var(fullfile(netpath,'net09-R'));
            thetas = General.load_var(fullfile(netpath,'net09-thetas'));
            phis = General.load_var(fullfile(netpath,'net09-phis'));
            
            f.Name = fname;
            ind = find(abs(pert0) >0,1,'first');
            node = nodes(end+1-ind);
            
            R_max = max(R);
            old_R = R(ind);
            theta = thetas(ind);phi=phis(ind);
            new_R =(1 + pert0(ind)*sfactor)*R_max;
            [node.XData,node.YData,node.ZData] = sph2cart(...
                theta,phi,new_R);
            
            A_ind = A(ind,:);
            num_links_to_draw = 5;num_links_drawn = 0;
            for i1 = 1:obj.N
                if A_ind(I(i1)) > 0
                    node2 = nodes(end+1-I(i1));
                    x = [node.XData, node2.XData];
                    y = [node.YData, node2.YData];
                    z = [node.ZData, node2.ZData];
                    plot3(x,y,z,'LineStyle','--',...
                        'Color',node2.Color,...
                        'LineWidth',0.1);
                    num_links_drawn = num_links_drawn+1;
                    if num_links_drawn > num_links_to_draw
                        break
                    end
                end
            end
            General.save_fig(f,fname,fullfile(obj.resultsPath,'figs'));
            
           % fig30b
           fname = 'fig30b';
           f.Name = fname;
           
           max_steps = 10;
           min_steps = 1;
           tau_A = floor(2*tau_A);
           if tau_A > 2
               num_steps = min(max_steps,tau_A);
           elseif tau_A < 2
               num_steps = max(min_steps,tau_A);
           else
               num_steps = 2;
           end
           R_step_size = (new_R - old_R)/num_steps;
           Rs = new_R:-R_step_size:old_R;
           [x,y,z] = sph2cart(theta,phi,Rs(2:end));
           orig_color = node.MarkerFaceColor;
           white_diff = [1 1 1] - orig_color;
           shades = 1:-1/(num_steps+1):1/(num_steps+1);
           white_diffs = white_diff'.*shades;
           new_colors = ([1;1;1] - white_diffs)';
           
           node.MarkerFaceColor = new_colors(end,:);
           node.MarkerEdgeColor = new_colors(end,:);
           for i1 = 1:num_steps
               plot3(x(i1),y(i1),z(i1),'o','MarkerSize',node.MarkerSize,...
                   'MarkerFaceColor',new_colors(end-i1,:),...
                   'MarkerEdgeColor',new_colors(end-i1,:));hold on;
           end
           General.save_fig(f,fname,fullfile(obj.resultsPath,'figs'));
        end
        %% fig31*
        function plot_network_dts_2(obj)
            %fig31a

            myfuncs = {@abs,@(x) log10(abs(x)*9+1)};figsufs = {'','-zlog'};
            basefigs = {'net09c.fig','net10a.fig'};
            sfactor = 1;
            network = obj.networkName;
            netpath = fullfile('networks/',network);            
%             sol_t = General.load_var(fullfile(obj.resultsPath,'obj_properties',...
%                 'single_node_pert_sol','sol_t_2_4_1_22_10.mat'));
%             sol_x = General.load_var(fullfile(obj.resultsPath,'obj_properties',...
%                 'single_node_pert_sol','sol_x_2_4_1_22_10.mat'));
%             sol_t = General.load_var(fullfile(obj.resultsPath,'obj_properties',...
%                 'single_node_pert_sol','sol_t_2_4_1_22_10_19_9_46_18_62.mat'));
%             sol_x = General.load_var(fullfile(obj.resultsPath,'obj_properties',...
%                 'single_node_pert_sol','sol_x_2_4_1_22_10_19_9_46_18_62.mat'));
%             sol_t = General.load_var(fullfile(obj.resultsPath,'obj_properties',...
%                 'solution_random_perts','sol_t_random_perturbations_1.mat'));
%             sol_x = General.load_var(fullfile(obj.resultsPath,'obj_properties',...
%                 'solution_random_perts','sol_x_random_perturbations_1.mat'));
            sol_t = General.load_var(fullfile(obj.resultsPath,'obj_properties',...
                'random_bin_perts','sol_t_pert1.mat'));
            sol_x = General.load_var(fullfile(obj.resultsPath,'obj_properties',...
                'random_bin_perts','sol_x_pert1.mat'));
            
            pert = (sol_x' - obj.steady_state')';max_pert = max(abs(pert),[],'all');
            pert = pert/max_pert;
            pert0 = pert(1,:);nodes_sub_inds = 1:size(pert0,2); %find(pert0~=0);
            for i2 = 1:length(figsufs)
                for i3 = 1:length(basefigs)
                    basefig = basefigs{i3};
                    fnamepre = 'fig31a';
                    fname = [fnamepre '-' num2str(i3) figsufs{i2}];
                    f = openfig(fullfile(netpath,'figs/',basefig));
                    colormap(parula);
                    cmap = colormap;cmap_size = size(cmap,1);
%                     set(gca,'Color','k');
                    f.Name = fname;
                    hax = f.Children;
                    nodes = hax.Children;
%                     nodes.MarkerFaceColor = 'w';nodes.MarkerEdgeColor = 'w';
                    nodes.ZData = myfuncs{i2}(pert0);
                    f.InvertHardcopy = 'off';
                    General.save_fig(f,fname,fullfile(obj.resultsPath,'figs'));
                    
                    %fig31b
                    markalpha = .2;
                    fnamepre = 'fig31b';
                    fname = [fnamepre '-' num2str(i3) figsufs{i2}];
                    f.Name = fname;
                    tau_As = General.load_var(fullfile(obj.resultsPath,'obj_properties',...
                        'single_node_pert_sol/','sol_x_2_4_1_22_10_tau_A.mat'));
                    tau_A = tau_As{1};
                    nodes.MarkerFaceAlpha = markalpha;nodes.MarkerEdgeAlpha = markalpha;
                    num_taus = 1;frames_per_tau = 256;
                    hold on;
                    X = nodes.XData(nodes_sub_inds);Y = nodes.YData(nodes_sub_inds);C = [1 1 1];%nodes.CData;
                    sz = nodes.SizeData;prev_ind = 1;same_ind_count = 1;C_ind = 1;
                    if length(sz) == 1
                        mysz = sz;
                    else
                        mysz = sz(nodes_sub_inds);
                    end
                    num_iter = num_taus*frames_per_tau;
                    for i1 = 1:num_iter
                        cur_t = i1*tau_A/frames_per_tau;
                        ind = find(sol_t<=cur_t,1,'last');
                        disp(['i1 = ' num2str(i1)]);
                        disp(['ind = ' num2str(ind)]);
                        
                        
                        if (ind ~= prev_ind) || i1 == num_iter
                            pert_prev_ind = pert(prev_ind,:);
                            
                            if ind==prev_ind && ind < size(pert,1)
%                                 same_ind_count = same_ind_count+1;
                                ind = prev_ind+1;
                                pert_ind = pert(ind,:);
                            elseif ind == size(pert,1)
                                pert_ind = pert_prev_ind - (pert_prev_ind - pert(prev_ind - 1,:));
                            else
                                pert_ind = pert(ind,:);
                            end
                            
                            
                            for i4 = 1:same_ind_count
                                C = cmap(C_ind,:);
                                z = pert_prev_ind - (i4*(pert_prev_ind-pert_ind)/same_ind_count);
                                s = scatter3(X,Y,myfuncs{i2}(z(nodes_sub_inds)),mysz,'MarkerFaceAlpha',markalpha,...
                                    'MarkerFaceColor',C,'MarkerEdgeColor',C,...
                                    'MarkerEdgeAlpha',markalpha);
                                disp(['C_ind = ' num2str(C_ind)]);
                                C_ind = C_ind + 1;
                                
                            end
                            same_ind_count = 1;
                            prev_ind = ind;
                        else
                            same_ind_count = same_ind_count + 1;
                        end
                        
                    end
                    cbh = colorbar('southoutside');cbh.Ticks = [0 .5 1];
                    cbh.TickLabels = {'$t=0$','$t=0.5\tau$','$t=1\tau$'};
                    cbh.TickLabelInterpreter = 'latex';
                    General.save_fig(f,fname,fullfile(obj.resultsPath,'figs'));
                end
            end
        end
        function plot_network_dts_3(obj)
            %fig32a
            fname = 'fig32a';
            subnetsize = 500;
            network = obj.networkName;
            netpath = fullfile('networks/',network);
            k = General.load_var(fullfile(netpath,'degree_vector'));
            [K,I] = sort(k,'descend');
            k_sub = k(I(1:subnetsize));
            sol_t = General.load_var(fullfile(obj.resultsPath,'obj_properties',...
                'single_node_pert_sol','sol_t_2_4_1_22_10.mat'));
            sol_x = General.load_var(fullfile(obj.resultsPath,'obj_properties',...
                'single_node_pert_sol','sol_x_2_4_1_22_10.mat'));
            pert = (sol_x' - obj.steady_state')';max_pert = max(abs(pert),[],'all');
            pert = pert/max_pert;
            pert0 = pert(1,:);
            f=figure('Name',fname,'NumberTitle','off');hax = axes;
%             f = figure(20);f.Name = fname;f.NumberTitle = 'off';
            hax = General.ThreeDim_network_layout_visualization(k_sub,log10(abs(pert0(I(1:subnetsize)))*9+1),subnetsize,fname,hax);
            General.save_fig(f,fname,fullfile(obj.resultsPath,'figs'));
            
            fnamepre = 'fig32b';
            fname = [fnamepre];% figsufs{i2}];
            f.Name = fname;
            tau_As = General.load_var(fullfile(obj.resultsPath,'obj_properties',...
                'single_node_pert_sol/','sol_x_2_4_1_22_10_tau_A.mat'));
            tau_A = tau_As{1};
%            nodes.MarkerFaceAlpha = .1;nodes.MarkerEdgeAlpha = .1;
            num_taus = 2;frames_per_tau = 3;
            hold(hax,'on');
            for i1 = 1:num_taus*frames_per_tau
                ind = find(sol_t>=i1*tau_A/frames_per_tau,1,'first');
                pert_ind = pert(ind,:);
                hax = General.ThreeDim_network_layout_visualization(k_sub,log10(abs(pert_ind(I(1:subnetsize)))*9+1),subnetsize,fname,hax);
            end
            General.save_fig(f,fname,fullfile(obj.resultsPath,'figs'));
        end
        %fig33*
        function obj = plot_Qs_distribution(obj)
            name = 'fig33a';desc = '$10^5$ Random binary perturbations $Q_{10}$ Distribution';
            f=figure('Name',name,'NumberTitle','off');
            Qs = General.load_var(fullfile(obj.resultsPath,'obj_properties/','Q10s_rand_binary_perts'));
            histogram(Qs);
            title({[name ' ' obj.scenarioName];obj.desc;desc},'interpreter','latex');
            xlabel('$Q_{10}$','Interpreter','latex');
%             ylabel('Interpreter','latex');
            legend;
            General.save_fig(f,name,fullfile(obj.resultsPath,'figs'));
        end
        function obj = plot_Qs_distribution_2(obj)
            name = 'fig33b';desc = '$Q_{10}$ Distribution';
            f=figure('Name',name,'NumberTitle','off');
            Qs1 = General.load_var(fullfile(obj.resultsPath,'obj_properties/','Q10s_single_node_perts.mat'));
            Qs2 = General.load_var(fullfile(obj.resultsPath,'obj_properties/','Q10s_double_node_perts.mat'));
            Qs3 = General.load_var(fullfile(obj.resultsPath,'obj_properties/','Q10s_triple_node_perts.mat'));
            Qs = [Qs1; Qs2; Qs3];
            h = histogram(Qs); h.Normalization = 'probability'; set(gca,'YScale','log');
            title({[name ' ' obj.scenarioName];obj.desc;desc},'interpreter','latex');
            xlabel('$Q_{10}$','Interpreter','latex');
%             ylabel('Interpreter','latex');
            legend;
            General.save_fig(f,name,fullfile(obj.resultsPath,'figs'));
        end
        
        %%
        function obj = set_networkNameSF1(obj)
            obj.networkName = 'SF1';
            obj.plot_localization2();
        end
        %%
        function decayTimes = find_decay_times(~,pert,sol_t)
            pert0 = pert(1,:);
            pert1 = exp(-1)*pert0;
            c = size(pert,2);
            decayTimes = zeros(1,c);
            for i=1:c
                perti = pert(:,i);
                pert1i = pert1(1,i);
                ind = find(perti < pert1i,1,'first');
                decayTimes(1,i) = sol_t(ind);
            end
        end
        function obj = save_fig(obj,f,name)
            if ~isfolder(fullfile(obj.resultsPath,'figs'))
                mkdir(fullfile(obj.resultsPath,'figs'));
            end
            try
                saveas(f,fullfile(obj.resultsPath,'figs',name),'fig');
            catch exception
                disp(exception.message);
            end
            try
                saveas(f,fullfile(obj.resultsPath,'figs',name),'png');
            catch exception
                disp(exception.message);
            end
            try
                saveas(f,fullfile(obj.resultsPath,'figs',name),'svg');
            catch exception
                disp(exception.message);
            end
        end
        function obj = save_obj(obj)
            save(fullfile(obj.resultsPath,[obj.scenarioName 'Obj.mat']),'obj','-v7.3');
        end
        function write_gephi_nodes_table(obj,sol_path,sol_x_filename,sol_t_filename,steady_state)
            x = General.load_var(fullfile(sol_path,sol_x_filename));
            p = (x' - steady_state')';
            t = General.load_var(fullfile(sol_path,sol_t_filename));
            [norms, norms_normed, thetas,H_A,tau_A,Q6] = EngineClass.calc_norms_thetas(p,t);
            General.save_var(norms,sol_path,[sol_x_filename '_norms']);
            General.save_var(thetas,sol_path,[sol_x_filename '_thetas']);
            General.save_var(norms_normed,sol_path,[sol_x_filename '_norms_normed']);
            %             [m1,i1] = max(thetas);
            p_normed = (p'./norms)';p_prop = (p'./sum(p'))';
            [m1,i1] = min(max(p_prop'));
            ts = [0 t(floor(i1/2)) t(i1)];
            for cur_t = ts
                filename = ['gephi_nodes_' sol_x_filename '_t_' num2str(floor(cur_t))];
                EngineClass.create_gephi_nodes_table(p_prop,t,cur_t,sol_path,filename);
                
                name = ['fig27a-t-' num2str(floor(cur_t))];
                figdesc = ['$\delta x(t)$ at $t = ' num2str(floor(cur_t)) '$'];
                f = figure('Name',name,'NumberTitle','off');
                ind = find(t>=cur_t,1);norm = norms(ind); theta = thetas(ind);
                [x, y] = pol2cart(theta,norm);
                c = compass(x,y);
                c1 = c(1);
                c1.LineWidth = 2;
                c1.Color = 'r';
                title({[name ' ' obj.scenarioName];obj.desc;figdesc},...
                    'Interpreter','latex');
                obj.save_fig(f,name);
                
            end
        end
        
    end
    methods (Static)
        function [] = save_var(var,path,folder_name,filename)
            varname = inputname(1);
            if ~isfolder(fullfile(path,folder_name))
                mkdir(fullfile(path,folder_name));
            end
            save(fullfile(path,folder_name,filename),'var','-v7.3');
        end
        function xout = test_fun(a,b)
            xout = a+b;
        end
        function [norms, norms_normed, thetas,H_A,tau_A,Q,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10] = calc_norms_thetas(x,ss,t,sys_half_life_amp,mu,k,xi)
            % Calc pert:
            x_sol = x;
            x = (x' - ss')';
            x_pert = x;
            % Calc pert/ss, or dlogx (dx/x)
            x = (x'./ss')';
            norms = vecnorm(x');
            norms_normed = norms/norms(1);
%             disp(['norms_normed(1) = ' num2str(norms_normed(1))]);
            ind = find(norms_normed <= norms_normed(1)/2,1,'first');
            H_A = t(ind);
            
            k_mu = k.^(-mu);
            abs_x_pert_init = abs(x_pert(1,:));
            x_pert_init = x_pert(1,:);
            abs_x_init = abs(x(1,:));
            x_init = x(1,:);
            %             Q = dot(abs_x_pert_init/norm(abs_x_pert_init),k_mu/norm(k_mu));
            Q = dot(x_pert_init/norm(x_pert_init),k_mu/norm(k_mu));
            %             Q = dot(abs_x_init/norm(abs_x_init),k_mu/norm(k_mu));
            Q2 = dot(x_pert_init/norm(x_pert_init,1),k_mu/norm(k_mu,1));
            Q3 = dot(abs(x_pert_init)/norm(x_pert_init,1),abs(k_mu)/norm(k_mu,1));
            Q4 = dot(abs(x_pert_init)/norm(x_pert_init,1),k_mu);
            Q5 = dot(abs(x_pert_init)/norm(x_pert_init,1),k_mu/norm(k_mu));
            Q6 = dot(abs(x_pert_init)/norm(x_pert_init,1),abs(k_mu)/max(k_mu));
            k_alpha = k.^(-xi - mu);
            Q7 = dot(abs(x_init),k_alpha)/(norm(x_init,1)*max(k_alpha));
            Q8 = dot(abs(x_init)/norm(x_init,1),k_alpha);
            Q9 = dot(abs(x_init)/norm(x_init,1),k_mu.*ss');
            k_xi = k.^(-xi);
            Q10 = dot(abs(x_init),k_alpha)/(dot(abs(x_init),k_xi)); %set4f
            x_normed = x'./norms;
            x_normed_init = x_normed(:,1);
            thetas = acos(dot(repmat(x_normed_init,1,size(x_normed,2)),x_normed));
            for i1 = 1:length(sys_half_life_amp)
                tau_A{i1} = H_A/sys_half_life_amp{i1};
            end
        end
        function write_norms_thetas_single_sol(path,x_filename,t_filename,ss,sys_half_life_amp,mu,k,sufstr,xi)
            sol_x = General.load_var(fullfile(path,x_filename));
            sol_t = General.load_var(fullfile(path,t_filename));
            [norms,norms_normed,thetas,H_A,tau_A,Q,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10] = EngineClass.calc_norms_thetas(sol_x,ss,sol_t,sys_half_life_amp,mu,k,xi);
            General.save_var(norms,path,[x_filename '_norms']);
            General.save_var(thetas,path,[x_filename '_thetas']);
            General.save_var(norms_normed,path,[x_filename '_norms_normed']);
            General.save_var(H_A,path,[x_filename '_H_A']);
            General.save_var(tau_A,path,[x_filename '_tau_A' sufstr]);
            General.save_var(Q,path,[x_filename '_Q']);
            General.save_var(Q2,path,[x_filename '_Q2']);
            General.save_var(Q3,path,[x_filename '_Q3']);
            General.save_var(Q4,path,[x_filename '_Q4']);
            General.save_var(Q5,path,[x_filename '_Q5']);
            General.save_var(Q6,path,[x_filename '_Q6']);
            General.save_var(Q7,path,[x_filename '_Q7']);
            General.save_var(Q8,path,[x_filename '_Q8']);
            General.save_var(Q9,path,[x_filename '_Q9']);
            General.save_var(Q10,path,[x_filename '_Q10']);
        end
        
        function write_norms_thetas_multi_sol(path,ss,sys_half_life_amp,mu,k,sufstr,xi)
            contents = dir(path);
            for j = 1:length(contents)
                content = contents(j);
                fname = content.name(1:end-4);
                if length(fname) > 2 && contains(fname,'_x_') && ~contains(fname,'_norms') && ~contains(fname,'_thetas') && ...
                        ~contains(fname,'_H_A') && ~contains(fname,'_Q') && ~contains(fname,'_tau_A') && ~contains(fname,'gephi')
                    tfname = strrep(fname,'x','t');
                    EngineClass.write_norms_thetas_single_sol(path,fname,tfname,ss,sys_half_life_amp,mu,k,sufstr,xi);
                end
            end
        end
        function str = array2str(arr)
            str = '';
            for n = arr
                str = [str num2str(n) '_'];
            end
            str = str(1:end-1);
        end
        function str = filename2legend(str_inp)
            str = strrep(str_inp,'_',',');
        end
        function create_gephi_nodes_table(solution_x,solution_t,time,filepath,filename)
            nodes = (1:size(solution_x,2))';
            ind = find(solution_t>=time,1);
            values = (solution_x(ind,:))';
            General.general_gephi_nodes_table(nodes,values,filepath,filename);
            %             T = array2table([nodes,nodes,values],'VariableNames',{'ID','Label','value'});
            %             writetable(T,fullfile(filepath,[filename '.csv']));
        end
        function [ind_bins_var,bins_edges,ind_bins_var_sizes] = set_bins_generic(numbins,values_vec,tol,cond_vec)
            ind_bins_var = cell(numbins,1);
            bins_edges = [];
            ind_bins_var_sizes = [];
            %             ind_bins_var = {};
            values_vec = values_vec(cond_vec);
            max_val = max(values_vec);
            min_val = min(values_vec);
            zero_inds = [];
            if min_val==0
                zero_inds = find(values_vec == 0);
                values_vec_no_zeros = values_vec;
                values_vec_no_zeros(zero_inds) = [];
                min_val = min(values_vec_no_zeros);
            end
            maxmin_val = max_val/min_val;
            c = maxmin_val ^ (1/numbins);
            n0 = log(min_val)/log(c);
            nn = log(max_val)/log(c);
            ns = n0:n0+numbins;
            s=0;
            for k = 2:numbins+1
                if k==2
                    ind = find(values_vec>=c^ns(k-1)-tol & values_vec<=c^ns(k));
                    ind = [zero_inds; ind];ind = sort(ind);
                    bins_edges(1) = c^ns(k-1);
                else
                    ind = find(values_vec>c^ns(k-1) & values_vec<=c^ns(k)+tol);
                    bins_edges(end+1) = c^ns(k-1);
                end
                %                 if ~isempty(ind)
                ind_bins_var{k-1,1} = ind;
                ind_bins_var_sizes(end+1) = size(ind,1);
                s = s + size(ind,1);
                %                 end
            end
            bins_edges(end+1) = c^ns(k);
            %Check bins
            num_vals = size(values_vec,1);
            if s == num_vals
                disp('Bins check passed')
            else
                disp('Bins check failed')
                disp(['s=' num2str(s)])
                disp(['num_vals = ' num2str(num_vals)])
            end
        end
        function [pert_norm, pert_pert0_dot_normed, half_life, sol_tau,ind_1,ind_2] = compute_pert_props(sol_t,sol_x,ss)
            pert = (sol_x' - ss')';
            pert0 = pert(1,:);
            pert_norm = vecnorm(pert');
            pert_normed = (pert'./pert_norm)';
            pert0_normed = pert_normed(1,:);
            pert0_normed_mat = pert0_normed' + zeros(size(pert'));
            pert_pert0_dot_normed = dot(pert_normed',pert0_normed_mat);
            ind1 = find(diff(pert_norm)<=0);half_life = [];ind_1 = ind1(1);
            for i = ind1
                ind2 = find(pert_norm < .5*pert_norm(i),1);
                if ~isempty(ind2)
                    half_life(end+1) = sol_t(ind2) - sol_t(i);
                    ind_2prev=i;
                else
                    ind_2 = ind_2prev;
                    break
                end
            end
            sol_tau = sol_t/half_life(1);
        end
        function concentration = compute_concentration(pert_t,I)
            concentration = sum(abs(pert_t(I)))/sum(abs(pert_t));
        end
        function conservation = compute_conservation(pert_t,pert_0)
            conservation = abs(dot(pert_t/norm(pert_t),pert_0/norm(pert_0)));
        end
        function [wmeans, wvars, columns_prop_ordered] = compute_wmeans_wvars(columns, weights)
            % columns is a matrix, we find the weighted mean and variance
            % of each column, according to the weight vector weights. Don't
            % forget to take absolute value of columns before calling
            % function if want mass distribution.
            [weights_ordered, weights_ordered_ind] = sort(weights);
            weights_mat = weights_ordered + zeros(size(columns));
            s = sum(columns);
            columns_prop_ordered = columns(weights_ordered_ind,:)./s;
            d = columns_prop_ordered.*weights_mat;clear weights_mat;
            wmeans = sum(d);clear d;
            wvars = sqrt(sum(columns_prop_ordered.*(weights_ordered-wmeans).^2));
        end
        function J = compute_J(Dii,Wij,C_D,C_W)
            J = C_W*(Wij - diag(diag(Wij))) + C_D*(diag(Dii));
        end
        function deg_vec = compute_degree_vector_weighted(matrix)
            deg_vec = sum(matrix,2);
        end
        function k_inn = compute_k_inn(matrix)
            ki = sum(matrix,2);
            tmp = matrix * ki;
            k_inn = tmp ./ ki;
        end
        function weighted_dot = compute_weighted_dot_product(v1,v2,weights_vector,weights_threshold,is_greater_than_threshold)
            if is_greater_than_threshold
                mask = weights_vector >= weights_threshold;
            else
                mask = weights_vector <= weights_threshold;
            end
            v1_masked = v1(mask);v1_masked_normed = v1_masked/norm(v1_masked);
            v2_masked = v2(mask);v2_masked_normed = v2_masked/norm(v2_masked);
            weighted_dot = dot(v1_masked_normed,v2_masked_normed);
        end
        function [weighted_dots,weighted_dots_means] = compute_weighted_dot_products(v1_set,v2_set,weights_vector,weights_thresholds,is_greater_than_threshold,vector_indices_for_averaging,is_first_n_indices)
            n = size(v1_set,2);
            m = size(weights_thresholds,1);
            o = size(vector_indices_for_averaging,1)-1;
            weighted_dots = zeros(n,m);
            weighted_dots_means = zeros(o,m);
            for j = 1:m
                weights_threshold = weights_thresholds(j,1);
                for i=1:n
                    v1 = v1_set(:,i);v2 = v2_set(:,i);
                    weighted_dots(i,j) = EngineClass.compute_weighted_dot_product(v1,v2,weights_vector,weights_threshold,is_greater_than_threshold);
                end
            end
            for k = 1:o
                vec_ind = vector_indices_for_averaging(k+1);
                if is_first_n_indices
                    weighted_dots_means(k,:) = mean(abs(weighted_dots(1:vec_ind,:)),1);
                else
                    weighted_dots_means(k,:) = mean(abs(weighted_dots(end-vec_ind+1:end,:)),1);
                end
            end
        end
        function [ind_bins_var,binned_vals] = set_bins_percentiles(numbins,values_vec)
            n = size(values_vec,1);
            bin_size = floor(n/numbins);
            [B,I] = sort(values_vec);
            %             sz = [bin_size,numbins];
            sz = [bin_size,n/bin_size];
            ind_bins_var = reshape(I,sz);
            binned_vals = reshape(B,sz);
        end
        function [ind_bins_var,bins_edges,ind_bins_var_sizes] = set_bins_nonlog(num_bins,values_vec)
            ind_bins_var = cell(num_bins,1);bins_edges = [];ind_bins_var_sizes = [];
            minval = min(values_vec); maxval = max(values_vec);
            range = maxval - minval;  step = range/num_bins;
            bins_edges(end+1) = minval;
            for k = 1:num_bins
                if k == 1
                    ind = find(values_vec >= minval + step*(k-1) & values_vec <= minval + step*k);
                else
                    ind = find(values_vec > minval + step*(k-1) & values_vec <= minval + step*k);
                end
                ind_bins_var{k,1} = ind;
                bins_edges(end+1) = minval + step*k;ind_bins_var_sizes(end+1,1) = length(ind);
            end
        end
        function [ordered_eigs] = reorder_eigvecs_nodes(eigvec_set, new_ordering)
            ordered_eigs = eigvec_set(new_ordering,:);
        end
        function M = create_random_matrix(size_in,isSymmetric)
            M = rand(size_in);
            if isSymmetric
                M = M - tril(M,-1) + tril(M',-1);
                if isequal(M,M')
                    disp('M is symmetric');
                else
                    disp('M is not symmetric');
                end
            end
        end
        function [] = check_eigen0()
            path = fullfile('tests','eigvec_test_1');
            sizes = [10, 100, 1000, 10000];
            EngineClass.save_var(sizes,path,'','sizes');
            for i1 = 1:size(sizes,2)
                s = sizes(i1);
                M = EngineClass.create_random_matrix(s,true);
                EngineClass.save_var(M,path,'matrices',['M-' num2str(s)]);
            end
        end
        function [] = check_eigen1()
            path = fullfile('tests','eigvec_test_1');
            mydata = load(fullfile(path,'sizes'),'var');sizes = mydata.var;
            n_max = 200;
            EngineClass.save_var(n_max,path,'','n_max');
            for i1 = 1:size(sizes,2)
                s = sizes(i1);
                disp(['s = ' num2str(s)]);
                mydata = load(fullfile(path,'matrices',['M-' num2str(s)]),'var');M=mydata.var;
                n = min(n_max,s);
                [v,d,flag] = eigs(M,n,'largestreal');
                if flag == 0
                    disp('Eigenvectors converged');
                else
                    disp('Eigenvectors did not converge');
                end
                EngineClass.save_var(v,path,'eigvecs',['M-' num2str(s) '-1-v']);
                EngineClass.save_var(d,path,'eigvals',['M-' num2str(s) '-1-d']);
            end
        end
        function [] = check_eigen2()
            path = fullfile('tests','eigvec_test_1');
            mydata = load(fullfile(path,'sizes'),'var');sizes = mydata.var;
            mydata = load(fullfile(path,'n_max'),'var');n_max = mydata.var;
            C_Ds = [2 5 10 100 1000];
            EngineClass.save_var(C_Ds,path,'','C_Ds');
            for i1 = 1:size(sizes,2)
                s = sizes(i1);
                mydata = load(fullfile(path,'matrices',['M-' num2str(s)]),'var');M=mydata.var;
                for i2 = 1:size(C_Ds,2)
                    C_D = C_Ds(i2);
                    disp(['C_D = ' num2str(C_D)]);
                    M2 = M - diag(diag(M)) + C_D*diag(diag(M));
                    n = min(n_max,s);
                    [v,d,flag] = eigs(M2,n,'largestreal');
                    if flag == 0
                        disp('Eigenvectors converged');
                    else
                        disp('Eigenvectors did not converge');
                    end
                    EngineClass.save_var(v,path,'eigvecs',['M-' num2str(s) '-' num2str(C_D) '-v']);
                    EngineClass.save_var(d,path,'eigvals',['M-' num2str(s) '-' num2str(C_D) '-d']);
                end
            end
        end
        
        function [] = check_eigen3()
            path = fullfile('tests','eigvec_test_1');
            mydata = load(fullfile(path,'sizes'),'var');sizes = mydata.var;
            mydata = load(fullfile(path,'C_Ds'),'var');C_Ds = mydata.var;
            mydata = load(fullfile(path,'n_max'),'var');n_max = mydata.var;
            for i1 = 1:size(sizes,2)
                disp(['i1 = ' num2str(i1)]);
                s = sizes(i1);
                n = min(s,n_max);
                mydata = load(fullfile(path,'eigvecs',['M-' num2str(s) '-' num2str(1) '-v']),'var');
                vecs1 = mydata.var;
                for i2 = 1:size(C_Ds,2)
                    disp(['i2 = ' num2str(i2)]);
                    C_D = C_Ds(i2);
                    MDots = zeros(n);
                    mydata = load(fullfile(path,'eigvecs',['M-' num2str(s) '-' num2str(C_D) '-v']),'var');
                    vecsC_D = mydata.var;
                    for i3 = 1:n
                        %disp(['i3 = ' num2str(i3)]);
                        vi3 = vecs1(:,i3);
                        for i4 = 1:n
                            %disp(['i4 = ' num2str(i4)]);
                            ui4 = vecsC_D(:,i4);
                            MDots(i3,i4) = dot(vi3,ui4);
                        end
                    end
                    EngineClass.save_var(MDots,path,'MDots',['M-' num2str(s) '-' num2str(C_D) '-MDots']);
                end
            end
        end
        function [] = check_eigen4()
            path = fullfile('tests','eigvec_test_1');
            mydata = load(fullfile(path,'sizes'),'var');sizes = mydata.var;
            mydata = load(fullfile(path,'C_Ds'),'var');C_Ds = mydata.var;
            for i1 = 1:size(sizes,2)
                s = sizes(i1);
                for i2 = 1:size(C_Ds,2)
                    C_D = C_Ds(i2);
                    mydata = load(fullfile(path,'MDots',['M-' num2str(s) '-' num2str(C_D) '-MDots']),'var');
                    MDots = mydata.var;
                    name = ['Eigenvectors Dot Product Comparison Matrix $S = ' num2str(s) '$, $C_D = ' num2str(C_D) '$'];
                    namef = ['fig1-M-' num2str(s) '-C_D-' num2str(C_D)];
                    f = figure('Name',namef,'NumberTitle','off');
                    image(abs(MDots),'CDatamapping','scaled');
                    colorbar;
                    ylabel('C_D = 1 Eigenvectors');
                    xlabel(['C_D = ' num2str(C_D) ' Eigenvectors']);
                    title({name;'$\mid v_i \cdot v_j \mid $'},'interpreter','latex');
                    EngineClass.save_fig_static(f,namef,path);
                end
            end
        end
        function [] = check_eigen5()
            path = fullfile('tests','eigvec_test_1');
            mydata = load(fullfile(path,'sizes'),'var');sizes = mydata.var;
            mydata = load(fullfile(path,'C_Ds'),'var');C_Ds = [1 mydata.var];
            for i1 = 1:size(sizes,2)
                s = sizes(i1);
                mydata = load(fullfile(path,'matrices',['M-' num2str(s)]),'var');M = mydata.var;
                namef = ['fig3-M-' num2str(s)];
                f = figure('Name',namef,'NumberTitle','off');
                image(M,'CDatamapping','scaled');
                colorbar;
                title(namef);
                EngineClass.save_fig_static(f,namef,path);
                namef = ['fig2-M-' num2str(s)];
                name = ['Eigenvalues Comparison $S = ' num2str(s) '$'];
                f = figure('Name',namef,'NumberTitle','off');
                legendStr = {};
                for i2 = 1:size(C_Ds,2)
                    C_D = C_Ds(i2);
                    mydata = load(fullfile(path,'eigvals',['M-' num2str(s) '-' num2str(C_D) '-d']),'var');
                    eigvals = diag(mydata.var);
                    plot(eigvals,'o');
                    hold on;
                    legendStr{end+1} = ['C_D = ' num2str(C_D)];
                end
                xlabel('index');
                ylabel('Eigenvalues');
                legend(legendStr);
                title(name,'Interpreter','latex');
                EngineClass.save_fig_static(f,namef,path);
            end
        end
        function [] = plot_image_static(M,name,namef,path,mytitle,xlabelstr,ylabelstr)
            f = figure('Name',name,'NumberTitle','off');
            image(abs(M),'CDatamapping','scaled');
            colorbar;set(gca,'ColorScale','log');
            title(mytitle);
            xlabel(xlabelstr,'Interpreter','latex');
            ylabel(ylabelstr);
            EngineClass.save_fig_static(f,namef,path);
        end
        function M2 = compute_M2(M,C_W,C_D)
            M2 = C_W*(M - diag(diag(M))) + C_D*(diag(diag(M)));
        end
        function [vecs,vals] = compute_eigs(M,n)
            [vecs,vals] = eigs(M,n,'largestreal');
        end
        function [e_v, e_l, e_sum] = compute_errors(vecs0,vals0,vecs1,vals1)
            v0 = vecs0(:,1);
            v1 = vecs1(:,1);
            vdot = dot(v0,v1);
            e_v = 1 - abs(vdot);
            l0 = vals0(1,1);
            l1 = vals1(1,1);
            e_l = l1 - l0;
            e_sum = abs(e_v) + abs(e_l);
        end
        function [e_vs,e_ls,e_sums] = compute_loop(A,M,C_D,C_W,n)
            %n=3;
            %A = magic(n);A = A - tril(A,-1) + tril(A',-1);
            %M = EngineClass.compute_M2(A,2,4);
            [vecs0,vals0] = EngineClass.compute_eigs(A,n);
            %[vecs1,vals1] = EngineClass.compute_eigs(M,n);
            %C_D = [-1:.5:2];C_W = [-1:.25:1.25];
            tic;
            e_vs = [];e_ls = []; e_sums = [];
            for i=1:size(C_D,2)
                disp(['i=' num2str(i)]);
                C_Di = C_D(i);
                parfor j=1:size(C_W,2)
                    disp(['j=' num2str(j)]);
                    M2 = EngineClass.compute_M2(M,C_Di,C_W(j));
                    [vecs,vals] = EngineClass.compute_eigs(M2,n);
                    [e_vs(j,i), e_ls(j,i), e_sums(j,i)] = EngineClass.compute_errors(vecs0,vals0,vecs,vals);
                end
            end
            toc;
            [X,Y] = meshgrid(C_D,C_W);
            path = fullfile('tests','eigvec_test_1');
            f=figure;surf(X,Y,e_vs);xlabel('C_D');ylabel('C_W');zlabel('e_v');
            title('Error in first eigenvector $1 - \mid v \cdot \hat{v} \mid$','Interpreter','latex');
            EngineClass.save_fig_static(f,'e_vs',path);
            f=figure;surf(X,Y,abs(e_ls));xlabel('C_D');ylabel('C_W');zlabel('|e_l|');
            title('Error in first eigenvalue $ \mid \lambda - \hat{\lambda} \mid $','Interpreter','latex');
            EngineClass.save_fig_static(f,'abs_e_ls',path);
            f=figure;surf(X,Y,e_sums);xlabel('C_D');ylabel('C_W');zlabel('|e_v| + |e_l|');
            title('Sum of errors in first eigenvector and eigenvalue');
            EngineClass.save_fig_static(f,'e_sum',path);
            figure;plot(C_D,e_vs(6,:));xlabel('C_D, C_W = .25');ylabel('e_vs');
            figure;plot(C_W,e_vs(:,4));xlabel('C_W, C_D = .5');ylabel('e_vs');
            figure;plot(C_D,abs(e_ls(6,:)));xlabel('C_D, C_W = .25');ylabel('e_ls');
            figure;plot(C_W,abs(e_ls(:,4)));xlabel('C_W, C_D = .5');ylabel('e_ls');
        end
        function [] = save_fig_static(f,name,path)
            if ~isfolder(fullfile(path,'figs'))
                mkdir(fullfile(path,'figs'));
            end
            try
                saveas(f,fullfile(path,'figs',name),'fig');
            catch exception
                disp(exception.message);
            end
            try
                saveas(f,fullfile(path,'figs',name),'png');
            catch exception
                disp(exception.message);
            end
        end
        function [M] = test_global()
            global M;
            M(end+1) = 1;
        end
    end
    
end

