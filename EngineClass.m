classdef EngineClass <  handle
    %ENGINECLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        scenarioName %String
        desc %String
        adjacencyMatrix %2D array
        initialValues %1D array same length as adjacency matrix
        maxTime %positive integer
        maxDerivative % positive double
        solverTimeStep % positive integer
        randSeed %positive integer
        f_M0 %function handle
        f_M1 %function handle
        f_M2 %function handle
        f_R %function handle
        f_dM0 %function handle, derivate of M0
        f_dM1 %derivative of M1
        f_dM2 %derivative of M2
        Dii % vector representing diagonal of Jacobian
        Wij % matrix representing Jacobian without diagonal
        Dii_asy
        Wij_asy
        Dii_num
        Wij_num
        M2_i_bigodot
        M2_k_bigodot
        perturbation
        half_life
        eigenvalues_ana
        eigenvectors_ana
        eigenvalues_asy
        eigenvectors_asy
        eigenvalues_num
        eigenvectors_num
        node_relax_time
        sys_relax_time
        node_half_life
        sys_half_life
        N
        k_nn
        
        difEqSolver %function handle
        absTol % positive double
        relTol % positive double
        solution_t = 0;
        solution_x ;
        resultsPath
        header
        degree_vector
        degree_vector_weighted
        steady_state
        steady_state_calculated
        mu = 0
        nu = 0
        rho = 0
        eta = 0
    end
    methods
        function obj = EngineClass(propertiesMap)
            %ENGINECLASS Construct an instance of this class
            %   Loop over keys of propertiesMap, assign each value to
            %   property of EngineClass
            for key = keys(propertiesMap)
                eval(['obj.' key{1} '=propertiesMap(key{1})']);
            end
            obj.solution_x = obj.initialValues';
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
            obj.calculate_degree();
            obj.calculate_degree_weighted();
            obj.set_dM0();
            obj.set_dM1();
            obj.set_dM2();
            obj.set_R();
            obj.set_const_functions();
            obj.set_N();
            obj.set_knn();
            obj.set_Dii_asy();
            obj.set_Wij_asy();
            obj.save_obj();
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
        function obj = solve(obj)
            %Solve the system
            tic;
            opts = odeset('RelTol',obj.relTol,'AbsTol',obj.absTol);
            odefun = @(t,x) (obj.f_M0(x) + (obj.adjacencyMatrix*obj.f_M2(x)).*obj.f_M1(x));
            t = 0; %initial step start time
            init = obj.initialValues; % set initial value for step
            setBreak = false; % stop while loop?
            while ~setBreak
                % check if this step passes maxTime and set stepEndTime
                if t + obj.solverTimeStep >= obj.maxTime
                    stepEndTime = obj.maxTime;
                    setBreak = true;
                else
                    stepEndTime = t + obj.solverTimeStep;
                end
                % run solver step
                display(['stepEndTime = ' num2str(stepEndTime)]);
                [sol_t,sol_x] = obj.difEqSolver(odefun,[t stepEndTime],init,opts);
                t = stepEndTime;
                init = sol_x(end,:);
                % append results to solution_t and solution_x
                obj.solution_t = [obj.solution_t;sol_t(2:end,:)];
                obj.solution_x = [obj.solution_x;sol_x(2:end,:)];
            end
            toc;
            obj.set_steady_state();
            obj.set_M2_i_bigodot();
            obj.set_steady_state_calculated();
            obj.set_Dii();
            obj.set_Wij();
            obj.save_obj();
        end
        function obj = set_steady_state(obj)
            obj.steady_state = obj.solution_x(end,:);
            disp('set_steady_state(obj): obj.steady_state = ');
            disp(obj.steady_state);
        end
        function obj = set_M2_i_bigodot(obj)
            obj.M2_i_bigodot = zeros(size(obj.initialValues));
            M2_x = obj.f_M2(obj.steady_state);
            neighbor_sum_M2_x = obj.adjacencyMatrix * M2_x';
            obj.M2_i_bigodot = neighbor_sum_M2_x ./ obj.degree_vector_weighted;
            disp('set_M2_i_bigdot(obj): obj.M2_i_bigdot = ');
            disp(obj.M2_i_bigodot);
        end
        function obj = set_steady_state_calculated(obj)
            syms x;
            R_inv = finverse(obj.f_R(x));
            x = 1./(obj.degree_vector_weighted .* obj.M2_i_bigodot);
            obj.steady_state_calculated = double(subs(R_inv));
            disp('set_steady_state_calculated(obj): obj.steady_state_calculated = ');
            disp(obj.steady_state_calculated);
        end
        function obj = set_Dii(obj)
            x = obj.steady_state';
            obj.Dii = double(obj.f_dM0(x)) + ...
                double(obj.f_dM1(x)).*double(obj.adjacencyMatrix*obj.f_M2(x));
            disp('set_Dii(obj): obj.Dii = ');
            %disp(obj.Dii);
        end
        function obj = set_Wij(obj)
            x = obj.steady_state;
            obj.Wij = double(obj.f_M1(x)).*...
                obj.adjacencyMatrix.*double(obj.f_dM2(x));
            disp('set_Wij(obj): obj.Wij = ');
            %disp(obj.Wij);
        end
        function obj = set_N(obj)
            obj.N = size(obj.adjacencyMatrix,1);
        end
        function obj = set_knn(obj)
            tmp = obj.adjacencyMatrix*obj.degree_vector_weighted;
            tmp = tmp./obj.degree_vector_weighted;
            obj.k_nn = sum(tmp)/obj.N;
        end
        function obj = set_Dii_asy(obj)
            obj.Dii_asy = -obj.k_nn^obj.eta*obj.degree_vector_weighted.^obj.mu;
        end
        function obj = set_Wij_asy(obj)
            obj.Wij_asy = obj.adjacencyMatrix .* (obj.degree_vector_weighted.^obj.nu * obj.degree_vector_weighted'.^obj.rho);
        end
        function obj = print_output(obj) %create output file
            T = array2table([obj.solution_t,obj.solution_x],'VariableNames',obj.header);
            if ~isfolder(obj.resultsPath)
                mkdir(obj.resultsPath)
            end
            %writetable(T,fullfile(obj.resultsPath,'output.csv'));
            save(fullfile(obj.resultsPath,'output.mat'),'T');
        end
        function obj = plot_results(obj, isDilute)
            if isDilute
                step = 100;
                suffix = '-diluted';
            else
                step = 1;
                suffix = [];
            end
%             ind_t = 1:step:length(obj.solution_t);
             ind_x = 1:step:size(obj.solution_x,2);
             ind_t = 1:step:length(obj.solution_t);
            t = obj.solution_t(ind_t);
            x = obj.solution_x(ind_t,ind_x);
            name = ['fig1' suffix];
            f = figure('Name',name,'NumberTitle','off');
            plot(t,x,'.-','MarkerSize',12);
            legend(obj.header(2:step:end));
            title({[name ' ' obj.scenarioName];obj.desc});
            xlabel(obj.header{1});
            ylabel('x');
            obj.save_fig(f,name);
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
        function obj = plot_steady_vs_degree(obj)
            name = 'fig3';
            f = figure('Name',name,'NumberTitle','off');
            plot(obj.degree_vector,obj.steady_state,'.','MarkerSize',12);
            title({[name ' ' obj.scenarioName];obj.desc});
            xlabel('k - node degree');
            ylabel('Steady State');
            obj.save_fig(f,name);
            name = 'fig3a';
            f = figure('Name',name,'NumberTitle','off');
            plot(log(obj.degree_vector),log(obj.steady_state),'.','MarkerSize',12);
            title({[name ' ' obj.scenarioName];obj.desc});
            xlabel('log(k) - node degree');
            ylabel('log(Steady State)');
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
            xlabel('log(k) - weighted node degree');
            ylabel('log(1 - Steady State)');
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
            name = 'fig4a';
            f = figure('Name',name,'NumberTitle','off');
            x = obj.degree_vector_weighted;
            y = obj.Dii;
            loglog(x,y,'.','MarkerSize',12);
            hold on;
            y1 = obj.Dii_asy;
            loglog(x,y1,'^','MarkerSize',10);
            xlabel('k_i - weighted');
            ylabel('J_{ii}');
            title({[name ' ' obj.scenarioName];obj.desc});
            legend('Analytic Jacobian', 'Asymptotic Jacobian');
            obj.save_fig(f,name);
            name = 'fig4b';
            f = figure('Name',name,'NumberTitle','off');
            x1 = x * x';
            y = obj.Wij;
            loglog(x1(:),y(:),'.','MarkerSize',12);
            hold on;
            y1 = obj.Wij_asy;
            loglog(x1(:),y1(:),'^','MarkerSize',10);
            xlabel('k_ik_j - weighted');
            ylabel('W_{ij}');
            title({[name ' ' obj.scenarioName];obj.desc});
            legend('Analytic Jacobian','Asymptotic Jacobian');
            obj.save_fig(f,name);
            name = 'fig4c';
            f = figure('Name',name,'NumberTitle','off');
            x2 = (x.^obj.nu) * (x.^obj.rho)';
            loglog(x2(:),y1(:),'.','MarkerSize',12);
            hold on;
            loglog(x2(:),y(:),'^','MarkerSize',10);
            xlabel(['k_i^{\nu}k_j^{\rho} - weighted, \nu = ' num2str(obj.nu)...
                ', \rho = ' num2str(obj.rho)]);
            ylabel('W_{ij}');
            title({[name ' ' obj.scenarioName];obj.desc});
            legend('Analytic Jacobian','Asymptotic Jacobian');
            obj.save_fig(f,name);
        end
        function obj = save_fig(obj,f,name)
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
        end
        function obj = save_obj(obj)
            save(fullfile(obj.resultsPath,[obj.scenarioName 'Obj.mat']),'obj');
        end
        % plot results
%             f = figure;
%             plot(obj.solution_t,obj.solution_x,'.-','MarkerSize',12);
%             legend(obj.header(2:end));
%             title(obj.scenarioName);
%             xlabel(obj.header{1});
%             ylabel('x');
%             text(1,.5,num2str(obj.adjacencyMatrix));
%             if ~isfolder(fullfile(obj.resultsPath,'figs'))
%                 mkdir(fullfile(obj.resultsPath,'figs'));
%             end
%             try
%                 saveas(f,fullfile(obj.resultsPath,'figs','fig1.fig'),'fig');
%             catch exception
%                 display(exception.message);
%             end
    end
    
end

