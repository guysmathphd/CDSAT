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
        Dii_ana % vector representing diagonal of Jacobian
        Wij_ana % matrix representing Jacobian without diagonal
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
        eigenvalues_asy_permuted
        eigenvectors_asy
        eigenvectors_asy_permuted
        eigenvalues_num
        eigenvectors_num
        node_relax_time
        sys_relax_time
        node_half_life
        sys_half_life
        N
        k_nn
        ki_nn
        bins
        binsij
        kbinned
        kijbinned
        Dii_anabinned
        Dii_asybinned
        Wij_anabinned
        Wij_asybinned
        numeigen
        eigvec_dist_comparison_mat_ana2asy
        eigvec_angle_comparison_mat_ana2asy
        eigvec_dist_comparison_mat_ana2asy_permuted
        eigvec_angle_comparison_mat_ana2asy_permuted
        permutation_eigvec_ana2asy
        
        
        difEqSolver %function handle
        absTol % positive double
        relTol % positive double
        solution_t = 0;
        solution_x ;
        solution_t_eigvecana;
        solution_t_eigvecasy;
        solution_x_eigvecana;
        solution_x_eigvecasy;
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
            obj.set_eig_asy();
            obj.set_bins();
            obj.set_kbinned();
            obj.set_Dii_asybinned();
            obj.set_Wij_asybinned();
            obj.save_obj();
        end
        function obj = set_bins(obj)
            obj.bins = cell(1,obj.numbins);
            obj.binsij = cell(1,obj.numbins);
            obj.kijbinned = zeros(obj.numbins,1);
            kmax = max(obj.degree_vector_weighted);
            kmin = min(obj.degree_vector_weighted);
            kmaxmin = kmax/kmin;
            c = kmaxmin ^ (1/obj.numbins);
            n0 = log(kmin)/log(c);
            nn = log(kmax)/log(c);
            ns = n0:n0+obj.numbins;
            %c = kmax^(1/obj.numbins);
            x = obj.degree_vector_weighted;
            x2 = (x.^obj.nu) * (x.^obj.rho)';
            x2max = max(max(x2));
            x2min = min(min(x2));
            x2maxmin = x2max/x2min;
            c1 = x2maxmin ^ (1/obj.numbins);
            m0 = log(x2min)/log(c);
            mm = log(x2max)/log(c);
            ms = m0:m0+obj.numbins;
            %c1 = x2max^(1/obj.numbins);
            l = 0;
            for k = 2:obj.numbins+1
                if k == 2
                    ind = find(obj.degree_vector_weighted>=c^ns(k-1) & obj.degree_vector_weighted<=c^ns(k));
                else
                    ind = find(obj.degree_vector_weighted>c^ns(k-1) & obj.degree_vector_weighted<=c^ns(k));
                end
                obj.bins{k-1} = ind;
                l = l + 1;
                disp(l);
                w = [];
                for i = 1:obj.N
                    %ki = obj.degree_vector_weighted(i);
                    for j = 1:obj.N
                        %kj = obj.degree_vector_weighted(j);
                        kk = x2(i,j);
                        %kk = ki^obj.nu * kj^obj.rho;
                        if obj.adjacencyMatrix(i,j) > 0 && (kk > c1^ms(k-1) || (k==2 && kk==c1^ms(k-1))) && kk <=c1^ms(k)
                            obj.binsij{k-1}(end+1,:) = [i j];
                            w(end+1) = kk;
                        end
                    end
                end
                if ~isempty(w)
                    obj.kijbinned(k-1,1) = mean(w);
                end
            end            
        end
        function obj = set_kbinned(obj)
            obj.kbinned = zeros(obj.numbins,1);
            for i=1:obj.numbins
                ind = obj.bins{i};
                if ~isempty(ind)
                    obj.kbinned(i,1) = mean(obj.degree_vector_weighted(ind));
                end
            end
        end
        function obj = set_Dii_anabinned(obj)
            obj.Dii_anabinned = zeros(obj.numbins,1);
            for i=1:obj.numbins
                ind = obj.bins{i};
                if ~isempty(ind)
                    obj.Dii_anabinned(i,1) = mean(obj.Dii_ana(ind));
                end
            end
        end
        function obj = set_Dii_asybinned(obj)
            obj.Dii_asybinned = zeros(obj.numbins,1);
            for i = 1:obj.numbins
                ind = obj.bins{i};
                if ~isempty(ind)
                    obj.Dii_asybinned(i,1) = mean(obj.Dii_asy(ind));
                end
            end
        end
        function obj = set_Wij_anabinned(obj)
            obj.Wij_anabinned = zeros(obj.numbins,1);
            for i = 1:obj.numbins
                ind = obj.binsij{i};
                if ~isempty(ind)
                    w = [];
                    for k = 1:size(ind,1)
                        r = ind(k,1);
                        c = ind(k,2);
                        w(end+1) = obj.Wij_ana(r,c);
                    end
                    obj.Wij_anabinned(i,1) = mean(w);
                end
            end
        end
        function obj = set_Wij_asybinned(obj)
            obj.Wij_asybinned = zeros(obj.numbins,1);
            for i = 1:obj.numbins
                ind = obj.binsij{i};
                if ~isempty(ind)
                    w = [];
                    for k = 1:size(ind,1)
                        r = ind(k,1);
                        c = ind(k,2);
                        w(end+1) = obj.Wij_asy(r,c);
                    end
                    obj.Wij_asybinned(i,1) = mean(w);
                end
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
        function obj = solve(obj,useEigVecAna,useEigVecAsy)
            %Solve the system
            tic;
            %nn=10;
            nn=1;
            eps = 1;
            ss = obj.steady_state;
            pert = 0;
            if useEigVecAsy
                pert = eps * obj.eigenvectors_asy_permuted(:,1:nn);
                inits = ss' + eps*pert;
                obj.solution_t_eigvecasy = cell(1,nn);
                obj.solution_t_eigvecasy = cell(1,nn);
                for i = 1:nn
                    obj.solution_t_eigvecasy{1,i} = 0;
                    obj.solution_x_eigvecasy{1,i} = inits(:,i)';
                end
            elseif useEigVecAna
                pert = eps * obj.eigenvectors_ana(:,1:nn);
                inits = ss' + eps*pert;
                obj.solution_t_eigvecana = cell(1,nn);
                obj.solution_t_eigvecana = cell(1,nn);
                for i = 1:nn
                    obj.solution_t_eigvecana{1,i} = 0;
                    obj.solution_x_eigvecana{1,i} = inits(:,i)';
                end
            else
                nn=1;
                pert = obj.initialValues;
                inits = ss' + pert;
            end
            
            opts = odeset('RelTol',obj.relTol,'AbsTol',obj.absTol);
            odefun = @(t,x) (obj.f_M0(x) + (obj.adjacencyMatrix*obj.f_M2(x)).*obj.f_M1(x));
            for ii = 1:nn
                disp(['ii = ' num2str(ii)]);
                t = 0; %initial step start time
                init = inits(:,ii);
                %init = obj.initialValues; % set initial value for step
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
                    if useEigVecAna
                        obj.solution_t_eigvecana{ii}(end+1:end+length(sol_t)-1,1) = sol_t(2:end,:);
                        obj.solution_x_eigvecana{ii}(end+1:end+size((sol_x),1)-1,:) = sol_x(2:end,:);
                    elseif useEigVecAsy
                        obj.solution_t_eigvecasy{ii}(end+1:end+length(sol_t)-1) = sol_t(2:end,:);
                        obj.solution_x_eigvecasy{ii}(end+1:end+size((sol_x),1)-1,:) = sol_x(2:end,:);
                    else
                        obj.solution_t = [obj.solution_t;sol_t(2:end,:)];
                        obj.solution_x = [obj.solution_x;sol_x(2:end,:)];
                    end
                end
            end
            toc;
            if ~useEigVecAna && ~useEigVecAsy
                obj.set_steady_state();
                obj.set_M2_i_bigodot();
                obj.set_steady_state_calculated();
                obj.set_Dii_ana();
                obj.set_Wij_ana();
                obj.set_eig_ana();
                obj.set_Dii_anabinned();
                obj.set_Wij_anabinned();
                obj.set_eigvec_comparison_mats(true,false);
                obj.set_permutation_eigvec_ana2asy();
                obj.set_eig_asy_permuted();
                obj.set_eigvec_comparison_mats(false,true);
            end
            obj.save_obj();
        end
        function obj = set_steady_state(obj)
            obj.steady_state = obj.solution_x(end,:);
            disp('set_steady_state(obj): obj.steady_state = ');
            %disp(obj.steady_state);
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
            obj.Wij_ana = obj.adjacencyMatrix .* (double(obj.f_M1(x))'*...
                double(obj.f_dM2(x)));                
            disp('set_Wij_ana(obj): obj.Wij_ana = ');
            %disp(obj.Wij_ana);
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
            obj.Wij_asy = obj.adjacencyMatrix .* (obj.degree_vector_weighted.^obj.nu * obj.degree_vector_weighted'.^obj.rho);
        end
        function obj = set_eig_ana(obj)
            J = obj.Wij_ana;
            J = J - diag(diag(J));
            J = J + diag(obj.Dii_ana);
            %[v,d] = eig(J);
            [v,d] = eigs(J,min(size(J,1),obj.numeigen),'largestreal');
            obj.eigenvalues_ana = diag(d);
            obj.eigenvectors_ana = v;
        end
        function obj = set_eig_asy(obj)
            J = obj.Wij_asy;
            J = J - diag(diag(J));
            J = J + diag(obj.Dii_asy);
            %[v,d] = eig(J);
            [v,d] = eigs(J,min(size(J,1),obj.numeigen),'largestreal');
            obj.eigenvalues_asy = diag(d);
            obj.eigenvectors_asy = v;
        end
        function obj = set_eigvec_comparison_mats(obj,computeRegular,computePermuted)
            e1 = obj.eigenvectors_ana;
            e2 = obj.eigenvectors_asy;
            e3 = obj.eigenvectors_asy_permuted;
            M = zeros(size(e1,2));
            NN = M;
            Mp = M;
            NNp = M;
            for i = 1:size(e1,2)
                if mod(i,100)==0
                    disp(i);
                end
                v1 = e1(:,i); % ith analytic eigenvector
                if computeRegular
                    M(i,:) = vecnorm(e2 - v1);
                end
                for j = 1:size(e2,2)
                    if computeRegular
                        v2 = e2(:,j);
                        NN(i,j) = acos(dot(v1,v2))*180/pi;
                    end
                    if computePermuted
                        v3 = e3(:,j);
                        NNp(i,j) = acos(dot(v1,v3))*180/pi;
                    end
                end
                if computePermuted
                    Mp(i,:) = vecnorm(e3-v1);
                end
            end
            obj.eigvec_dist_comparison_mat_ana2asy = M;
            obj.eigvec_angle_comparison_mat_ana2asy = NN;
            obj.eigvec_dist_comparison_mat_ana2asy_permuted = Mp;
            obj.eigvec_angle_comparison_mat_ana2asy_permuted = NNp;
            if ~isequal(real(NN),NN)
                disp('Warning: eigvec_angle_comparison_mat has complex values');
            end
        end
        function obj = set_permutation_eigvec_ana2asy(obj)
            NN = obj.eigvec_angle_comparison_mat_ana2asy;
            NN = abs(90 - NN);
            [M,I] = max(NN,[],2);
            M = M < 70;
            ii = 1:size(I,1);
            I(M) = ii(M);
            [val,ind] = unique(I);
            if isequal(ind,ii)
                disp('I is a permutation')
            end
            obj.permutation_eigvec_ana2asy = I';
        end
        function obj = set_eig_asy_permuted(obj)
            obj.eigenvalues_asy_permuted = obj.eigenvalues_asy(obj.permutation_eigvec_ana2asy);
            obj.eigenvectors_asy_permuted = obj.eigenvectors_asy(:,obj.permutation_eigvec_ana2asy);
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
        function obj = plot_results2(obj)
            n=1;
            angles_ana = cell(1,n);
            angles_asy = cell(1,n);
            step = 100;
            ss = obj.steady_state';
            ss = ss/norm(ss);
            for i = 1:n
                x1 = obj.solution_x_eigvecana{1,i};
                x2 = obj.solution_x_eigvecasy{1,i};
                n1 = size(x1,1);
                n2 = size(x2,1);
                for j = 1:n1
                    v1 = x1(j,:);
                    v1 = v1/norm(v1);
                    angles_ana{1,i}(j,1) = real(acos(dot(ss,v1)))*180/pi;
                end
                for j = 1:n2
                    v2 = x2(j,:);
                    v2 = v2/norm(v2);
                    angles_asy{1,i}(j,1) = real(acos(dot(ss,v2)))*180/pi;
                end
            end
            name = 'fig2a';
            f = figure('Name',name,'NumberTitle','off');
            
            CM = jet(n);
            hold on;
            legendStr = {};
            for i = 1:n
                plot(obj.solution_t_eigvecana{1,i}(1:step:end,1),angles_ana{1,i}(1:step:end,1),'.','Color',CM(i,:),'MarkerSize',12);
                legendStr{end+1} = ['x_0 = ss + \epsilon * v_{' num2str(i) ',ana'];
            end
            for i = 1:n
                plot(obj.solution_t_eigvecasy{1,i}(1:step:end,1),angles_asy{1,i}(1:step:end,1),'^','Color',CM(i,:));
                legendStr{end+1} = ['x_0 = ss + \epsilon * v_{' num2str(i) ',asy'];
            end
            figdesc = 'Angle from Steady State';
            title({[name ' ' obj.scenarioName];figdesc});
            xlabel('Time');
            ylabel('\theta_{x,ss}');
            legend(legendStr);
            obj.save_fig(f,name);
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
            y = obj.Dii_ana;
            loglog(x,y,'.','MarkerSize',12);
            hold on;
            y1 = obj.Dii_asy;
            loglog(x,y1,'^','MarkerSize',10);
            xlabel('k_i - weighted');
            ylabel('J_{ii}');
            title({[name ' ' obj.scenarioName];obj.desc});
            legend('Analytic Jacobian', 'Asymptotic Jacobian');
            obj.save_fig(f,name);
            %%%%
            name = 'fig4b';
            f = figure('Name',name,'NumberTitle','off');
            x1 = x * x';
            y = obj.Wij_ana;
            loglog(x1(:),y(:),'.','MarkerSize',12);
            hold on;
            y1 = obj.Wij_asy;
            loglog(x1(:),y1(:),'^','MarkerSize',10);
            xlabel('k_ik_j - weighted');
            ylabel('W_{ij}');
            title({[name ' ' obj.scenarioName];obj.desc});
            legend('Analytic Jacobian','Asymptotic Jacobian');
            obj.save_fig(f,name);
            %%%%
            name = 'fig4c';
            f = figure('Name',name,'NumberTitle','off');
            x2 = (x.^obj.nu) * (x.^obj.rho)';
            loglog(x2(:),y(:),'.','MarkerSize',12);
            hold on;
            loglog(x2(:),y1(:),'^','MarkerSize',10);
            xlabel(['k_i^{\nu}k_j^{\rho} - weighted, \nu = ' num2str(obj.nu)...
                ', \rho = ' num2str(obj.rho)]);
            ylabel('W_{ij}');
            title({[name ' ' obj.scenarioName];obj.desc});
            legend('Analytic Jacobian','Asymptotic Jacobian');
            obj.save_fig(f,name);
            %%%%
            name = 'fig4d';
            figdesc = 'Dii Real and Theoretical binned';
            f = figure('Name',name,'NumberTitle','off');
            x = obj.kbinned;
            y1 = obj.Dii_anabinned;
            y2 = obj.Dii_asybinned;
            loglog(x,y1,'s','MarkerSize',12);
            hold on;
            loglog(x,y2,'^','MarkerSize',10);
            xlabel('k_i - weighted');
            ylabel('J_{ii}');
            title({[name ' ' obj.scenarioName];obj.desc;figdesc});
            legend('Analytic Jacobian','Asymptotic Jacobian');
            obj.save_fig(f,name);
            %%%%
            name = 'fig4e';
            figdesc = 'Wij Real and Theoretical binned';
            f = figure('Name',name,'NumberTitle','off');
            x = obj.kijbinned;
            y1 = obj.Wij_anabinned;
            y2 = obj.Wij_asybinned;
            loglog(x,y1,'s','MarkerSize',12);
            hold on;
            loglog(x,y2,'^','MarkerSize',10);
            xlabel(['k_i^{\nu}k_j^{\rho} - weighted, \nu = ' num2str(obj.nu)...
                ', \rho = ' num2str(obj.rho)]);
            ylabel('W_{ij}');
            title({[name ' ' obj.scenarioName];obj.desc;figdesc});
            legend('Analytic Jacobian','Asymptotic Jacobian');
            obj.save_fig(f,name);
        end
        function plot_eigenvalues(obj)
            name = 'fig5a';
            figdesc = 'Jacobian Eigenvalues';
            f = figure('Name',name,'NumberTitle','off');
            plot(real(obj.eigenvalues_ana),'*-');
            hold on;
            plot(real(obj.eigenvalues_asy),'^-');
            xlabel('n');
            ylabel('real(\lambda_n)');
            title({[name ' ' obj.scenarioName];obj.desc;figdesc});
            legend('Analytic Jacobian','Asymptotic Jacobian');
            obj.save_fig(f,name);
        end
        function plot_eigenvectors(obj,usePermuted)
            if usePermuted
                e2 = obj.eigenvectors_asy_permuted;
                suf = '-Permuted';
                dist = obj.eigvec_dist_comparison_mat_ana2asy_permuted;
                angle = obj.eigvec_angle_comparison_mat_ana2asy_permuted;
            else
                e2 = obj.eigenvectors_asy;
                suf = '';
                dist = obj.eigvec_dist_comparison_mat_ana2asy;
                angle = obj.eigvec_angle_comparison_mat_ana2asy;
            end
            e1 = obj.eigenvectors_ana;                        
            %%% fig6a
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
            %%% fig6c
            name = ['fig6c' suf];
            figdesc = ['Jacobian Eigenvectors Distance Comparison Matrix' suf];
            f = figure('Name',name,'NumberTitle','off');
            image(dist,'CDatamapping','scaled');
            colorbar;
            ylabel('Analytic Eigenvectors');
            xlabel(['Asymptotic Eigenvectors' suf]);
            title({[name ' ' obj.scenarioName];obj.desc;figdesc;'$\mid v-\hat{v} \mid$'},'interpreter','latex');
            obj.save_fig(f,name);
            %%% fig6d
            name = ['fig6d' suf];
            figdesc = ['Jacobian Eigenvectors Angle Comparison Matrix' suf];
            f = figure('Name',name,'NumberTitle','off');
            image(angle,'CDatamapping','scaled');
            colorbar;
            ylabel('Analytic Eigenvectors');
            xlabel(['Asymptotic Eigenvectors' suf]);
            title({[name ' ' obj.scenarioName];obj.desc;figdesc;'$\theta_{v,\hat{v}} [deg]$'},'interpreter','latex');
            obj.save_fig(f,name);
        end
        function obj = plot_eigenvectors2(obj,isFirst10,isKinn,isPermuted,isLogLog)
            %%% fig7a,fig7b, fig7c,fig7d,fig7e,fig7f -permuted
            if isKinn == true
                fnumab = 'fig7b';
                fnumcd = 'fig7d';
                x = obj.ki_nn;
                figdescab = 'Eigenvector Elements vs Node Average Nearest Neighbor Degree';
                figdesccd = 'Eigenvector Elements Mean vs Node Average Nearest Neighbor Degree';
                xlab = 'k_{nn,i} - weighted';
            else
                fnumab = 'fig7a';
                fnumcd = 'fig7c';
                x = obj.degree_vector_weighted;
                figdescab = 'Eigenvector Elements vs Node Degree';
                figdesccd = 'Eigenvector Elements Mean vs Node Degree';
                xlab = 'k_i - weighted';
            end
            if isFirst10 == true
                suffix = '-first10';
                n = 10;
            else
                n = size(obj.eigenvectors_ana,2);
                suffix = '';
            end
            if isPermuted
                suffix2 = '-permuted';
                eigvecasy = obj.eigenvectors_asy_permuted;
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
            CM = jet(n);
            %% fig7a fig7b
            name = [fnumab suffix suffix2 suffix3];
            f = figure('Name',name,'NumberTitle','off');
            hold on;
            for i=1:n
                %myplot(x,abs4log(obj.eigenvectors_ana(:,i)),'.','MarkerSize',18,'Color',CM(i,:));
                plot(log10(x),log10(abs4log(obj.eigenvectors_ana(:,i))),'.','MarkerSize',18,'Color',CM(i,:));
            end
            for i=1:n
                %myplot(x,abs4log(eigvecasy(:,i)),'^','MarkerSize',12,'Color',CM(i,:));
                plot(log10(x),log10(abs4log(eigvecasy(:,i))),'^','MarkerSize',12,'Color',CM(i,:));
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
            name = [fnumcd suffix suffix2 suffix3];
            f = figure('Name',name,'NumberTitle','off');
            hold on;
            myplot(x,abs4log(mean(obj.eigenvectors_ana(:,1:n),2)),'.','MarkerSize',18);
            myplot(x,abs4log(mean(eigvecasy(:,1:n),2)),'^','MarkerSize',12);
            xlabel(xlab);
            ylabel('$\mid \bar{v_i} \mid$','Interpreter','Latex');
            legend('Analytic Eigenvectors Mean', 'Asymptotic Eigenvectors Mean');
            title({[name ' ' obj.scenarioName];obj.desc;figdesccd;'v_i^{(j)} = i^{th} element of j^{th} eigenvector'});
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
            save(fullfile(obj.resultsPath,[obj.scenarioName 'Obj.mat']),'obj','-v7.3');
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

