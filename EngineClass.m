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
        perturbations
        pert_eigvec_ana_1
        pert_eigvec_asy_1
        num_perturbations = 1;
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
        numbins
        bins
        binsij
        kbinned
        kijbinned
        Dii_anabinned
        Dii_asybinned
        Wij_anabinned
        Wij_asybinned
        numeigen=10;
        numeigenplots = 5;
        eigvec_dist_comparison_mat_ana2asy
        eigvec_angle_comparison_mat_ana2asy
        eigvec_dist_comparison_mat_ana2asy_permuted
        eigvec_angle_comparison_mat_ana2asy_permuted
        permutation_eigvec_ana2asy
        eps = [1];
        eps_adjusted;
        isInitsLegit;
        epsThreshold = .01;
        perturbation_factor = .1;
        
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
        steady_state
        steady_state_calculated
        mu = 0
        nu = 0
        rho = 0
        eta = 0
        isEngineSet = false
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
            obj.set_kinn();
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
            ww = [];
            for k = 2:obj.numbins+1
                if k == 2
                    ind = find(obj.degree_vector_weighted>=c^ns(k-1) & obj.degree_vector_weighted<=c^ns(k)+1e-13);
                else
                    ind = find(obj.degree_vector_weighted>c^ns(k-1) & obj.degree_vector_weighted<=c^ns(k)+1e-13);
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
                        if obj.adjacencyMatrix(i,j) > 0 && (kk > c1^ms(k-1) || (k==2 && kk>c1^ms(k-1)-1e-13)) && kk <=c1^ms(k) + 1e-13
                            obj.binsij{k-1}(end+1,:) = [i j];
                            w(end+1) = kk;
                        end
                    end
                end
                if ~isempty(w)
                    obj.kijbinned(k-1,1) = mean(w);
                end
                ww = [ww, w];
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
        function obj = set_perturbations(obj,isOnlyPositive)
            ss = obj.steady_state;
            k = obj.num_perturbations;
            n = size(ss,2);
            seed = obj.randSeed;
            rng(seed);
            obj.perturbations = {};
            obj.perturbations{1} = zeros(n,k);
            isSSzero = norm(ss)<obj.absTol;
            ss = ones(size(ss)).*max(isSSzero,ss);
            for j = 1:k
                percent = rand(1,n)*obj.perturbation_factor;
                sign = rand(1,n);
                sign = (sign > (~isOnlyPositive/2))*2 - 1;
                obj.perturbations{1}(:,j) = ss.*percent.*sign;
            end
%             numvec = 5;
%             obj.perturbations{2}(:,1) = sum(obj.eigenvectors_ana(:,1:numvec),2);
%             obj.perturbations{3}(:,1) = sum(obj.eigenvectors_asy_permuted(:,1:numvec),2);
        end
        function X = set_pert_eigvec_1(obj,eigenvectors)
            n = obj.numeigenplots;
            V = eigenvectors(:,1:n)';
            B = ones(n,1)/sqrt(n);
            X = linsolve(V,B);
        end
        function obj = solve(obj,pertType,epsIndStart,pertIndStart,isAdjustEpsilon)
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
                    sol_t_var_str = 'obj.solution_t';
                    sol_x_var_str = 'obj.solution_x';
                    stop_cond_str = 'size(obj.solution_t{1,pertInd,epsInd},1)>100 && max(abs(obj.solution_x{1,pertInd,epsInd}(end,:)-obj.solution_x{1,pertInd,epsInd}(end-100,:)))<obj.absTol';
                case 2 % perturbations are eigenvectors
                    perts = [obj.eigenvectors_ana(:,1:nn), obj.eigenvectors_asy_permuted(:,1:nn),...
                        sum(obj.eigenvectors_ana(:,1:nn),2),sum(obj.eigenvectors_asy_permuted(:,1:nn),2)...
                        obj.pert_eigvec_ana_1,obj.pert_eigvec_asy_1];
                    sol_t_var_str = 'obj.solution_t_eigvec';
                    sol_x_var_str = 'obj.solution_x_eigvec';
                    stop_cond_str = 'max(abs(obj.solution_x_eigvec{1,pertInd,epsInd}(end,:)-obj.steady_state))<obj.absTol*100';
                case 3 % perturbations are random and small relative to steady state
                    perts = obj.perturbations{1};
                    sol_t_var_str = 'obj.solution_t_perturbations';
                    sol_x_var_str = 'obj.solution_x_perturbations';
                    stop_cond_str = 'max(abs(obj.solution_x_perturbations{1,pertInd,epsInd}(end,:)-obj.steady_state))<obj.absTol*100';
                otherwise
            end
            numperts = size(perts,2);
            numeps = length(eps_vals);
            eval([sol_t_var_str ' = {};']);
            eval([sol_x_var_str ' = {};']);
            for epsInd = epsIndStart:numeps
                for pertInd = pertIndStart:numperts
                    eps_val = eps_vals(epsInd);
                    obj.isInitsLegit(epsInd,pertInd) = false;
                    while ~obj.isInitsLegit(epsInd,pertInd)
                        disp(['eps = ' num2str(eps_val)]);
                        if eps_val < obj.epsThreshold
                            disp('eps too small, breaking');
                            break
                        end
                        pert = eps_val * perts(:,pertInd);
                        init = ss' + pert;
                        if isAdjustEpsilon
                            if any(init<0 | init > 1,'all')
                                disp('Adjusting Epsilon');
                                eps_val = eps_val*.9;
                                disp(['eps = ' num2str(eps_val)]);
                            else
                                obj.isInitsLegit(epsInd,pertInd) = true;
                            end
                        else
                            obj.isInitsLegit(epsInd,pertInd) = true;
                        end
                    end
                    obj.eps_adjusted(epsInd,pertInd) = eps_val;
                    disp(['Final eps= ' num2str(eps_val)]);
                    inits{1,epsInd}(:,pertInd) = init;
                    init = init';
                    
                    eval([sol_t_var_str '{1,pertInd,epsInd} = 0;']);
                    eval([sol_x_var_str '{1,pertInd,epsInd} = init;']);
                end
            end
            opts = odeset('RelTol',obj.relTol,'AbsTol',obj.absTol);
            clear odefun;
            odefun = @(tt,x) (obj.f_M0(x) + (obj.adjacencyMatrix*obj.f_M2(x)).*obj.f_M1(x));
            for epsInd = epsIndStart:numeps
                disp(['epsInd = ' num2str(epsInd)]);
                disp(['eps_var = ' num2str(obj.eps_adjusted(epsInd))]);
                for pertInd = pertIndStart:numperts
                    if obj.isInitsLegit(epsInd,pertInd)
                        disp(['pertInd = ' num2str(pertInd)]);
                        init = inits{1,epsInd}(:,pertInd);
                        t = 0;
                        %%%%%%%%%
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
                            eval([sol_t_var_str '{1,pertInd,epsInd}(end+1:end+length(sol_t)-1,1)=sol_t(2:end,:);']);
                            eval([sol_x_var_str '{1,pertInd,epsInd}(end+1:end+size(sol_x,1)-1,:)=sol_x(2:end,:);']);
                            eval(['setBreak = ' stop_cond_str ';']);
                            if setBreak
                                disp('setBreak = true');
                            end
                        end
                    else
                        disp(['epsInd = ' num2str(epsInd) ', pertInd = ' num2str(pertInd)]);
                        disp(['Did not solve, no valid inits' ]);
                    end
                end
            end
            toc;
            if pertType == 1
                obj.solution_t = obj.solution_t{1};
                obj.solution_x = obj.solution_x{1};
                obj.set_steady_state();
                obj.set_M2_i_bigodot();
                obj.set_steady_state_calculated();
                obj.set_Dii_ana();
                obj.set_Wij_ana();
                obj.set_eig_ana();
                obj.set_Dii_anabinned();
                obj.set_Wij_anabinned();
                obj.set_eigvec_comparison_mats(true,false);
                obj.set_permutation_eigvec_ana2asy(40);
                obj.set_eig_asy_permuted();
                obj.set_eigvec_comparison_mats(false,true);
                obj.pert_eigvec_ana_1 = obj.set_pert_eigvec_1(obj.eigenvectors_ana);
                obj.pert_eigvec_asy_1 = obj.set_pert_eigvec_1(obj.eigenvectors_asy_permuted);
                if norm(obj.steady_state) < obj.absTol
                    isOnlyPositive = true;
                else
                    isOnlyPositive = false;
                end
                obj.set_perturbations(isOnlyPositive);
            elseif pertType==2
                obj.split_solution_eigvec();
            end
            obj.save_obj();
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
            disp('set_eigvec_comparison_mats');
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
        function obj = set_permutation_eigvec_ana2asy(obj,thresh)
            NN = obj.eigvec_angle_comparison_mat_ana2asy;
            NN = abs(90 - NN);
            [M,I] = max(NN,[],2);
            M = M < thresh;
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
            ss = obj.steady_state';
            t1 = obj.solution_t;ind_t1 = 1:step:length(t1);t1 = t1(ind_t1);
            x1 = obj.solution_x(ind_t1,:); p1 = (x1'-ss)';
            t2 = obj.solution_t_perturbations{1};ind_t2 = 1:step:length(t2);
            t2 = t2(ind_t2);x2 = obj.solution_x_perturbations{1};x2 = x2(ind_t2,:);
            p2 = (x2'-ss)';
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
                        if obj.isInitsLegit(epsInd,i)
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
                        if obj.isInitsLegit(epsInd,n+i)
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
                            if ~isempty(figdesc{figInd1,figInd}) && ~isempty(obj.solution_t_eigvecana{1,i,epsInd})
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
                                        if obj.isInitsLegit(epsInd,i)
                                            epsStr = num2str(obj.eps_adjusted(epsInd,i));
                                            myPlot{myPlotInd}(obj.solution_t_eigvecana{1,i,epsInd}(1:step:end,1),angles_ana{figInd1,i,figInd}(1:step:end,1),'.','Color',CM(i,:),'MarkerSize',12);
                                            hold on;
                                            legendStr{end+1} = ['$x_0 = ' ssstr epsStr ' * v_{' num2str(i) ',ana}$'];
                                        end
                                    end
                                    for i = 1:n
                                        if obj.isInitsLegit(n+i)
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
                        if figInd == 4 && myPlotInd == 1
                            myPlotInd = 2;
                        else
                            break
                        end
                    end
                end
            end
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
            a = obj.adjacencyMatrix;
            a_uw = (a>0); %uw = unweighted
            g_uw = graph(a_uw);
            p = {};
            
            p{1} = plot(g);
            for i = 1:obj.numbins
                inds = obj.bins;
                
            end
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

