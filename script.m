myKeys = {'scenarioName','adjacencyMatrix','initialValues','maxTime','maxDerivative','solverTimeStep','randSeed','f_M0','f_M1','f_M2','difEqSolver','absTol','relTol','resultsPath'};
%%
myKeys = {'scenarioName','desc','adjacencyMatrix','initialValues','maxTime','maxDerivative','solverTimeStep','randSeed','f_M0','f_M1','f_M2','difEqSolver','absTol','relTol','resultsPath'};
%%
myKeys = {'scenarioName','desc','adjacencyMatrix','initialValues','maxTime','maxDerivative','solverTimeStep','randSeed','f_M0','f_M1','f_M2','difEqSolver','absTol','relTol','resultsPath','mu','nu','rho','eta'};
%%
myKeys = {'scenarioName','desc','adjacencyMatrix','initialValues','maxTime','maxDerivative','solverTimeStep','randSeed','f_M0','f_M1','f_M2','difEqSolver','absTol','relTol','resultsPath','mu','nu','rho','eta','numbins'};
%%
myKeys = {'scenarioName','desc','adjacencyMatrix','initialValues','maxTime','maxDerivative','solverTimeStep','randSeed','f_M0','f_M1','f_M2','difEqSolver','absTol','relTol','resultsPath','mu','nu','rho','eta','numbins','numeigen','isEngineSet'};
%%
myKeys = {'scenarioName','desc','adjacencyMatrix','initialValues','maxTime','maxDerivative','solverTimeStep','randSeed','f_M0','f_M1','f_M2','difEqSolver','absTol','relTol','resultsPath','mu','nu','rho','eta','numbins','numeigen','isEngineSet','init_condition','stop_condition_2'};
%%
myKeys = {'scenarioName','desc','adjacencyMatrix','initialValues','maxTime','maxDerivative','solverTimeStep','randSeed','f_M0','f_M1','f_M2','difEqSolver','absTol','relTol','resultsPath','mu','nu','rho','eta','numbins','numeigen','isEngineSet','init_condition_str','stop_condition_2'};
%% 
% A simple system with two independent nodes whose solution is x = Cexp(-t)
A = eye(2);
x0 = [5; -5];
maxT = 50;
maxDt = 0.01;
timeStep = 10;
seed = 1;
m0 = @(x) (0);
m1 = @(x) -x;
m2 = @(x) [1;1];
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
resultsPath = fullfile('tests','test1','results');
myValues = {'test1',A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath};
props = containers.Map(myKeys,myValues);
eng = EngineClass(props);
display(eng)
disp('now solving');
eng.solve();
eng.print_output();
eng.plot_results(false);
eng.plot_steady_vs_degree();

%% test2SIS1 - one equation
A = 1;
x0 = .1;
maxT = 50;
maxDt = 0.01;
timeStep = 10;
seed = 1;
m0 = @(x) (-2*x);
m1 = @(x) (1-x);
m2 = @(x) (4*x);
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
resultsPath = fullfile('tests','test2SIS1','results');
myValues = {'test2SIS',A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
disp('now solving');
SIS.solve();
SIS.print_output();
%% test3SIS2 - a two node dependent system
name = 'test3SIS2';
A = ones(2,2);
x0 = [.1; 5];
maxT = 50;
maxDt = 0.01;
timeStep = 10;
seed = 1;
m0 = @(x) (-2*x);
m1 = @(x) (1-x);
m2 = @(x) (4*x);
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
resultsPath = fullfile('tests',name,'results');
myValues = {name,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
disp('now solving');
SIS.solve();
SIS.print_output();

%% test4SIS3 10 nodes, non-symmetric, edges uniformly distributed
name = 'test4SIS3';
seed = 1;
rng(seed);
A = randi([0 1], 10,10);
x0 = rand(1,10)';
maxT = 50;
maxDt = 0.01;
timeStep = 10;
m0 = @(x) (-2*x);
m1 = @(x) (1-x);
m2 = @(x) (4*x);
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
resultsPath = fullfile('tests',name,'results');
myValues = {name,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
disp('now solving');
SIS.solve();
SIS.print_output();

%% test5SIS4 10 nodes, symmetric A, edges uniformly distributed
name = 'test5SIS4';
seed = 1;
rng(seed);
A = randi([0 1], 10,10);
A = triu(A,1) + tril(A');
x0 = rand(1,10)';
maxT = 50;
maxDt = 0.01;
timeStep = 10;
m0 = @(x) (-2*x);
m1 = @(x) (1-x);
m2 = @(x) (4*x);
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
resultsPath = fullfile('tests',name,'results');
myValues = {name,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
disp('now solving');
SIS.solve();
SIS.print_output();

%% test6SIS5 4 nodes, node 1 and 2 connected, 3 and 4 connected
name = 'test6SIS5';
seed = 1;
rng(seed);
A = [0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0];
x0 = [.1, .9, 1.5, -.5]';
maxT = 50;
maxDt = 0.01;
timeStep = 10;
m0 = @(x) (-2*x);
m1 = @(x) (1-x);
m2 = @(x) (4*x);
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
resultsPath = fullfile('tests',name,'results');
myValues = {name,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
%display(SIS)
disp('now solving');
SIS.solve();
SIS.print_output();
SIS.plot_results(false);
SIS.plot_steady_vs_degree;

%% test7VolterraLotka1
name = 'test7VolterraLotka';
seed = 1;
rng(seed);
A = [0 -2; 1 0];
x0 = [1; 3];
maxT = 50;
maxDt = 0.01;
timeStep = 10;
m0 = @(x) ([2*x(1); -x(2)]);
m1 = @(x) (x);
m2 = @(x) (x);
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
resultsPath = fullfile('tests',name,'results');
myValues = {name,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
disp('now solving');
SIS.solve();
SIS.print_output();

%% test8SIS6 1000 nodes, symmetric A, edges uniformly distributed
name = 'test8SIS6';
seed = 1;
rng(seed);
A = randi([0 1], 1000,1000);
A = triu(A,1) + tril(A');
x0 = rand(1,1000)';
maxT = 50;
maxDt = 0.01;
timeStep = 10;
m0 = @(x) (-2*x);
m1 = @(x) (1-x);
m2 = @(x) (4*x);
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
resultsPath = fullfile('tests',name,'results');
myValues = {name,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
disp('now solving');
SIS.solve();
SIS.print_output();
SIS.plot_results(true);
SIS.plot_steady_vs_degree();
SIS.calculate_degree_weighted();
%% %% test9SIS7 - a two node dependent system with no self loops
name = 'test9SIS7';
A = ones(2,2)-eye(2);
x0 = [.1; 5];
maxT = 50;
maxDt = 0.01;
timeStep = 10;
seed = 1;
m0 = @(x) (-2*x);
m1 = @(x) (1-x);
m2 = @(x) (4*x);
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
resultsPath = fullfile('tests',name,'results');
myValues = {name,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
disp('now solving');
SIS.solve();
SIS.print_output();
SIS.plot_results(false);
SIS.plot_steady_vs_degree;

%% test10SIS8 10 nodes, symmetric A, edges uniformly distributed
name = 'test10SIS8';
desc = '10 node random';
seed = 1;
rng(seed);
mask = ones(10,10) - eye(10);
A = randi([0 1], 10,10);
A = triu(A,1) + tril(A');
A = A.*mask;
x0 = rand(1,10)';
maxT = 50;
maxDt = 0.01;
timeStep = 10;
m0 = @(x) (-2*x);
m1 = @(x) (1-x);
m2 = @(x) (4*x);
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
resultsPath = fullfile('tests',name,'results');
myValues = {name,desc,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
disp('now solving');
SIS.solve();
SIS.print_output();
SIS.plot_results(false);
SIS.plot_steady_vs_degree();
% SIS.set_M2_i_bigodot;
% SIS.set_steady_state_calculated;

%% test11SIS9 Same as test10SIS8, A loaded from file
name = 'test11SIS9';
 seed = 1;
 rng(seed);
mask = ones(10,10) - eye(10);
A = randi([0 1], 10,10);
A = triu(A,1) + tril(A');
A = A.*mask;
load(fullfile('networks','A8.mat'));
x0 = rand(1,10)';
maxT = 50;
maxDt = 0.01;
timeStep = 10;
m0 = @(x) (-2*x);
m1 = @(x) (1-x);
m2 = @(x) (4*x);
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
resultsPath = fullfile('tests',name,'results');
myValues = {name,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
disp('now solving');
SIS.solve();
SIS.print_output();
%% %% test12SIS10 A = SF1
name = 'test12SIS10';
desc = 'SF1';
seed = 1;
rng(seed);
load(fullfile('networks','SF1.mat'));
x0 = rand(1,6000)';
maxT = 50;
maxDt = 0.01;
timeStep = 10;
m0 = @(x) (-2*x);
m1 = @(x) (1-x);
m2 = @(x) (4*x);
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
resultsPath = fullfile('tests',name,'results');
myValues = {name,desc,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
disp('now solving');
SIS.solve();
SIS.print_output();
SIS.plot_results(true);
SIS.calculate_degree();
SIS.calculate_degree_weighted();
SIS.set_steady_state();
SIS.plot_steady_vs_degree();
SIS.set_dM0;
SIS.set_dM1;
SIS.set_dM2;
SIS.set_const_functions();
SIS.set_Dii();
SIS.set_Wij;
SIS.mu = 1; SIS.nu = -1; SIS.rho = 0; SIS.eta = 0;
SIS.plot_jacobian;

%% %% %% test13SIS11 A = SF1>0 (SF1 projected onto {0,1})
name = 'test13SIS11';
seed = 1;
rng(seed);
load(fullfile('networks','SF1.mat'));
A = double(A>0);
x0 = rand(1,6000)';
maxT = 50;
maxDt = 0.01;
timeStep = 10;
m0 = @(x) (-2*x);
m1 = @(x) (1-x);
m2 = @(x) (4*x);
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
resultsPath = fullfile('tests',name,'results');
myValues = {name,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
SIS.set_R;
disp('now solving');
SIS.solve();
SIS.print_output();
SIS.plot_results(true);
SIS.calculate_degree();
SIS.calculate_degree_weighted();
SIS.set_steady_state();
SIS.mu = 1; SIS.nu = -1; SIS.rho = 0; SIS.eta = 0;
SIS.plot_steady_vs_degree();
%% %% %% test14SIS12 A = SF2
name = 'test14SIS12';
seed = 1;
rng(seed);
load(fullfile('networks','SF2.mat'));
x0 = rand(1,6000)';
maxT = 50;
maxDt = 0.01;
timeStep = 10;
m0 = @(x) (-2*x);
m1 = @(x) (1-x);
m2 = @(x) (4*x);
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
resultsPath = fullfile('tests',name,'results');
myValues = {name,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
disp('now solving');
SIS.solve();
SIS.print_output();
SIS.plot_results(true);
SIS.plot_steady_vs_degree();
%% test15SIS13 3 nodes, node 1 and 2, 2 and 3, 3 and 1 connected
name = 'test15SIS13';
desc = '3 nodes fully connected';
seed = 1;
rng(seed);
A = ones(3,3) - eye(3);
x0 = [.1, .9, 1.5]';
maxT = 50;
maxDt = 0.01;
timeStep = 10;
m0 = @(x) (-2*x);
m1 = @(x) (1-x);
m2 = @(x) (4*x);
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
resultsPath = fullfile('tests',name,'results');
myValues = {name,desc,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
%display(SIS)
disp('now solving');
SIS.solve();
SIS.print_output();
SIS.plot_results(false);
SIS.plot_steady_vs_degree;
SIS.plot_jacobian;
%% %% %% test16REG1 A = SF1
name = 'test16REG1';
desc = 'SF1';
seed = 1;
rng(seed);
load(fullfile('networks','SF1.mat'));
x0 = rand(1,6000)';
maxT = 50;
maxDt = 0.01;
timeStep = 10;
m0 = @(x) (-x);
m1 = @() (1);
m2 = @(x) (x./(1+x));
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
resultsPath = fullfile('tests',name,'results');
myValues = {name,desc,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath};
props = containers.Map(myKeys,myValues);
REG = EngineClass(props);
%display(SIS)
disp('now solving');
REG.solve();
REG.print_output();
REG.plot_results(true);
REG.plot_steady_vs_degree();
%REG.mu = 1; REG.nu = -1; REG.rho = 0; REG.eta = 0;
REG.save_obj();
REG.plot_jacobian;
REG.set_Dii_ana();
REG.set_Wij_ana();
REG.set_eig_ana();
REG.set_Dii_anabinned();
REG.set_Wij_anabinned();
REG.mu=0;REG.nu=0;REG.rho=-2;REG.eta=0;
REG.numbins = 15;
            REG.set_N();
            REG.set_knn();
            REG.set_Dii_asy();
            REG.set_Wij_asy();
            REG.set_eig_asy();
            REG.set_bins();
            REG.set_kbinned();
            REG.set_Dii_asybinned();
            REG.set_Wij_asybinned();
REG.set_eigvec_comparison_mats(true,false);
REG.set_permutation_eigvec_ana2asy();
REG.set_eig_asy_permuted();
REG.set_eigvec_comparison_mats(false,true);
REG.solve(true,false,true,1,1,false);
REG.solve(false,true,true,1,1,false);
%% test17SIS14 A = SF1 absTol, relTol = 1e-12
name = 'test17SIS14';
desc = 'SF1';
seed = 1;
rng(seed);
load(fullfile('networks','SF1.mat'));
x0 = rand(1,6000)';
maxT = 50;
maxDt = 0.01;
timeStep = 10;
m0 = @(x) (-2*x);
m1 = @(x) (1-x);
m2 = @(x) (4*x);
difEqSolver = @ode45;
absTol = 1e-12;
relTol = 1e-12;
mu = 1;
nu = -1;
rho = 0;
eta = 0;
resultsPath = fullfile('tests',name,'results');
myValues = {name,desc,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath,mu,nu,rho,eta};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
disp('now solving');
SIS.solve();
SIS.print_output();
SIS.plot_results(true);
SIS.plot_steady_vs_degree();
%SIS.mu = 1; SIS.nu = -1; SIS.rho = 0; SIS.eta = 0;
SIS.plot_jacobian;
%             SIS.set_N();
%             SIS.set_knn();
%             SIS.set_Dii_asy();
%             SIS.set_Wij_asy();
%             SIS.save_obj();
% SIS.set_eig_ana();
% SIS.set_eig_asy();
SIS.plot_eigenvalues();
SIS.plot_eigenvectors();
% SIS.save_obj();
%% test18 knn calculation
name = 'test18';
desc = 'knn_calculation';
seed = 1;
rng(seed);
A = eye(3);
x0 = rand(1,3)';
maxT = 50;
maxDt = 0.01;
timeStep = 10;
m0 = @(x) (-2*x);
m1 = @(x) (1-x);
m2 = @(x) (4*x);
difEqSolver = @ode45;
absTol = 1e-12;
relTol = 1e-12;
mu = 1;
nu = 1;
rho = 1;
eta = 1;
resultsPath = fullfile('tests',name,'results');
myValues = {name,desc,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath,mu,nu,rho,eta};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
%% test19 engine set knn calculation
name = 'test19';
desc = '';
seed = 1;
rng(seed);
A = ones(3)-eye(3);
x0 = rand(1,3)';
maxT = 50;
maxDt = 0.01;
timeStep = 10;
m0 = @(x) (-2*x);
m1 = @() (1);
m2 = @(x) (x.^2);
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
mu = 0;
nu = 0;
rho = 2;
eta = 4;
resultsPath = fullfile('tests','engset1',name,'results');
myValues = {name,desc,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath,mu,nu,rho,eta};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
disp('now solving');
SIS.solve();
SIS.print_output();
SIS.plot_results(false);
SIS.plot_steady_vs_degree();
SIS.plot_jacobian;
%% test20 engine set knn calculation
name = 'test20';
n = 100;
desc = '';
seed = 1;
rng(seed);
A = ones(n)-eye(n);
x0 = rand(1,n)';
maxT = 50;
maxDt = 0.01;
timeStep = 10;
m0 = @(x) (-(n-1)*x);
m1 = @() (1);
m2 = @(x) (x.^2);
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
mu = 0;
nu = 0;
rho = 2;
eta = 4;
resultsPath = fullfile('tests','engset1',name,'results');
myValues = {name,desc,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath,mu,nu,rho,eta};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
disp('now solving');
SIS.solve();
SIS.print_output();
SIS.plot_results(false);
SIS.plot_steady_vs_degree();
SIS.plot_jacobian;
%% test21 J sanity test
name = 'test21';
n = 6000;
desc = '';
seed = 1;
rng(seed);
load(fullfile('networks','SF1.mat'));
x0 = rand(1,n)';
maxT = 5;
maxDt = 0.01;
timeStep = 10;
m0 = @(x) (-(n-1)*x);
m1 = @() (1);
m2 = @(x) (x.^2);
difEqSolver = @ode45;
absTol = 1e-2;
relTol = 1e-2;
mu = 0;
nu = 0;
rho = 2;
eta = 4;
resultsPath = fullfile('tests','engset1',name,'results');
myValues = {name,desc,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath,mu,nu,rho,eta};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
disp('now solving');
SIS.solve();
SIS.print_output();
SIS.plot_results(false);
SIS.plot_steady_vs_degree();
SIS.plot_jacobian;
%% test makePropertiesMaps
myKeys = {'key1','key2','key3'};
v1 = num2cell(1:4);
v2 = {'key2value1','key2value22','key2value333'};
v3 = {eye(1),ones(2,2),zeros(3,3)};
myValues = {v1,v2,v3};
propMapsDefs = containers.Map(myKeys,myValues);
propertiesMaps=makePropertiesMaps(propMapsDefs);
%% test22SIS15 A = SF1 absTol, relTol = 1e-5
name = 'test22SIS15';
desc = 'SF1';
seed = 1;
rng(seed);
load(fullfile('networks','SF1.mat'));
x0 = rand(1,6000)';
maxT = 50;
maxDt = 0.01;
timeStep = 10;
m0 = @(x) (-2*x);
m1 = @(x) (1-x);
m2 = @(x) (4*x);
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
mu = 1;
nu = -1;
rho = 0;
eta = 0;
numbins = 15;
resultsPath = fullfile('tests',name,'results');
myValues = {name,desc,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath,mu,nu,rho,eta,numbins};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
disp('now solving');
SIS.numeigen = 50;
SIS.eps = [.1 .5];
SIS.maxTime = 5;
SIS.solve(false,false,true,1,1,false);
SIS.solve(true,false,true,3,1,false);
obj.solve(false,true,true,3,2,false);
% SIS.print_output();
% SIS.plot_results(true);
% SIS.plot_steady_vs_degree();
%SIS.mu = 1; SIS.nu = -1; SIS.rho = 0; SIS.eta = 0;
% SIS.plot_jacobian;
%             SIS.set_N();
%             SIS.set_knn();
%             SIS.set_Dii_asy();
%             SIS.set_Wij_asy();
%             SIS.save_obj();
% SIS.set_eig_ana();
% SIS.set_eig_asy();
% SIS.plot_eigenvalues();
% SIS.plot_eigenvectors();
% SIS.save_obj();
%% test23SIS16 Eigenvalues Check
isEngineSet = {true};
f = [.1, 1, 5, 10, 50,100];
g = f;
ind = 1;
name = {'test23SIS16'};
m0={};
m2={};
str1 = 'm0{ind} = @(x) (-';
str2 = 'm2{ind} = @(x) (';
for kf=1:length(f)
    eval([str1 num2str(f(kf)) '*x)']);
    eval([str2 num2str(g(kf)) '*x)']);
    ind = ind+1;
end
desc = {'SF1'};
seed = {1};
rng(seed{1});
load(fullfile('networks','SF1.mat'));
x0 = {rand(1,6000)'};
maxT = {10};
maxDt = {0.01};
timeStep = {10};
m1 = {@(x) (1-x)};
difEqSolver = {@ode45};
absTol = {1e-5};
relTol = {1e-5};
mu = {1};
nu = {-1};
rho = {0};
eta = {0};
numbins = {15};
numeigen = {200};
resultsPath={};
resultsPath{1} = fullfile('tests','test23SIS16set1');
myValues = {name,desc,{A},x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath,mu,nu,rho,eta,numbins,numeigen,isEngineSet};
propMapsDefs = containers.Map(myKeys,myValues);
propertiesMaps = makePropertiesMaps(propMapsDefs);
% SISset = EngineSet(propertiesMaps);
path = 'tests\test23SIS16set1';
name = 'test23SIS16set1';
SISset = EngineSet(path,name);
SISset.set_engine_paths();
SISset.setSolve(32);
SISset.setPlot_results(32,true);
SISset.setFunction(@plot_eigenvalues,32);
SISset.set_eigenvalues(1,999);
%%%%%% Run these:
SISset.setSolve(1,false,true,1,1,true);
SISset.setSolve(1,true,true,1,1,true);
%%%%%%
SISset.plot_setEigenvalues();
%% test24ECO1 ecological dynamics B=K=1 type II mutualistic interactions
isEngineSet = false;
B = 1;
K = 1;
h = 1;
ind = 1;
name = 'test24ECO1';
m0 = @(x) (B*x).*(1-x/K);
m1 = @(x) (x);
m2 = @(x) x.^h ./ (1+x.^h);
desc = 'Ecological Dynamics on SF1';
seed = 1;
rng(seed);
load(fullfile('networks','SF1.mat'));
x0 = rand(1,6000)';
maxT = 1000;
maxDt = {0.01};
timeStep = .2;
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
mu = 1;
nu = 1;
rho = -2;
eta = 0;
numbins = 15;
numeigen = 200;
resultsPath= fullfile('tests',name,'results');
myValues = {name,desc,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath,mu,nu,rho,eta,numbins,numeigen,isEngineSet};
props = containers.Map(myKeys,myValues);
ECO = EngineClass(props);
ECO.solve(false,true,1,1,true);
ECO.solve(true,true,1,1,false);
ECO.plot_results(true);
ECO.plot_results2(true);
ECO.plot_steady_vs_degree();
ECO.plot_jacobian();
ECO.plot_eigenvalues();
ECO.plot_eigenvectors(true);
ECO.plot_eigenvectors2(false,false,true,true);
ECO.plot_eigenvectors2(true,false,true,true);
ECO.plot_eigenvectors2(true,true,true,true);
ECO.plot_random_states();
%propertiesMaps = makePropertiesMaps(propMapsDefs);
% SISset = EngineSet(propertiesMaps);
% path = 'tests\test23SIS16set1';
% name = 'test23SIS16set1';
% SISset = EngineSet(path,name);
% SISset.set_engine_paths();
% SISset.setSolve(32);
% SISset.setPlot_results(32,true);
% SISset.setFunction(@plot_eigenvalues,32);
% SISset.set_eigenvalues(1,999);
% %%%%%% Run these:
% SISset.setSolve(1,false,true,1,1,true);
% SISset.setSolve(1,true,true,1,1,true);
% %%%%%%
% SISset.plot_setEigenvalues();
%% test25BIO1 biochemical dynamics F=4, W=6, B=7, U=2, Q=3
isEngineSet = false;
F = 4;
W = 6;
B = 7;
U = 2;
Q = 3;
Btilde = 4.2; %QB/(U+Q) = 21/5 = 4.2
ind = 1;
name = 'test25BIO1';
m0 = @(x) (F - W*x);
m1 = @(x) (-Btilde*x);
m2 = @(x) (x);
desc = 'Biochemical Dynamics on SF1';
seed = 1;
rng(seed);
load(fullfile('networks','SF1.mat'));
x0 = rand(1,6000)';
maxT = 1000;
maxDt = {0.01};
timeStep = .2;
difEqSolver = @ode45;
absTol = 1e-5;
relTol = 1e-5;
mu = 1;
nu = -1;
rho = 0;
eta = -1;
numbins = 15;
numeigen = 200;
resultsPath= fullfile('tests',name,'results');
myValues = {name,desc,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath,mu,nu,rho,eta,numbins,numeigen,isEngineSet};
props = containers.Map(myKeys,myValues);
BIO = EngineClass(props);
BIO.solve(false,true,1,1,true);
BIO.solve(true,true,1,1,false);
BIO.plot_results(true);
BIO.plot_results2(true);
BIO.plot_steady_vs_degree();
BIO.plot_jacobian();
BIO.plot_eigenvalues();
BIO.plot_eigenvectors(true);
BIO.plot_eigenvectors2(false,false,true,true);
BIO.plot_eigenvectors2(true,false,true,true);
BIO.plot_eigenvectors2(true,true,true,true);
BIO.plot_random_states();
%%
sis.absTol = 1e-14;
sis.relTol = 1e-14;
sis.perturbation_factor = .01;
sis.solve(1,1,1,true);
sis.perturbation_factor = .001;
sis.set_perturbations(true);
sis.solve(3,1,1,true);
sis.plot_results2(true);
%%
obj.set_single_node_pert_struct(ones(6000,1));
%%

eps_vector = .1*obj.steady_state';
struct2 = set_single_node_pert_struct(obj,eps_vector,'single_node_perts_eps_0p1_ss');

%%
obj.plot_pert_approx(5631);
%%

[eigvals_asy_v2_set, eigvecs_asy_v2_set, Dii_asy_set, Wij_asy_set] = obj.set_eig_v2_sets(obj.Dii_asy, obj.Wij_asy, obj.C_D_set, obj.C_W_set, obj.numeigen,1);
%%
obj.plot_pert_approx(eigvecs_asy_v2_set,eigvals_asy_v2_set);
%%
[C_D, guesses1, guesses2, errors1, errors2, eigvals,eigvecs] = obj.find_best_C_D(obj.eigenvalues_ana(1,1),276);
obj.C_D_v3 = C_D;
obj.guesses1_v3 = guesses1;
obj.guesses2_v3 = guesses2;
obj.errors1_v3 = errors1;
obj.errors2_v3 = errors2;
obj.eigvals_v3 = eigvals;
obj.eigvecs_v3 = eigvecs;
%%
obj8.solve(2,1,1,1);
obj8.solve(3,1,1,1);
%%
obj.plot_pert_approx()
obj.plot_eigenvalues();
%%
obj1.solve(1,1,1,1);
obj1.solve(2,1,1,1);
obj1.solve(3,1,1,1);
%%
obj.save_obj();
clear obj;
% tests = {'test16REG1','test24ECO1','test25BIO1'};
tests = {'test16REG1','test24ECO1'};

for i=1:length(tests)
    disp(['i = ' num2str(i)]);
    test = tests{i};
    testobj = [test 'Obj.mat'];
    mypath = fullfile('tests',test,'results',testobj);
    load(mypath);
%     obj.numeigen = 6000;
%     [obj.kbinned,obj.kbinned_mins,obj.kbinned_maxs] = obj.set_binned_vals(obj.degree_vector_weighted,obj.bins);
%     obj.set_eig_asy();
%     obj.set_eig_ana();
%     obj.set_eigvec_comparison_mats(true,false);
%     obj.set_permutation_eigvec_ana2asy(0);
%     obj.set_eig_asy_permuted();
%     obj.set_eigvec_comparison_mats(false,true);
%     obj.set_weighted_dot_products();
%     obj.set_eig_ana_ordered_nodes();
%     obj.plot_eigenvectors(true);
%     obj.plot_eigenvectors4();

tol = 1e-13;
disp('obj.bins = obj.set_bins_generic');
obj.bins = obj.set_bins_generic(obj.numbins,obj.degree_vector_weighted,tol,true(obj.N,1));
disp('obj.kbinned = obj.set_binned_vals'); [obj.kbinned,obj.kbinned_mins,obj.kbinned_maxs] = obj.set_binned_vals(obj.degree_vector_weighted,obj.bins);
obj.eigenvectors_ana_binned_k = obj.set_binned_vals(obj.eigenvectors_ana,obj.bins);
obj.eigenvectors_ana_binned_kinn = obj.set_binned_vals(obj.eigenvectors_ana,obj.binskinn);
obj.eigenvectors_asy_permuted_binned_k = obj.set_binned_vals(obj.eigenvectors_asy_permuted,obj.bins);
obj.eigenvectors_asy_permuted_binned_kinn = obj.set_binned_vals(obj.eigenvectors_asy_permuted,obj.binskinn);
obj.plot_eigenvectors5();
    obj.save_obj();
    clear obj;
end
%% test16REG1
% obj.absTol = 1e-18;obj.relTol = 1e-20; obj.solverTimeStep =
% .2;obj.maxTime=1000;
obj.solve(1,1,1,1);
obj.plot_eigenvectors(true);
obj.plot_eigenvectors4();
    obj.plot_eigenvectors5();
    %%
    %% test26REG2 A = SF1 a=1/2
name = 'test26REG2';
desc = 'SF1';
seed = 1;
rng(seed);
load(fullfile('networks','SF1.mat'));
x0 = rand(1,6000)';
maxT = 1000;
maxDt = 0.01;
timeStep = .2;
m0 = @(x) (-x.^(1/2));
m1 = @() (1);
m2 = @(x) (x./(1+x));
difEqSolver = @ode45;
absTol = 1e-14;
relTol = 1e-14;
resultsPath = fullfile('tests',name,'results');
mu = -1; nu = 0; rho = -4; eta = 0;
myKeys = {'scenarioName','desc','adjacencyMatrix','initialValues','maxTime','maxDerivative','solverTimeStep','randSeed','f_M0','f_M1','f_M2','difEqSolver','absTol','relTol','resultsPath','mu','nu','rho','eta','numbins','numeigen','isEngineSet','init_condition','stop_condition_2'};
myValues = {name,desc,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath,mu,nu,rho,eta,15,6000,false,@(x) (true),[]};
props = containers.Map(myKeys,myValues);
REG = EngineClass(props);
%display(SIS)
disp('now solving');
REG.solve(1,1,1,1);
%%
% REG.print_output();
REG.plot_results(true);
REG.plot_steady_vs_degree();
REG.save_obj();
%%
obj.solve(1,1,1,2);
obj.solve_eigvec_pert_max_hub(2);
obj.save_obj();
nt = length(obj.solution_t);
ind = 1;step=50000;
for i = 1:step:nt
    str = ['solution_t_' num2str(ind)];ind = ind+1;
    endind = min(i + step-1,nt);
        obj.save_var(obj.solution_t(i:endind),fullfile(obj.resultsPath,'obj_properties'),'solution',str);
    start = i+1;
end
ind = 1;step=500;
for i = 1:step:nt
    str = ['solution_x_' num2str(ind)];ind = ind+1;
    endind = min(i + step-1,nt);
        obj.save_var(obj.solution_x(i:endind,:),fullfile(obj.resultsPath,'obj_properties'),'solution',str);
    start = i+1;
end
obj.solution_t = [];obj.solution_x = [];
obj1.save_obj();
obj.save_var(obj.solution_x,obj.resultsPath,'obj_properties','solution_x');

obj.save_var(sol_t_ana,fullfile(obj.resultsPath,'obj_properties'),'eigvec_pert_max_hub','sol_t_ana');
obj.save_var(sol_x_ana,fullfile(obj.resultsPath,'obj_properties'),'eigvec_pert_max_hub','sol_x_ana');

obj.load_solution('eigvec_pert_max_hub', 'sol_', '_ana', false, 'obj.solution_t_eigvecana{1,2,1}','obj.solution_x_eigvecana{1,2,1}');
obj.load_solution('eigvec_pert_max_hub', 'sol_', '_asy', false, 'obj.solution_t_eigvecasy{1,2,1}','obj.solution_x_eigvecasy{1,2,1}');
obj.isInitsLegit=[1,1,1,1,1,1,1];
obj.eps_adjusted=[.1,.1,.1,.1,.1,.1,.1];
obj.plot_results2();
%%
obj.save_obj();
obj.maxTime = 1000;
obj.solverTimeStep = .5;
obj.solve(2,1,1,2);
obj.solve_eigvec_pert_max_hub(1);
obj.solve_eigvec_pert_max_hub(2);
% obj.solve(3,1,1,2);

%%
for i=1:5
    sol_t = obj.solution_t_eigvecana{i};
    sol_x = obj.solution_x_eigvecana{i};
    obj.save_var(sol_t,obj.resultsPath,'obj_properties',['sol_t_ana_v' num2str(i)]);
    obj.save_var(sol_x,obj.resultsPath,'obj_properties',['sol_x_ana_v' num2str(i)]);
end
for i=1:5
    sol_t = obj.solution_t_eigvecasy{i};
    sol_x = obj.solution_x_eigvecasy{i};
    obj.save_var(sol_t,obj.resultsPath,'obj_properties',['sol_t_asy_v' num2str(i)]);
    obj.save_var(sol_x,obj.resultsPath,'obj_properties',['sol_x_asy_v' num2str(i)]);
end
%% EngineSet2 test27

name = 'test27set2';
folderNames = {'test23SIS16-8','test16REG1','test24ECO1','test25BIO1','test26REG2'};
legendNames = {'SIS','Regulatory a = 1','Ecological','Biochemical','Regulatory a = .5'};
obj = EngineSet2(folderNames,name);
obj.batchFunction(@solve_eigvec_pert_max_hub_1,[3 4]);
% obj.batchFunction(@set_graph_object,1:4);
obj.plotLocalization1();

%%
obj.initialValues = obj.steady_state';
obj.solve(1,1,1,2);
%%
alpha = -1:.01:1;
[X,Y] = meshgrid(alpha);
gamma = (((Y-X).^2).*X + 1 + X - Y)/(Y - X);
% surf(X,Y,gamma);
% surf(X,Y,sign(gamma));
figure;image(sign(gamma),'CDatamapping','scaled');colorbar;
%%
x = .1:.1:2;
y1 = 1./x;
y2 = 1./x - 1;
y3 = x;
y4 = 1 - x + x.^2 - x.^3 + x.^4 - x.^5 + x.^6 - x.^7;
y5 = 1./(x+1);
figure;plot(x,y1,x,y2,x,y3,x,y4,x,y5);axis equal;legend;
%%
SF1 = Network('SF1',A);
SF1.set_degree_perts();
%%
obj.networkName = 'SF1';
obj.solve_degree_weighted_perts();
%%
load('C:\49_CDSAT\tests\test27set2\test27set2Obj.mat');
 obj.batchFunction(@plot_localization3,[1:5]);
% obj.batchFunction(@solve_degree_weighted_perts,[1:5]);
obj.batchFunction(@solve_eigvec_pert_max_hub_1,[2 5]);
 obj.batchFunction(@plot_localization2,[5:5]);
 obj.batchFunction(@plot_eigenvalues2,[1:5]);
 obj.batchFunction(@solve_single_node_perts_batch,[1:5]);
 obj.batchFunction(@rename_files,[1:5]);
 obj.batchFunction(@plot_single_node_pert,[1:5]);