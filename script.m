myKeys = {'scenarioName','adjacencyMatrix','initialValues','maxTime','maxDerivative','solverTimeStep','randSeed','f_M0','f_M1','f_M2','difEqSolver','absTol','relTol','resultsPath'};
%%
myKeys = {'scenarioName','desc','adjacencyMatrix','initialValues','maxTime','maxDerivative','solverTimeStep','randSeed','f_M0','f_M1','f_M2','difEqSolver','absTol','relTol','resultsPath'};
%%
myKeys = {'scenarioName','desc','adjacencyMatrix','initialValues','maxTime','maxDerivative','solverTimeStep','randSeed','f_M0','f_M1','f_M2','difEqSolver','absTol','relTol','resultsPath','mu','nu','rho','eta'};
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
SIS.set_Dii;
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
REG.mu = 1; REG.nu = -1; REG.rho = 0; REG.eta = 0;
REG.save_obj();
REG.plot_jacobian;
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
resultsPath = fullfile('tests',name,'results');
myValues = {name,desc,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
disp('now solving');
SIS.solve();
SIS.print_output();
SIS.plot_results(true);
SIS.plot_steady_vs_degree();
SIS.mu = 1; SIS.nu = -1; SIS.rho = 0; SIS.eta = 0;
SIS.plot_jacobian;
            SIS.set_N();
            SIS.set_knn();
            SIS.set_Dii_asy();
            SIS.set_Wij_asy();
            SIS.save_obj();
%% test19 engine set knn calculation
name = 'test19';
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
resultsPath = fullfile('tests','engset1',name,'results');
myValues = {name,desc,A,x0,maxT,maxDt,timeStep,seed,m0,m1,m2,difEqSolver,absTol,relTol,resultsPath,mu,nu,rho,eta};
props = containers.Map(myKeys,myValues);
SIS = EngineClass(props);
display(SIS)
% disp('now solving');
% SIS.solve();
% SIS.print_output();
% SIS.plot_results(true);
% SIS.plot_steady_vs_degree();
% SIS.plot_jacobian;
%%%% test18 knn calculation
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
%% test makePropertiesMaps
myKeys = {'key1','key2','key3'};
v1 = num2cell(1:4);
v2 = {'key2value1','key2value22','key2value333'};
v3 = {eye(1),ones(2,2),zeros(3,3)};
myValues = {v1,v2,v3};
propMapsDefs = containers.Map(myKeys,myValues);
propertiesMaps=makePropertiesMaps(propMapsDefs);


