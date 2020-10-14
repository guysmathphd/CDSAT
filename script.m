myKeys = {'scenarioName','adjacencyMatrix','initialValues','maxTime','maxDerivative','solverTimeStep','randSeed','f_M0','f_M1','f_M2','difEqSolver','absTol','relTol','resultsPath'};
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
display(eng)
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
display(eng)
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
display(eng)
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
display(eng)
disp('now solving');
SIS.solve();
SIS.print_output();

%% test6SIS5 4 nodes, node 1 and 2 connected, 3 and 4 connected
name = 'test6SIS5';
seed = 1;
rng(seed);
A = [1 1 0 0; 1 1 0 0; 0 0 1 1; 0 0 1 1];
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
display(eng)
disp('now solving');
SIS.solve();
SIS.print_output();

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
display(eng)
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
display(eng)
disp('now solving');
SIS.solve();
SIS.print_output();


