classdef EngineClass <  handle
    %ENGINECLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        scenarioName %String   
        adjacencyMatrix %2D array
        initialValues %1D array same length as adjacency matrix
        maxTime %positive integer
        maxDerivative % positive double
        solverTimeStep % positive integer
        randSeed %positive integer
        f_M0 %function handle
        f_M1 %function handle
        f_M2 %function handle
        difEqSolver %function handle
        absTol % positive double
        relTol % positive double
        solution_t = 0;
        solution_x ;
        resultsPath
        header
        
        
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
        end
        
        function obj = solve(obj)
            %Solve the system
            tic;
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
                 [sol_t,sol_x] = obj.difEqSolver(odefun,[t stepEndTime],init);
                 t = stepEndTime;
                 init = sol_x(end,:);
                 % append results to solution_t and solution_x
                 obj.solution_t = [obj.solution_t;sol_t(2:end,:)];
                 obj.solution_x = [obj.solution_x;sol_x(2:end,:)];
            end
            toc;
            % plot results
            f = figure;
            plot(obj.solution_t,obj.solution_x,'.-','MarkerSize',12);
            legend(obj.header(2:end));
            title(obj.scenarioName);
            xlabel(obj.header{1});
            ylabel('x');
            text(1,.5,num2str(obj.adjacencyMatrix));
            if ~isfolder(fullfile(obj.resultsPath,'figs'))
                mkdir(fullfile(obj.resultsPath,'figs'));
            end
            saveas(f,fullfile(obj.resultsPath,'figs','fig1.fig'),'fig');
                 
        end
        
        function obj = print_output(obj) %create output file
            T = array2table([obj.solution_t,obj.solution_x],'VariableNames',obj.header);
            if ~isfolder(obj.resultsPath)
                mkdir(obj.resultsPath)
            end
            writetable(T,fullfile(obj.resultsPath,'output.csv'));
        end
    end
end

