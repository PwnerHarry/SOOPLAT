% The pseudo global configuration of SOOPLAT
% The uniform Platform for Evolutionary Single-objective Optimization
classdef GLOBAL < handle
    properties (SetAccess=private, GetAccess=public)
        algorithm % name of the algorithm
        evaluated % FE consumption counter
        bestIndividual % the current best individual, should be updated by the method evaluate()
        bestFitness % the current best fitness, should be updated by the method evaluate()
        trace % updates whenever a better bestIndividual is found, should be updated by the method evaluate()
        runtime % total runtime
        evaluation % total FEs
        identifier % name
        problem
    end
    properties (SetAccess=public, GetAccess=public)
        additionals
    end
    methods
        function obj = GLOBAL(varargin)
            additionals = [];
            arginStrings = {'additionals', 'dimension', 'evaluation', 'problem', 'algorithm'};
            IsString = find(cellfun(@ischar, varargin));
            [~, Loc] = ismember(varargin(IsString), cellfun(@(S)['-', S], arginStrings, 'UniformOutput', false));
            for i = find(Loc ~= 0)
                eval([arginStrings{Loc(i)}, ' = varargin{', num2str(IsString(i) + 1), '};']);
            end
            obj.algorithm = algorithm;
            obj.trace = [0, inf];
            obj.problem = feval(problem, dimension);
            obj.evaluated = 0;
            obj.evaluation = evaluation;
            obj.runtime.start = datetime('now');
            obj.runtime.end = [];
            obj.runtime.total = [];
            obj.identifier = getIdentifier();
            obj.bestFitness = inf;
            obj.additionals = additionals;
        end
        function fit = evaluateIndividuals(obj, Individuals)
            fit = evaluate(obj, cat(1, Individuals.solution));
            for i = 1: numel(Individuals)
                Individuals(i).evaluated = true;
                Individuals(i).fitness = fit(i);
            end
        end
        function reduceTrace(obj, max_points)
            % reduce the trace to evenly sampled curve to reduce size
            % TODO: add block reduce feature to lower the requirement for
            % memory
            ratios = linspace(0, 1, max_points)';
            D = pdist2(obj.trace(:, 1), ratios);
            sample_index = [];
            for i = 1: max_points
                [~, I] = min(D(:, i));
                sample_index = [sample_index, I];
            end
            obj.trace = obj.trace(sample_index, :);
        end
        function flag = terminated(obj)
            if obj.evaluated >= obj.evaluation || obj.bestFitness <= obj.problem.idealfitness
                flag = true;
            else
                flag = false;
            end
        end
        function fit = evaluate(obj, Population)
            if istable(Population)
                Population = table2array(Population);
            end
            [~, ~, violated_index] = obj.boundaryCheck(Population);
            refinedPopulation = Population(violated_index == 0, :);
            prefit = feval(obj.problem.functionhandle, refinedPopulation);
            fit = NaN(size(Population, 1), 1);
            fit(violated_index == 0) = prefit;
            fit(fit <= obj.problem.idealfitness) = 0;
            obj.evaluated = obj.evaluated + length(prefit);
            if obj.evaluated >= obj.evaluation && obj.terminated == false
                % obj.runtime.end = datetime('now');
                % obj.runtime.total = obj.runtime.end - obj.runtime.start;
                if obj.trace(end, 1) ~= 1
                    obj.trace = [obj.trace; [1, obj.trace(end, 2)]];
                end
            end
            [minfit, id] = min(fit);
            if isnan(obj.bestFitness) || minfit < obj.bestFitness
                obj.bestFitness = minfit;
                obj.bestIndividual = Population(id, :);
                obj.trace = [obj.trace; [obj.evaluated / obj.evaluation, obj.bestFitness]];
            end
        end
        function [less_pos, greater_pos, violated_index] = boundaryCheck(obj, Population)
            NP = size(Population, 1);
            xl = repmat(obj.problem.lowerbound, NP, 1);% check the lower bound
            less_pos = Population < xl;
            xu = repmat(obj.problem.upperbound, NP, 1);% check the upper bound
            greater_pos = Population > xu;
            violated_index = any(less_pos + greater_pos, 2)';
        end
        function draw(obj)
            warning('GLOBAL.draw() is deprecated, please use renderCurve()');
            if isempty(obj.trace)
                return;
            end
            persistent lastCurve background oldProblem;
            if ~isempty(oldProblem)
                if ~strcmp(oldProblem, obj.problem.name)
                    lastCurve = [];
                    background = [];
                end
            end
            if obj.terminated
                if isempty(lastCurve)
                    return;
                end
                hold on;
                delete(lastCurve);
                background(numel(background) + 1).curve = semilogy(obj.trace(:, 1), obj.trace(:, 2) - obj.problem.idealfitness, 'LineWidth', 2);
                if isempty(obj.algorithm)
                    algoname = "unknown";
                else
                    algoname = string(obj.algorithm);
                end
                background(end).algoname = algoname;
                legend(cat(2, background.curve), mat2cell(cat(1, background.algoname), ones(1, numel(background)))');
                %set(background(end).legend, 'Fontname', 'Book Antiqua', 'FontWeight', 'bold', 'FontSize', 12);
                axis([0 1 0 inf]);
                lastCurve = [];
                set(gca, 'YScale', 'log');
                drawnow;
                return;
            end
            if ~isempty(lastCurve)
                hold on;
                delete(lastCurve);
                lastCurve = semilogy(obj.trace(:, 1), obj.trace(:, 2) - obj.problem.idealfitness, 'LineWidth', 2);
                hold off;
            else
                hold on;
                lastCurve = semilogy(obj.trace(:, 1), obj.trace(:, 2) - obj.problem.idealfitness, 'LineWidth', 2);
            end
            axis([0 1 0 inf]);
            if isempty(obj.algorithm)
                algoname = "unknown";
            else
                algoname = string(obj.algorithm);
            end
            currentLegend = sprintf("%s %.2f%% (%.2e)", algoname, 100 * obj.evaluated / obj.evaluation, obj.bestFitness);
            if isempty(background)
                L = legend(lastCurve, currentLegend);
            else
                L = legend([cat(2, background.curve), lastCurve], mat2cell([cat(1, background.algoname); currentLegend], ones(1, numel(background) + 1))');
            end
            set(L, 'Fontname', 'Book Antiqua', 'FontWeight', 'bold', 'FontSize', 12);
            set(gca, 'YScale', 'log');
            drawnow;
            oldProblem = obj.problem.name;
        end
        function msdata = mileStone(obj, milestones)
            j = 1;
            msdata = NaN(size(milestones));
            for i = 1: numel(milestones)
                while milestones(i) / obj.evaluation > obj.trace(j, 1) && j < size(obj.trace, 1)
                    j = j + 1;
                end
                msdata(i) = obj.trace(j, 2);
            end
            msdata = [milestones; msdata]';
        end
    end
end
