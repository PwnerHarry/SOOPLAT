function renderCurve(varargin)
if nargin == 1
    Global = varargin{1};
    draw_flag = true;
elseif nargin == 2
    Global = varargin{1};
    draw_flag = varargin{2};
end
if isempty(Global.trace)
    return;
end
persistent lastCurve background oldProblem;
if ~isempty(oldProblem)
    if ~strcmp(oldProblem, Global.problem.name)
        lastCurve = [];
        background = [];
        figure;
    end
end
if Global.terminated
    if isempty(lastCurve)
        return;
    end
    hold on;
    delete(lastCurve);
    background(numel(background) + 1).curve = semilogy(Global.trace(:, 1), Global.trace(:, 2) - Global.problem.idealfitness, 'LineWidth', 2);
    if isempty(Global.algorithm)
        algoname = "unknown";
    else
        algoname = string(Global.algorithm);
    end
    background(end).algoname = algoname;
    legend(cat(2, background.curve), mat2cell(cat(1, background.algoname), ones(1, numel(background)))');
    axis([0, 1, 0, inf]);
    lastCurve = [];
    set(gca, 'YScale', 'log');
    if draw_flag
        drawnow;
    end
    return;
end
if isempty(lastCurve)
    hold on;
    lastCurve = semilogy(Global.trace(:, 1), abs(Global.trace(:, 2) - Global.problem.idealfitness), 'LineWidth', 2);
else
    hold on;
    delete(lastCurve);
    lastCurve = semilogy(Global.trace(:, 1), abs(Global.trace(:, 2) - Global.problem.idealfitness), 'LineWidth', 2);
    hold off;
end
axis([0, 1, 0, inf]);
if isempty(Global.algorithm)
    algoname = "unknown";
else
    algoname = string(Global.algorithm);
end
currentLegend = sprintf("%s %.2f%% (%.2e)", algoname, 100 * Global.evaluated / Global.evaluation, Global.bestFitness);
if isempty(background)
    L = legend(lastCurve, currentLegend);
else
    L = legend([cat(2, background.curve), lastCurve], mat2cell([cat(1, background.algoname); currentLegend], ones(1, numel(background) + 1))');
end
set(L, 'Fontname', 'Book Antiqua', 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'YScale', 'log');
oldProblem = Global.problem.name;
if draw_flag
    drawnow;
end
end