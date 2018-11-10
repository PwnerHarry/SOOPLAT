clear all;
mode = 'test';
algorithms = {'CCCMAES'};
suite = 'LSO08';
function_index = 1;
runtimes = 1;
dimension = 1000;
evaluation = dimension * 5e3;
DIR = 'Data';
mkdir(DIR);
Globals = {};
pointer = 1;
for R = 1: runtimes
    for i = function_index
        for j = 1: numel(algorithms)
            Globals{pointer} = GLOBAL('-algorithm', algorithms{j}, '-evaluation', evaluation, '-dimension', dimension, '-problem', sprintf('%sF%d', suite, i));
            pointer = pointer + 1;
        end
    end
    pause(1);
end
pointer = 1;
drawnow;
error_list = [];
while pointer <= numel(Globals)
    if strcmp(mode, 'benchmark')
        try
            feval(Globals{pointer}.algorithm, Globals{pointer});
        catch ME
            error_list = [error_list, pointer];
        end
    elseif strcmp(mode, 'test')
        feval(Globals{pointer}.algorithm, Globals{pointer});
    end
    G = Globals{pointer};
    renderCurve(G);
    filename = sprintf('[%s]%s_%s_D%d.mat', G.identifier, G.algorithm, G.problem.name, G.problem.dimension);
    save(fullfile(DIR, filename), 'G');
    pointer = pointer + 1;
end
error_Globals = Globals(error_list);