function CSO(Global)
phi = 0.15;
if Global.problem.dimension >= 5000
    swarmSize = 1500;
elseif Global.problem.dimension >= 2000
    swarmSize = 1000;
elseif Global.problem.dimension >= 1000
    swarmSize = 500;
elseif Global.problem.dimension >= 500
    swarmSize = 250;
else
    swarmSize = 100;
end
halfSize = ceil(swarmSize / 2);
p = initPopulation([], swarmSize, Global);
fitness = Global.evaluate(p);
v = zeros(swarmSize, Global.problem.dimension);
gen = 0;
while ~Global.terminated
    gen = gen + 1;
    fprintf('CSO: CYCLE %d\n', gen);
    % generate random pairs and do pairwise competitions
    rlist = randperm(swarmSize);
    rpairs = [rlist(1: ceil(swarmSize / 2)); rlist(floor(swarmSize / 2) + 1:swarmSize)]';
    mask = (fitness(rpairs(:, 1)) > fitness(rpairs(:, 2)));
    losers = mask .* rpairs(:, 1) + ~mask .* rpairs(:, 2);
    winners = ~mask .* rpairs(:, 1) + mask .* rpairs(:, 2);
    % losers learn from winners
    randco1 = rand(halfSize, Global.problem.dimension);
    randco2 = rand(halfSize, Global.problem.dimension);
    randco3 = rand(halfSize, Global.problem.dimension);
    v(losers,:) = randco1.*v(losers,:) + randco2.*(p(winners,:) - p(losers,:)) + phi * randco3 .* (ones(halfSize, 1) * mean(p) - p(losers,:));
    p(losers,:) = p(losers,:) + v(losers,:);
    for i = 1: halfSize
        p(losers(i),:) = max(p(losers(i),:), Global.problem.lowerbound);
        p(losers(i),:) = min(p(losers(i),:), Global.problem.upperbound);
    end
    fitness(losers, :) = Global.evaluate(p(losers, :));
    Global.draw();
end
end

