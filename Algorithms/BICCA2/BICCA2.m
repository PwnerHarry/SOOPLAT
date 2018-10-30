function BICCA2(Global)
global LS_RESTART_FLAG CC_RESTART_FLAG BICCA2_LS_STAGNANT_FLAG;
recordGroupings('-Global', Global);
LS_RESTART_FLAG = true;
CC_RESTART_FLAG = true;
BICCA2_LS_STAGNANT_FLAG = false;
adapter = BOLTZADA({'LS', 'CC', 'GS'}, 5);
cycles = 0;
NP = 50;
[~, SOA] = sort(rand(50, Global.problem.dimension));
X = SOA ./ NP .* repmat(Global.problem.upperbound-Global.problem.lowerbound, NP, 1) + repmat(Global.problem.lowerbound, NP, 1);
Global.evaluate(X);
patterns = [];
while ~Global.terminated
    startFEs = Global.evaluated;
    choice = adapter.decide();
    if strcmp(choice, 'LS')
        cratio = LS(Global);
    elseif strcmp(choice, 'CC')
        [cratio, patterns] = CC(Global, patterns);
    elseif strcmp(choice, 'GS')
        cratio = GS(Global);
    end
    adapter.update(choice, cratio / (Global.evaluated - startFEs));
end
fprintf('cycles: %d\n', cycles);