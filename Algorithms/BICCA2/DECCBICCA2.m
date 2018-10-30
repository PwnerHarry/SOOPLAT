function DECCBICCA2(Global)
global CC_RESTART_FLAG CC_OPTIMIZER;
recordGroupings('-Global', Global);
CC_RESTART_FLAG = true;
CC_OPTIMIZER = 'SaNSDE';
patterns = [];
NP = 15;
pop = initPopulation([], NP, Global);
val = Global.evaluate(pop);
while ~Global.terminated
    [cratio, patterns] = CC(Global, patterns);
end
end