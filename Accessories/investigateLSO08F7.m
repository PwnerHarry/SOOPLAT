function investigateLSO08F7(Global)
% Global = GLOBAL('-algorithm', 'BICCA2', '-problem', 'LSO08F7', '-evaluation', 1e99, '-dimension', 1000);
while true
    ACoDE('-initPop', Global.bestIndividual, '-dims', 1: 1000, '-NP', 500, '-maxFEs', 1e6, '-Global', Global);
end
end
