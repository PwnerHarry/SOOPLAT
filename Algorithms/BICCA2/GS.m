function cratio = GS(Global)
GNP = 50;
FEs = 25 * Global.problem.dimension;
fitness_before_GS = Global.bestFitness;
fprintf('GS: ');
LSHADE('-initPop', Global.bestIndividual, '-initFit', Global.bestFitness, '-contextFitness', Global.bestFitness, '-contextVector', Global.bestIndividual, '-initPop', Global.bestIndividual, '-maxFEs', FEs, '-NP', GNP, '-dims', 1: Global.problem.dimension, '-Global', Global);
fitness_after_GS = Global.bestFitness;
cratio = abs((fitness_before_GS - fitness_after_GS) / fitness_before_GS);
hold on;
fill([0, 1, 1, 0], [fitness_after_GS, fitness_after_GS, fitness_before_GS, fitness_before_GS], [1, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', '0.2');
drawnow;
hold off;
fprintf('CRATIO = %.2f%%\n', 100 * cratio);
end

