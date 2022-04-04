# SOOPLAT: PLATform for Single Objective Optimization

Increasing interest in the research of single objective optimization algorithms calls for an experimental platform that could provide researchers with convenient and universal interfaces of the benchmark functions and easy gathering statistics. Here, I propose a prototype Single Objective Optimization PLATform (SOOPLAT) to meet the needs.

SOOPLAT is inspired by PlatEMO [1] of Ye Tian and colleagues, which provided me great flexibility and convenience when I was composing the papers [2, 3]. Therefore, I have decided to make use of my time composing a convenient open-source platform for researchers for single objective optimization. The development of SOOPLAT began shortly after a chat with Daniel Molina, chair of the IEEE Task Force on Large-Scale Global Optimization (www.tflsgo.org), as both of us agreed that a platform is in great need.

I have to apologize that since I am currently busy with my studies and research, I cannot hereby provide a detailed users’ manual for this platform, which I may provide when I am done with my recent business. Since my code is clear to comprehend, please read the existing code of the algorithms I have implemented such as CSO [4], DECC-DG2 [5], etc., and you will certainly know how easy-to-use SOOPLAT is. Also, please contribute to the source code. You are the most welcome!

## References

[1] Y. Tian, et. al., “PlatEMO: A MATLAB platform for evolutionary multi-objective optimization [educational forum],” IEEE Computational Intelligence Magazine, vol. 12, no. 4, pp. 73–87, 2017.

[2] H.Ge and M. Zhao, et. al., “A many-objective evolutionary algorithm with two interacting processes: Cascade clustering and reference point incremental learning,” IEEE Transactions on Evolutionary Computation, Accepted and Online, 2018.

[3] M. Zhao and H. Ge, et. al., “A many-objective evolutionary algorithm with fast clustering and reference point redistribution,” in IEEE Congress on Evolutionary Computation, 2018, pp. 1–6.

[4] R. Cheng and Y. Jin, “A competitive swarm optimizer for large scale optimization,” IEEE Transactions on Cybernetics, vol. 45, no. 2, pp. 191–204, 2015.

[5] M. N. Omidvar, et. al., “DG2: A faster and more accurate differential grouping for large-scale black-box optimization,” IEEE Transactions on Evolutionary Computation, vol. 21, no. 6, pp. 929–942, 2017.

# REMINDER
1. run initialize.m to load the dependencies into MATLAB.
