/*
 * =====================================================================================
 *
 *       Filename:  CCVIL.cpp
 *
 * Description:  The main source file of CCVIL algorithm's implementation, it includes
 *    				two stages:
 *    					1) Learning Stage
 *    					2) Optimization Stage
 *
 *        Version:  1.0
 *        Created:  02/24/2011 07:56:20 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wenxiang Chen (http://cs-chen.net), chenwx.ustc@gmail.com
 *        Company:  Nature Inspired Computation and Application Laboratory (NICAL), USTC
 *
 * =====================================================================================
 */

#include "CCVIL.h"

CCVIL::CCVIL(RunParameter* runParam){
	param = runParam;
	p = randPerm(param->dimension);

	Rng::seed(param->initRandomSeed);

	lookUpGroup = new unsigned[param->dimension];
	impreciseGroup = new bool[param->dimension];

	// MaxFitEval = runParam->fitnessCheckPoint[(*runParam).fitnessCheckPoint.size()-1];
	MaxFitEval = param->dimension * 5000;
	cout<<"Max Fitness Evaluation = "<<MaxFitEval<<endl;

	if (param->learnStrategy==0){
		lowerThreshold = runParam->lowerThreshold;
		upperThreshold = min(round(MaxFitEval*param->learnPortion/(runParam->dimension*((1+1)*(3)+1))), 800.0); 
	}
	else if (param->learnStrategy==5){
		lowerThreshold = 4 * param->dimension* (1+log(param->dimension)/log(2)); 
		upperThreshold = 4 * (log(1-0.95)/log(1 - 2/(double)param->dimension)) * (1+log(param->dimension)/log(2)); 
	}
	else if (param->learnStrategy>=4){
		//		lowerThreshold = 3*round(log(1-0.995)/log(1-1/(double)param->dimension));
		//		upperThreshold = 3*(param->dimension)*((param->dimension)-1)/2; 
		lowerThreshold = 3000000; 
		upperThreshold = 3000000; 
	}
	cout<<"Lower threshold = "<<lowerThreshold<<", Upper threshold = "<<upperThreshold<<endl;
	bestCand = new IndividualT<double>(ChromosomeT<double>(param->dimension));
}

CCVIL::~CCVIL(){
	delete bestCand;
	//	delete[] p;
	delete[] impreciseGroup;
	delete[] lookUpGroup;
}

void CCVIL::run(){
	struct timeval start, end;
	long seconds, useconds;    
	double mtime;
	//	printf ( "lookUpGroup\n" );
	//	printArray(lookUpGroup, runParam->dimension);

	// check if folders "result" and "trace" exist or not
	mkdir ("result", O_CREAT|S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
	mkdir ("trace", O_CREAT|S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	if (param->learnStrategy >= 1 && param->learnStrategy <= 3){
		getPriorInterStage(); 
	}

	// in result folder
	string resultStr("result/resF");
	resultStr += itos(fp->getID());
	if (param->learnStrategy >= 1 && param->learnStrategy <= 3){
		resultStr += "-P";
		resultStr += itos (floor(param->knownGroupPercent[0]*100));
	}
	resultStr += "-S";
	resultStr += itos(param->learnStrategy);
	resultStr += "-D";
	resultStr += itos(param->dimension);
	resultStr += ".txt";
	printf("resultStr = %s", resultStr.c_str());
	printf("\n");
	resultFP = fopen(resultStr.c_str(), "w");

	string timeStr("result/timeF");
	timeStr += itos(fp->getID());
	if (param->learnStrategy >= 1 && param->learnStrategy <= 3){
		timeStr += "-P";
		timeStr += itos (floor(param->knownGroupPercent[0]*100));
	}
	timeStr += "-S";
	timeStr += itos(param->learnStrategy);
	timeStr += "-D";
	timeStr += itos(param->dimension);
	timeStr += ".txt";
	printf("timeStr = %s", timeStr.c_str());
	printf("\n");
	timeFP = fopen(timeStr.c_str(), "w");

	for (unsigned i=0; i < param->numOfRun; i++){
		printf ( "\n\n\n========================== F %d, Run %d ========================\n\n\n", fp->getID(), i+1 );

		//		if ( param->learnStrategy >= 1 && param->learnStrategy<=3 ){
		/************************* re-initialize the sampling points *************************/
		groupInfo.clear();
		// initialize the groupInfo
		for (unsigned j = 0; j<param->dimension; j++){
			vector<unsigned> tempVec;
			tempVec.push_back(j);
			lookUpGroup[j] = j;
			groupInfo.push_back(tempVec);
		}
		//		}

		samplingPoints.clear();
		for (unsigned j=0 ; j<= param->samplingPoint; j++){
			samplingPoints.push_back(j*param->samplingInterval);
		}

		/************************* in trace folder *************************/
		// store grouping information
		string groupStr("trace/groupF");
		groupStr += itos(fp->getID());
		groupStr += "-R";
		groupStr += itos(i+1);
		if (param->learnStrategy >= 1 && param->learnStrategy <= 3){
			groupStr += "-P";
			groupStr += itos (floor(param->knownGroupPercent[0]*100));
		}
		groupStr += "-S";
		groupStr += itos(param->learnStrategy);
		groupStr += "-D";
		groupStr += itos(param->dimension);
		groupStr += ".txt";
		printf("groupStr = %s", groupStr.c_str());
		printf("\n");
		groupFP = fopen(groupStr.c_str(), "w");

		string groupFesStr("trace/groupFesF");
		groupFesStr += itos(fp->getID());
		groupFesStr += "-R";
		groupFesStr += itos(i+1);
		if (param->learnStrategy >= 1 && param->learnStrategy <= 3){
			groupFesStr += "-P";
			groupFesStr += itos (floor(param->knownGroupPercent[0]*100));
		}
		groupFesStr += "-S";
		groupFesStr += itos(param->learnStrategy);
		groupFesStr += "-D";
		groupFesStr += itos(param->dimension);
		groupFesStr += ".txt";
		printf("groupFesStr = %s", groupFesStr.c_str());
		printf("\n");
		groupFesFP = fopen(groupFesStr.c_str(), "w");

		string fesStr("trace/fesF");
		fesStr += itos(fp->getID());
		fesStr += "-R";
		fesStr += itos(i+1);
		if (param->learnStrategy >= 1 && param->learnStrategy <= 3){
			fesStr += "-P";
			fesStr += itos (floor(param->knownGroupPercent[0]*100));
		}
		fesStr += "-S";
		fesStr += itos(param->learnStrategy);
		fesStr += "-D";
		fesStr += itos(param->dimension);
		fesStr += ".txt";
		printf("fesStr = %s", fesStr.c_str());
		printf("\n");
		fesFP = fopen(fesStr.c_str(), "w");

		string valStr("trace/valF");
		valStr += itos(fp->getID());
		valStr += "-R";
		valStr += itos(i+1);
		if (param->learnStrategy >= 1 && param->learnStrategy <= 3){
			valStr += "-P";
			valStr += itos (floor(param->knownGroupPercent[0]*100));
		}
		valStr += "-S";
		valStr += itos(param->learnStrategy);
		valStr += "-D";
		valStr += itos(param->dimension);
		valStr += ".txt";
		printf("valStr = %s", valStr.c_str());
		printf("\n");
		valFP = fopen(valStr.c_str(), "w");

		fes = 0;
		bestFit = DBL_MAX;
		/* algorithm runing part: start */
		gettimeofday(&start, NULL);

		groupRec.push_back(groupInfo.size());
		groupFesRec.push_back(fes);

		//		printf ( "Vector, size of interaction array =%d\n", (fp->getInterArray()).size() );
		//		printVector(fp->getInterArray());

		if (param->learnStrategy == 0){
			learningStage();
		}else if (param->learnStrategy==4 || param->learnStrategy==6){
			sampleLearnStage();
		}else if (param->learnStrategy==5){
			binSearchLearnStage(); 
		}else if (param->learnStrategy==7){
			RandomSampleGenDef(); 
		}else if (param->learnStrategy==8){
			RandomWalkGenDef(); 
		}else if (param->learnStrategy==9){
			BinSearchRandWalk(); 
		}

		if (param->learnStrategy >= 1 && param->learnStrategy <= 3){
			// initialize bestCand
			(*bestCand)[0].initialize(fp->getMinX(), fp->getMaxX());
		}

		if (param->performOpt == 1){
			optimizationStage();
		}

		gettimeofday(&end, NULL);
		/* algorithm runing part: end */

		seconds  = end.tv_sec  - start.tv_sec;
		useconds = end.tv_usec - start.tv_usec;

		mtime = (((seconds) * 1000 + useconds/1000.0) + 0.5)/1000;
		printf ( "Result = %.8e, Running Time = %fs\n", bestFit, mtime);

		resultRec.push_back(bestFit);
		timeRec.push_back(mtime);

		printf ( "\n\n\n========================================================\n\n\n" );

		for (unsigned i=0; i<groupRec.size(); i++){
			fprintf(groupFP, "%d\n", groupRec[i]);
		}
		fclose(groupFP);
		groupRec.clear();

		for (unsigned i=0; i<groupFesRec.size(); i++){
			fprintf(groupFesFP, "%d\n", groupFesRec[i]);
		}
		fclose(groupFesFP);
		groupFesRec.clear();

		for (unsigned i=0; i<fesRec.size(); i++){
			fprintf(fesFP, "%d\n", fesRec[i]);
		}
		fclose(fesFP);
		fesRec.clear();

		for (unsigned i=0; i<valRec.size(); i++){
			fprintf(valFP, "%.8e\n", valRec[i]);
		}
		fclose(valFP);
		valRec.clear();
	}

	// delete all file pointers

	// result
	for (unsigned i=0; i<resultRec.size(); i++){
		fprintf(resultFP, "%.8e\n", resultRec[i]);
	}
	fclose(resultFP);
	resultRec.clear();

	// time
	for (unsigned i=0; i<timeRec.size(); i++){
		fprintf(timeFP, "%.8e\n", timeRec[i]);
	}
	fclose(timeFP);
	timeRec.clear();
}

/* 
 * procedure of learning stage, update the "groupInfo"
 */
void CCVIL::learningStage(){

	cout<<"Learning Stage ... "<<endl;
	cycle = 1; 
	int lastCycleIndex = -1;
	bool learnStageFlag, needCapture, isSameGroup = false, separableFunc = true; // assume every benchmark function is separable at the first beginning

	(*bestCand)[0].initialize(fp->getMinX(), fp->getMaxX());

	//	printf ( "Best Cand\n" );
	//	printPopulation((*bestCand));

	popGenerate(true);

	//	printf ( "Compute Learn Stage Flag = %d\n", !(cycle > upperThreshold || (cycle > lowerThreshold && groupInfo.size() == param->dimension) || groupInfo.size() == 1));

	while ( (learnStageFlag = !(cycle > upperThreshold 
					|| (cycle > lowerThreshold && groupInfo.size() == param->dimension) 
					|| groupInfo.size() == 1))==true ){
		// start a new cycle
		//		printf("===================================================\n=================== New Cycle %d ===================\n===================================================\n", cycle);

		for (unsigned i=0; i<param->dimension; i++) {
			//			printf("\n=========================================================================\nPhase = %d, Cycle = %d, impresice in last cycle: %d\n",i, cycle, impreciseGroup[lastCycleIndex]);
			if (i == 0){
				// start a new phase for each cycle, each phase checking one dimension, i.e., p[i]
				//	printf("Population Re-initialization & Issue Random Permutation\n");
				popInit();

				//				printf ( "beforeInitBestCand\n" );
				//				print2Dvector(pop);
				//				printPopulation((*bestCand));

				//				if (cycle==1){
				//					initBestCand(learnStageFlag);
				//				}

				//				printf ( "afterInitBestCand\n" );
				//				printPopulation((*bestCand));

				p = randPerm(param->dimension);

				//				printf("Random Permutation p:\n");
				//				printArray(p, param->dimension);

				lastCycleIndex = -1;

				// re-initialize the vector to all false, when starting a new cycle
				for (unsigned j=0; j<param->dimension; j++){
					impreciseGroup[j] = false;
				}
			}


			needCapture = groupInfo.size()!=1 && ((cycle<= lowerThreshold) ||(separableFunc == false && cycle <= upperThreshold)) && lastCycleIndex!=-1;
			//			printf("Need Capture ? %d: between %d and %d\n", needCapture, p[i], lastCycleIndex);

			if (lastCycleIndex!=-1){
				//	to decide whether current dimesion are in the same group with last dimension
				isSameGroup = sameGroup(p[i], lastCycleIndex);
				//				printf("In the Same group? %d\n", isSameGroup);
			}

			if (lastCycleIndex == -1 || isSameGroup == false){
				// if current dimension and last optimized index are in the different group, then optimize on current dimension
				//				printf ( "Index = %d\n", p[i] );
				JADECC(p[i], true); 	/*learnStage = true*/
			} else {
				needCapture = false;
			}

			// begin interaction capture process
			if ( needCapture == true && impreciseGroup[lastCycleIndex] == false ){
				captureInter( p[i],lastCycleIndex );
				separableFunc = false;
			}

			// record last optimized dimension
			if ( lastCycleIndex==-1 || isSameGroup == false ){
				lastCycleIndex = p[i];
			}

			if (i == param->dimension-1){
				delete[] p;
			}
		}
		printf("F %d, Learning Cycle =%d, fes = %ld, GroupAmount = %d, BestVal = %.8e \n",
				fp->getID(), 	cycle, 	fes,	(int)groupInfo.size(), bestCand->fitnessValue());

		cycle++;
	}

	//	printf ( "Before sorting, Group info\n" );
	//	print2Dvector(groupInfo);
	//
	//	printf ( "Look up group table\n" );
	//	printArray(lookUpGroup, param->dimension);

	sortGroupInfo();

	//	printf ( "After sorting, Group info\n" );
	//	print2Dvector(groupInfo);
	//
	//	printf ( "Look up group table\n" );
	//	printArray(lookUpGroup, param->dimension);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CCVIL::binSearchLearnStage
 *  Description:  locate the interaction with binary search approach, which intends to 
 *  							reduce the runtime complexity from O(l*l) to O( l * lg l)
 * =====================================================================================
 */
	void
CCVIL::binSearchLearnStage ()
{
	int interIndex; 
	// intialize groupInfo to {{1}, {2}, ..., {N}}
	groupInfo.clear();
	// initialize the groupInfo
	for (unsigned j = 0; j<param->dimension; j++){
		vector<unsigned> tempVec;
		tempVec.push_back(j);
		lookUpGroup[j] = j;
		groupInfo.push_back(tempVec);
	}

	// lowerThreshold: stop criteria if the problem is found separable, i.e., groupInfo.size() == dimension.
	// upperThreshold: stop criteria if the problem is found non-separable, i.e., groupInfo.size() < dimension.
	while (fes < upperThreshold && (fes<=lowerThreshold || groupInfo.size()<param->dimension)) {
		unsigned indexI = floor(Rng::uni()*param->dimension);
		//		printf ( "index I = %d\n", indexI );

		IndividualT<double> indiv1(ChromosomeT<double>(param->dimension));
		indiv1[0].initialize(fp->getMinX(), fp->getMaxX()); 
		//		printf ( "indiv 1\n" );
		//		printPopulation(indiv1); 

		IndividualT<double> indiv0(ChromosomeT<double>(param->dimension));
		indiv0[0].initialize(fp->getMinX(), fp->getMaxX()); 
		//		printf ( "indiv 0\n" );
		//		printPopulation(indiv0); 

		unsigned groupI = lookUpGroup[indexI]; 
		//		printf ( "group Index that indexI belongs to: %d\n", groupI );

		IndividualT<double> indiv2(indiv0);
		for (unsigned i=0; i < groupInfo[groupI].size(); i++){
			indiv2[0][groupInfo[groupI][i]] = indiv1[0][groupInfo[groupI][i]]; 
		}
		//		printf ( "indiv2\n" );
		//		printPopulation(indiv2); 

		if ( testInteraction(indiv1, indiv2, indexI) ){
			//			printf ( "fes delta not equal\n");
			interIndex = findInteractPosition(indiv1, indiv2, indexI); 
		}else {
			//			printf ( "fes delta equal\n" );
			interIndex = -1; 
		}

		if (interIndex != -1){
			unsigned group1 = lookUpGroup[indexI];
			unsigned group2 = lookUpGroup[interIndex];
			combineGroup( group1, group2 );
			printf ("%d\t&\t%d:\t%ld\t%d\n", indexI, interIndex, fes, groupInfo.size());
		}

		if (groupInfo.size()==1){
			printf ( "Converge to one single group\n" );
			break; 
		}
	}

	(*bestCand)[0].initialize(fp->getMinX(), fp->getMaxX());
	sortGroupInfo();
}		/* -----  end of function CCVIL::binSearchLearnStage  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CCVIL::BinSearchRandWalk
 *  Description:  Random Walk with Binary Search which helps to accelerate the speed of 
 *  							learning
 * =====================================================================================
 */
	void
CCVIL::BinSearchRandWalk ()
{
	int interIndex; 
	// intialize groupInfo to {{1}, {2}, ..., {N}}
	groupInfo.clear();
	// initialize the groupInfo
	for (unsigned j = 0; j<param->dimension; j++){
		vector<unsigned> tempVec;
		tempVec.push_back(j);
		lookUpGroup[j] = j;
		groupInfo.push_back(tempVec);
	}

	IndividualT<double> localBest(ChromosomeT<double>(param->dimension));
	localBest[0].initialize(fp->getMinX(), fp->getMaxX()); 
	// lowerThreshold: stop criteria if the problem is found separable, i.e., groupInfo.size() == dimension.
	// upperThreshold: stop criteria if the problem is found non-separable, i.e., groupInfo.size() < dimension.
	while (fes < upperThreshold && (fes <= lowerThreshold || groupInfo.size()<param->dimension)) {
		vector<	IndividualT<double> > improveIndivVec; 

		unsigned indexI = floor(Rng::uni()*param->dimension);
		IndividualT<double> indiv0(ChromosomeT<double>(param->dimension));
		indiv0[0].initialize(fp->getMinX(), fp->getMaxX()); 

		unsigned groupI = lookUpGroup[indexI]; 
		//		printf ( "group Index that indexI belongs to: %d\n", groupI );

		for ( unsigned i=0; i < groupInfo[groupI].size(); i++ ){
			indiv0[0][groupInfo[groupI][i]] = localBest[0][groupInfo[groupI][i]]; 
		}

		//		printf ( "indiv2\n" );
		//		printPopulation(indiv2); 

		if ( TestInterWalk(indiv0, indexI, localBest, improveIndivVec) ){
			//			printf ( "fes delta not equal\n");
			interIndex = findInterPosWalk(localBest, indiv0, indexI, improveIndivVec); 
		} else {
			//			printf ( "fes delta equal\n" );
			interIndex = -1; 
		}

		if (interIndex != -1){
			unsigned group1 = lookUpGroup[indexI];
			unsigned group2 = lookUpGroup[interIndex];
			combineGroup( group1, group2 );
			printf ("%d\t&\t%d:\t%ld\t%d\n", indexI, interIndex, fes, groupInfo.size());
		}

		if (groupInfo.size()==1) {
			printf ( "Converge to one single group\n" );
			break; 
		}

		// update the localBest by selecting the best Individual from the improveIndivVec
		//		if (improveIndivVec.size()>0){
		//			printf ( "Number of improved individual vector = %d\n",  improveIndivVec.size());
		//			unsigned bestIndivIndex=0; 
		//			for (unsigned i=0; i< improveIndivVec.size(); i++){
		//				if (improveIndivVec[i].getFitness()<localBest.getFitness()&&improveIndivVec[i].getFitness()<improveIndivVec[bestIndivIndex].getFitness()){
		//					bestIndivIndex = i; 
		//					printf ( "Update i = %d\n" , i);
		//				}	
		//			}
		//			localBest = improveIndivVec[bestIndivIndex]; 
		//		}
	}

	(*bestCand)[0].initialize(fp->getMinX(), fp->getMaxX());
	sortGroupInfo();
}		/* -----  end of function CCVIL::BinSearchRandWalk  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CCVIL::findInteractPosition
 *  Description:  find the interaction position using binary search approach
 *  							Translation from 's' in the original paper to indiv-like variable:
 *
 *  							indiv1 -> s
 *  							indiv2 -> s'
 *  							indiv3 -> s_2
 * =====================================================================================
 */
	int	
CCVIL::findInteractPosition ( IndividualT<double> indiv1, IndividualT<double> indiv2, unsigned indexI )
{
	vector<unsigned> diffDim; 

	//	printf("===================================================\n");
	//	printf ( "Find Interacting Position\n" );

	//	printf ( "indiv 1\n" );
	//	printPopulation(indiv1); 

	//	printf ( "indiv 2\n" );
	//	printPopulation(indiv2); 

	// find the dimensions where indiv1 & indiv2 differ
	for (unsigned i=0; i<param->dimension; i++){
		if (indiv1[0][i] != indiv2[0][i]){
			diffDim.push_back(i);
		}
	}
	//	printf ( "diff Dim vector\n" );
	//	printVector(diffDim); 

	if (diffDim.size()==1){
		return diffDim[0]; 
	} 

	IndividualT<double> indiv3(indiv2); 
	for (unsigned i=0; i<diffDim.size()/2; i++){ 
		indiv3[0][diffDim[i]] = indiv1[0][diffDim[i]]; 
	} //	printf ( "diffDim.size/2 = %d\n", diffDim.size()/2 );

	//	printf ( "indiv 3\n" );
	//	printPopulation(indiv3); 

	if ( testInteraction(indiv1, indiv3, indexI)){
		return findInteractPosition(indiv1, indiv3, indexI); 
	}else{
		return findInteractPosition(indiv2, indiv3, indexI); 
	}

	//	printf("===================================================\n");
}		/* -----  end of function CCVIL::findInteractPosition  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CCVIL::findInterPosWalk
 *  Description:  Random Walk version of find Interaction position position
 * =====================================================================================
 */
	int
CCVIL::findInterPosWalk ( IndividualT<double> &localBest, IndividualT<double> indiv2, unsigned indexI, vector< IndividualT<double> > &improveIndivVec )
{
	vector<unsigned> diffDim; 

	// find the dimensions where indiv1 & indiv2 differ
	for (unsigned i=0; i<param->dimension; i++){
		if (localBest[0][i] != indiv2[0][i]){
			diffDim.push_back(i);
		}
	}

	if (diffDim.size()==1){
		return diffDim[0]; 
	} else if (diffDim.size() == 0){
		printf ( "Size of diffDim = 0\n" );
		exit(EXIT_FAILURE); 
	}

	IndividualT<double> indiv3(indiv2); 
	for (unsigned i=0; i<diffDim.size()/2; i++){ 
		indiv3[0][diffDim[i]] = localBest[0][diffDim[i]]; 
	} 

	if ( TestInterWalk( indiv3, indexI, localBest,improveIndivVec  )){
		return findInterPosWalk(localBest, indiv3, indexI,improveIndivVec ); 
	}else{
		return findInterPosWalk(indiv2, indiv3, indexI,improveIndivVec ); 
	}
}		/* -----  end of function CCVIL::findInterPosWalk  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CCVIL::testInteraction
 *  Description:  
 * =====================================================================================
 */
	bool
CCVIL::testInteraction (  IndividualT<double> indiv1, IndividualT<double> indiv2, unsigned indexI  )
{
	double randI = Rng::uni() * (fp->getMaxX() - fp->getMinX()) + fp->getMinX(); 
	//	printf ( "mutated gene = %.16f\n", randI );

	IndividualT<double> indiv1_sub(indiv1);
	indiv1_sub[0][indexI] = randI; 
	//	printf ( "indiv 1 sub\n" );
	//	printPopulation(indiv1_sub);

	IndividualT<double> indiv2_sub(indiv2);
	indiv2_sub[0][indexI] = randI;
	//	printf ( "indiv 2 sub\n" );
	//	printPopulation(indiv2_sub);

	indiv1.setFitness( fp->compute(indiv1[0]) );
	fes += 1; 
	sampleInfo(indiv1.getFitness());
	if (indiv1.getFitness()<bestFit){
		bestFit = indiv1.getFitness(); 
	}

	indiv2.setFitness( fp->compute(indiv2[0]) );
	fes += 1; 
	sampleInfo(indiv2.getFitness());
	if (indiv2.getFitness()<bestFit){
		bestFit = indiv2.getFitness(); 
	}

	indiv1_sub.setFitness( fp->compute(indiv1_sub[0]) );
	fes += 1; 
	sampleInfo(indiv1_sub.getFitness());
	if (indiv1_sub.getFitness()<bestFit){
		bestFit = indiv1_sub.getFitness(); 
	}

	indiv2_sub.setFitness( fp->compute(indiv2_sub[0]) );
	fes += 1; 
	sampleInfo(indiv2_sub.getFitness());
	if (indiv2_sub.getFitness()<bestFit){
		bestFit = indiv2_sub.getFitness(); 
	}

	double fesIndiv1Delta = indiv1.getFitness() - indiv1_sub.getFitness(); 
	double fesIndiv2Delta = indiv2.getFitness() - indiv2_sub.getFitness(); 
	//	printf ( "fes: indiv1 =\t\t%.16f\tindiv2 =\t%.16f\nfes: indiv1_sub =\t%.16f\tindiv2_sub =\t%.16f\n",indiv1.getFitness(),indiv2.getFitness(),indiv1_sub.getFitness(),indiv2_sub.getFitness());
	//	printf ( "fes delta: indiv1 = %.30f, indiv2 = %.30f\n", fesIndiv1Delta, fesIndiv2Delta);
	//
	if ( abs(fesIndiv1Delta - fesIndiv2Delta)/max(abs(indiv1.getFitness()), abs(indiv2.getFitness())) >1e-10) {
		//		printf ( "**************** Index %d ******************************\n", indexI );
		//		printf ( "fes: indiv1 =\t\t%.16f\tindiv2 =\t%.16f\nfes: indiv1_sub =\t%.16f\tindiv2_sub =\t%.16f\n",indiv1.getFitness(),indiv2.getFitness(),indiv1_sub.getFitness(),indiv2_sub.getFitness());
		//		printf ( "fes delta: indiv1 = %.30f, indiv2 = %.30f\n", fesIndiv1Delta, fesIndiv2Delta);
		//		printf ( "the difference = %.30f, differ rate = %.10f, differ rate = %.10f\n", fesIndiv1Delta - fesIndiv2Delta, abs((fesIndiv1Delta - fesIndiv2Delta)/max(abs(indiv1.getFitness()), abs(indiv2.getFitness()))), abs(fesIndiv1Delta - fesIndiv2Delta)/max(abs(indiv1.getFitness()), abs(indiv2.getFitness())) );
		//		printf ( "*******************************************************\n\n" );
		return true; 
	}else{
		return false; 
	}
}		/* -----  end of function CCVILL::testInteraction  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CCVIL::TestInterWalk
 *  Description:  add a new parameter, localBest. It is updated inside the function
 * =====================================================================================
 */
	bool
CCVIL::TestInterWalk ( IndividualT<double> indiv2, unsigned indexI, IndividualT<double> &localBest )
{
	double randI = Rng::uni() * (fp->getMaxX() - fp->getMinX()) + fp->getMinX(); 

	IndividualT<double> indiv1(localBest); 

	IndividualT<double> indiv1_sub(indiv1);
	indiv1_sub[0][indexI] = randI; 

	IndividualT<double> indiv2_sub(indiv2);
	indiv2_sub[0][indexI] = randI;

	indiv1.setFitness( fp->compute(indiv1[0]) );
	fes += 1; 
	sampleInfo(indiv1.getFitness());
	if (indiv1.getFitness()<bestFit){
		bestFit = indiv1.getFitness(); 
	}

	indiv2.setFitness( fp->compute(indiv2[0]) );
	fes += 1; 
	sampleInfo(indiv2.getFitness());
	if (indiv2.getFitness()<bestFit){
		bestFit = indiv2.getFitness(); 
	}

	indiv1_sub.setFitness( fp->compute(indiv1_sub[0]) );
	fes += 1; 
	sampleInfo(indiv1_sub.getFitness());
	if (indiv1_sub.getFitness()<bestFit) {
		bestFit = indiv1_sub.getFitness(); 
	}

	indiv2_sub.setFitness( fp->compute(indiv2_sub[0]) );
	fes += 1; 
	sampleInfo(indiv2_sub.getFitness());
	if (indiv2_sub.getFitness()<bestFit){
		bestFit = indiv2_sub.getFitness(); 
	}

	// update the localBest for producing random candidate afterwards
	// TODO See whether update focuses on one dimension will help
	if (indiv1.getFitness() > indiv2.getFitness()){
		localBest = indiv2 ; 
	}

	if (indiv1.getFitness() > indiv1_sub.getFitness()){
		// indiv1_sub is better
		localBest[0][indexI] = indiv1_sub[0][indexI]; 
	}

	double fesIndiv1Delta = indiv1.getFitness() - indiv1_sub.getFitness(); 
	double fesIndiv2Delta = indiv2.getFitness() - indiv2_sub.getFitness(); 
	//	printf ( "fes: indiv1 =\t\t%.16f\tindiv2 =\t%.16f\nfes: indiv1_sub =\t%.16f\tindiv2_sub =\t%.16f\n",indiv1.getFitness(),indiv2.getFitness(),indiv1_sub.getFitness(),indiv2_sub.getFitness());
	//	printf ( "fes delta: indiv1 = %.30f, indiv2 = %.30f\n", fesIndiv1Delta, fesIndiv2Delta);
	//

	//	for (unsigned i=0; i<param->dimension; i++){
	//		if (indiv1[0][i] != localBest[0][i]){
	//			printf ( "Local Best Updatedm\n" );
	//			break; 
	//		}
	//	}		

	if ( abs(fesIndiv1Delta - fesIndiv2Delta)/max(abs(indiv1.getFitness()), abs(indiv2.getFitness())) >1e-10) {
		//		printf ( "**************** Index %d ******************************\n", indexI );
		//		printf ( "fes: indiv1 =\t\t%.16f\tindiv2 =\t%.16f\nfes: indiv1_sub =\t%.16f\tindiv2_sub =\t%.16f\n",indiv1.getFitness(),indiv2.getFitness(),indiv1_sub.getFitness(),indiv2_sub.getFitness());
		//		printf ( "fes delta: indiv1 = %.30f, indiv2 = %.30f\n", fesIndiv1Delta, fesIndiv2Delta);
		//		printf ( "the difference = %.30f, differ rate = %.10f, differ rate = %.10f\n", fesIndiv1Delta - fesIndiv2Delta, abs((fesIndiv1Delta - fesIndiv2Delta)/max(abs(indiv1.getFitness()), abs(indiv2.getFitness()))), abs(fesIndiv1Delta - fesIndiv2Delta)/max(abs(indiv1.getFitness()), abs(indiv2.getFitness())) );
		//		printf ( "*******************************************************\n\n" );
		return true; 
	}else{
		return false; 
	}
}		/* -----  end of function CCVIL::TestInterWalk  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CCVIL::TestInterWalk
 *  Description:  add a new parameter, localBest. It is updated inside the function
 * =====================================================================================
 */
	bool
CCVIL::TestInterWalk ( IndividualT<double> indiv2, unsigned indexI, IndividualT<double> &localBest, vector< IndividualT<double> > &improveIndivVec)
{
	double randI = Rng::uni() * (fp->getMaxX() - fp->getMinX()) + fp->getMinX(); 

	IndividualT<double> indiv1(localBest); 

	IndividualT<double> indiv1_sub(indiv1);
	indiv1_sub[0][indexI] = randI; 

	IndividualT<double> indiv2_sub(indiv2);
	indiv2_sub[0][indexI] = randI;

	indiv1.setFitness( fp->compute(indiv1[0]) );
	fes += 1; 
	sampleInfo(indiv1.getFitness());
	if (indiv1.getFitness()<bestFit){
		bestFit = indiv1.getFitness(); 
	}

	indiv2.setFitness( fp->compute(indiv2[0]) );
	fes += 1; 
	sampleInfo(indiv2.getFitness());
	if (indiv2.getFitness()<bestFit){
		bestFit = indiv2.getFitness(); 
	}

	indiv1_sub.setFitness( fp->compute(indiv1_sub[0]) );
	fes += 1; 
	sampleInfo(indiv1_sub.getFitness());
	if (indiv1_sub.getFitness()<bestFit) {
		bestFit = indiv1_sub.getFitness(); 
	}

	indiv2_sub.setFitness( fp->compute(indiv2_sub[0]) );
	fes += 1; 
	sampleInfo(indiv2_sub.getFitness());
	if (indiv2_sub.getFitness()<bestFit){
		bestFit = indiv2_sub.getFitness(); 
	}

	// update the localBest for producing random candidate afterwards
	// TODO See whether update focuses on one dimension will help
	if (indiv1.getFitness() > indiv2.getFitness()){
		improveIndivVec.push_back(indiv2); 
	}

	if (indiv1.getFitness() > indiv1_sub.getFitness()){
		// indiv1_sub is better
		improveIndivVec.push_back(indiv1_sub); 
	}

	double fesIndiv1Delta = indiv1.getFitness() - indiv1_sub.getFitness(); 
	double fesIndiv2Delta = indiv2.getFitness() - indiv2_sub.getFitness(); 
	//	printf ( "fes: indiv1 =\t\t%.16f\tindiv2 =\t%.16f\nfes: indiv1_sub =\t%.16f\tindiv2_sub =\t%.16f\n",indiv1.getFitness(),indiv2.getFitness(),indiv1_sub.getFitness(),indiv2_sub.getFitness());
	//	printf ( "fes delta: indiv1 = %.30f, indiv2 = %.30f\n", fesIndiv1Delta, fesIndiv2Delta);
	//

	//	for (unsigned i=0; i<param->dimension; i++){
	//		if (indiv1[0][i] != localBest[0][i]){
	//			printf ( "Local Best Updatedm\n" );
	//			break; 
	//		}
	//	}		

	if ( abs(fesIndiv1Delta - fesIndiv2Delta)/max(abs(indiv1.getFitness()), abs(indiv2.getFitness())) >1e-10) {
		//		printf ( "**************** Index %d ******************************\n", indexI );
		//		printf ( "fes: indiv1 =\t\t%.16f\tindiv2 =\t%.16f\nfes: indiv1_sub =\t%.16f\tindiv2_sub =\t%.16f\n",indiv1.getFitness(),indiv2.getFitness(),indiv1_sub.getFitness(),indiv2_sub.getFitness());
		//		printf ( "fes delta: indiv1 = %.30f, indiv2 = %.30f\n", fesIndiv1Delta, fesIndiv2Delta);
		//		printf ( "the difference = %.30f, differ rate = %.10f, differ rate = %.10f\n", fesIndiv1Delta - fesIndiv2Delta, abs((fesIndiv1Delta - fesIndiv2Delta)/max(abs(indiv1.getFitness()), abs(indiv2.getFitness()))), abs(fesIndiv1Delta - fesIndiv2Delta)/max(abs(indiv1.getFitness()), abs(indiv2.getFitness())) );
		//		printf ( "*******************************************************\n\n" );
		return true; 
	}else{
		return false; 
	}
}		/* -----  end of function CCVIL::TestInterWalk  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CCVIL::sampleLearnStage
 *  Description:  With the new sampling strategy in the learning and wanna see the 
 *  							efficacy of naive sampling approach
 * =====================================================================================
 */
void CCVIL::sampleLearnStage (  ) {
	unsigned indexI=0, indexJ=0, group1, group2; // the amount of testing the interaction
	double  randi_2,  randj_2, diff; 

	// generate groupInfo according to interPartArray
	groupInfo.clear();
	// initialize the groupInfo
	for (unsigned j = 0; j<param->dimension; j++){
		vector<unsigned> tempVec;
		tempVec.push_back(j);
		lookUpGroup[j] = j;
		groupInfo.push_back(tempVec);
	}

	// 	testTimes = MaxFitEval*(param->learnPortion)/(double)3; 
	//	printf ( "test times = %d\n", testTimes );
	//	localMaxFit = MaxFitEval*(param->learnPortion); 
	//	printf ( "Local Max Fitness Evaluation for Learning Stage = %d\n", localMaxFit );

	// each individual serves for one test
	IndividualT<double> tempIndiv(ChromosomeT<double>(param->dimension));
	tempIndiv[0].initialize(fp->getMinX(), fp->getMaxX()); 
	IndividualT<double> indiv1_1, indiv2_1, indiv1_2, indiv2_2;

	tempIndiv.setFitness( fp->compute(tempIndiv[0]) );
	fes += 1; 
	sampleInfo(tempIndiv.getFitness());
	if (tempIndiv.getFitness()<bestFit){
		bestFit = tempIndiv.getFitness(); 
	}

	// indiv0 -> geneVal1
	// indiv1 -> geneVal2
	// indiv2 -> cooperation and fitness evalution
	while (fes < upperThreshold && (fes<=lowerThreshold || groupInfo.size()<param->dimension)) {

		indiv1_1 = tempIndiv;
		indiv2_1 = tempIndiv;
		indiv1_2 = tempIndiv;
		indiv2_2 = tempIndiv;

		indexI = floor(Rng::uni()*param->dimension);
		indexJ = floor(Rng::uni()*param->dimension);

		while ( (indexI==indexJ || lookUpGroup[indexI]==lookUpGroup[indexJ])&& groupInfo.size()!=1 ){
			indexI = floor(Rng::uni()*param->dimension);
			indexJ = floor(Rng::uni()*param->dimension);
		}

		if (groupInfo.size()==1){
			printf ( "Converge to one single group\n" );
			break; 
		}

		//		printf ( "randi = %d, randj = %d\n", indexI, indexJ);


		randi_2 = Rng::uni() * (fp->getMaxX() - fp->getMinX()) + fp->getMinX(); 
		randj_2 = Rng::uni() * (fp->getMaxX() - fp->getMinX()) + fp->getMinX(); 

		indiv2_1[0][indexJ] = randj_2;
		indiv1_2[0][indexI] = randi_2;

		indiv2_2[0][indexI] = randi_2;
		indiv2_2[0][indexJ] = randj_2;

		indiv1_1.setFitness( tempIndiv.getFitness() );
		indiv2_1.setFitness( fp->compute(indiv2_1[0]) );
		indiv1_2.setFitness( fp->compute(indiv1_2[0]) );
		indiv2_2.setFitness( fp->compute(indiv2_2[0]) );

		fes+=3; 

		// for minimization problem
		if (param->learnStrategy == 4){
			if ( indiv1_1.getFitness() > indiv2_1.getFitness() ){
				tempIndiv = indiv2_1; 
				sampleInfo(indiv2_1.getFitness());
				if (tempIndiv.getFitness()<bestFit){
					bestFit = tempIndiv.getFitness(); 
				}
			} else {
				sampleInfo(indiv1_1.getFitness());
			}
		}else if (param->learnStrategy == 6){
			if ( Rng::uni()<0.5 ){
				tempIndiv = indiv2_1; 
				sampleInfo(indiv2_1.getFitness());
				if (tempIndiv.getFitness()<bestFit){
					bestFit = tempIndiv.getFitness(); 
				}
			} else {
				sampleInfo(indiv1_1.getFitness());
			}
		}

		//		if ( indiv1_2.getFitness() > indiv2_2.getFitness() ){
		//			tempIndiv = indiv2_2; 
		//		}else{
		//			tempIndiv = indiv1_2;
		//		}

		diff = (indiv1_1.getFitness()-indiv2_1.getFitness()) * (indiv1_2.getFitness()-indiv2_2.getFitness()); 

		//		printf ( "diff = %f\n", diff );

		if ( diff < 0 ){ 
			group1 = lookUpGroup[indexI];
			group2 = lookUpGroup[indexJ];
			printf ("%d\t&\t%d:\t%ld\t%d\n", indexI, indexJ, fes, groupInfo.size());
			combineGroup( group1, group2 );
		}
	}

	sortGroupInfo();

	// utilizing the information and store best found value in bestCand
	// (*bestCand)=tempIndiv; 

	(*bestCand)[0].initialize(fp->getMinX(), fp->getMaxX());

	//	printf ( "After sorting, Group info\n" );
	//	print2Dvector(groupInfo);
}		/* -----  end of function CCVIL::sampleLearnStage  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CCVIL::RandomSampleGenDef
 *  Description:  Random Sampling 
 * =====================================================================================
 */
	void
CCVIL::RandomSampleGenDef (  )
{
	unsigned indexI=0, indexJ=0; 
	double randi; 

	// generate groupInfo according to interPartArray
	groupInfo.clear();
	// initialize the groupInfo
	for (unsigned j = 0; j<param->dimension; j++){
		vector<unsigned> tempVec;
		tempVec.push_back(j);
		lookUpGroup[j] = j;
		groupInfo.push_back(tempVec);
	}


	while (fes < upperThreshold && (fes<=lowerThreshold || groupInfo.size()<param->dimension)) {

		// each individual serves for one test
		IndividualT<double> tempIndiv(ChromosomeT<double>(param->dimension));
		tempIndiv[0].initialize(fp->getMinX(), fp->getMaxX()); 
		IndividualT<double> indiv(tempIndiv);

		indexI = floor(Rng::uni()*param->dimension);
		indexJ = floor(Rng::uni()*param->dimension);

		while ( (indexI==indexJ || lookUpGroup[indexI]==lookUpGroup[indexJ])&& groupInfo.size()!=1 ){
			indexI = floor(Rng::uni()*param->dimension);
			indexJ = floor(Rng::uni()*param->dimension);
		}

		randi = Rng::uni() * (fp->getMaxX() - fp->getMinX()) + fp->getMinX(); 
		indiv[0][indexJ] = randi; 

		if(testInteraction(tempIndiv, indiv, indexI)){
			unsigned group1 = lookUpGroup[indexI];
			unsigned group2 = lookUpGroup[indexJ];
			combineGroup( group1, group2 );
			printf ("%d\t&\t%d:\t%ld\t%d\n", indexI, indexJ, fes, groupInfo.size());
		}

		if (groupInfo.size()==1){
			printf ( "Converge to one single group\n" );
			break; 
		}
	}

	(*bestCand)[0].initialize(fp->getMinX(), fp->getMaxX());
	sortGroupInfo();
}		/* -----  end of function CCVIL::RandomSampleGenDef  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CCVIL::RandomWalkGenDef
 *  Description:  Random Walk with Generalized Definition
 * =====================================================================================
 */
	void
CCVIL::RandomWalkGenDef (  )
{
	unsigned indexI=0, indexJ=0; 
	double randi; 

	// generate groupInfo according to interPartArray
	groupInfo.clear();
	// initialize the groupInfo
	for (unsigned j = 0; j<param->dimension; j++){
		vector<unsigned> tempVec;
		tempVec.push_back(j);
		lookUpGroup[j] = j;
		groupInfo.push_back(tempVec);
	}

	IndividualT<double> localBest(ChromosomeT<double>(param->dimension));
	localBest[0].initialize(fp->getMinX(), fp->getMaxX()); 

	while (fes < upperThreshold && (fes<=lowerThreshold || groupInfo.size()<param->dimension)) {

		// each individual serves for one test
		IndividualT<double> indiv(localBest);

		indexI = floor(Rng::uni()*param->dimension);
		indexJ = floor(Rng::uni()*param->dimension);

		while ( (indexI==indexJ || lookUpGroup[indexI]==lookUpGroup[indexJ])&& groupInfo.size()!=1 ){
			indexI = floor(Rng::uni()*param->dimension);
			indexJ = floor(Rng::uni()*param->dimension);
		}

		randi = Rng::uni() * (fp->getMaxX() - fp->getMinX()) + fp->getMinX(); 
		indiv[0][indexJ] = randi; 


		IndividualT<double> tempIndiv(localBest);
		if(TestInterWalk( indiv, indexI, localBest)){
			unsigned group1 = lookUpGroup[indexI];
			unsigned group2 = lookUpGroup[indexJ];
			combineGroup( group1, group2 );
			printf ("%d\t&\t%d:\t%ld\t%d\n", indexI, indexJ, fes, groupInfo.size());
		}

		if (groupInfo.size()==1){
			printf ( "Converge to one single group\n" );
			break; 
		}
	}

	(*bestCand)[0].initialize(fp->getMinX(), fp->getMaxX());
	sortGroupInfo();
}		/* -----  end of function CCVIL::RandomWalkGenDef  ----- */

/*
 * procedure of optimization stage, based on the groupInfo to group the entire population
 */
void CCVIL::optimizationStage(){
	unsigned groupAmount = groupInfo.size(), innerImprove;
	bool learnStageFlag = false;


	printf("\nOptimization Stage, groupAmount = %d\n", groupAmount);

	cycle = 0;
	//	learnStageFlag = false, generate new population with regard to the groupInfo
	popGenerate(learnStageFlag);

	popInit();
	//	popInitZeros();

	groupCR = new double[groupAmount];
	groupF = new double[groupAmount];
	failCounter = new unsigned[groupAmount];

	// initialize the control parameters for each group
	for (unsigned i = 0; i< groupAmount; i++){
		groupF[i] = 0.5;
		groupCR[i] = 0.5;
		failCounter[i] = 0;
	}

	double lastCycleBestVal = 0, improveRate = 0;

	while (fes<MaxFitEval){
		++cycle;
		//		printf("===================================================\n\n\n\n===================================================\n");
		//		printf ( "Cycle %d\n", cycle );

		//		for (unsigned i=0; i<pop.size(); i++){
		//			printf("%d:\tGroupSize = %d,\tPopSize = %d\n", i, groupInfo[i].size(), pop[i].size());
		//		}
		//		printf("F %d, Optimization Cycle =%d, GroupAmount = %d, fes = %ld, Improved FES = %f, BestVal = %.8e\n",
		//				fp->getID(), 	cycle, 			(int)groupInfo.size(), fes, improveRate, bestCand->fitnessValue());

		for (unsigned i=0; i<groupAmount && fes<MaxFitEval; i++) {
			//			printf ( "Phase = %d\n" , i);
			//			printf ( "*****************************************************\n" );
			if (failCounter[i] <= param->failThreshold){
				//only optimize on current group, if no a single improvement in the past "failThreshold" successive cycle
				innerImprove = JADECC(i,learnStageFlag);
				if (innerImprove==0){
					failCounter[i] = 0;
				}else{
					failCounter[i] += innerImprove;
				}
			}

			if (sum(failCounter,groupAmount)>=((param->failThreshold+1)*groupAmount)){
				printf ( "*** Restart as no group can be optimized ***\n" );
				//				exit(EXIT_SUCCESS);
				popSizeVary(3.0);
				popInit();
				//				(*bestCand)[0].initialize(fp->getMinX(), fp->getMaxX());

				//				printf ( "beforeInitBestCand\n" );
				//				print2Dvector(pop);
				//				printPopulation((*bestCand));

				//				initBestCand(learnStageFlag);

				//				printf ( "afterInitBestCand\n" );
				//				printPopulation((*bestCand));

				// initialize the control parameters for each group
				for (unsigned i = 0; i< groupAmount; i++){
					//					groupF[i] = 0.5;
					//					groupCR[i] = 0.5;
					failCounter[i] = 0;
				}
			}
		}

		improveRate = abs((lastCycleBestVal-bestCand->fitnessValue())/bestCand->fitnessValue());

		printf ( "Fail group counter = %d\n", countFailGroupNum() );
		printf("F %d, Optimization Cycle =%d, GroupAmount = %d, groupCR = %f, groupF = %f, fes = %ld, Improved FES = %f, BestVal = %.8e\n",
				fp->getID(), 	cycle, 	(int)groupInfo.size(), mean(groupCR,groupAmount), mean(groupF,groupAmount), fes, improveRate, bestCand->fitnessValue());

		if (cycle>1 && improveRate<0.01 && fes< MaxFitEval){
			printf ( "*** Restart as non-improvement ***\n" );
			//			exit(EXIT_SUCCESS);
			popSizeVary(3.0);
			popInit();
			//			(*bestCand)[0].initialize(fp->getMinX(), fp->getMaxX());	

			//			printf ( "beforeInitBestCand\n" );
			//			print2Dvector(pop);
			//			printPopulation((*bestCand));

			initBestCand(learnStageFlag );

			//			printf ( "afterInitBestCand\n" );
			//			printPopulation((*bestCand));

			// initialize the control parameters for each group
			for (unsigned i = 0; i< groupAmount; i++){
				//				groupF[i] = 0.5;
				//				groupCR[i] = 0.9;
				failCounter[i] = 0;
			}
		}

		lastCycleBestVal = bestCand->fitnessValue();
	}

	delete[] groupCR;
	delete[] groupF;
	delete[] failCounter;
}

/* 
 * JADECC: the internal optimizer of CCVIL
 * Optimize on one specific dimension for each run
 *
 * index: the ID of group
 */
unsigned CCVIL::JADECC(unsigned index, bool learnStageFlag){
	/*********************************************************** 
	 * LB: 	Lower Bound
	 * UB: 	Upper Bound
	 * D: 	Dimension of current group
	 * NP:	Number of Population
	 * G:	Generation Limit for optimizing current group
	 ************************************************************/

	/***************************** Parameters Setting **************************/
	unsigned D = pop[index][0][0].size(), NP = pop[index].size(), G, g, *r1, *r2;
	int LB = fp->getMinX(), UB = fp->getMaxX();
	double Fm, CRm, c = param->c, p = param->p, preBestVal, preValInBestCand = 0;
	double *F, *CR;
	// f_rec stands for the fitness value improvement, which is necessary for the adaptation strategy in rJADE
	vector<unsigned> vecIndex;
	Archive* archive = NULL;

	if (learnStageFlag == true){
		G = 1;
	} else if (pop.size() == 1){
		G = INT_MAX;
	} else{
		G = min(D+5, (unsigned)500);
	}

	//  build up the vector of index, depending on whether it is learning stage or not
	if (learnStageFlag == true){
		vecIndex.push_back(index);
	}else{
		vecIndex = groupInfo[index];
	}

	//	printf("LB = %d, UB = %d, D = %d, NP = %d, G = %d, index = %d\n", LB, UB, D, NP, G, index);
	if ( learnStageFlag == true ){
		Fm = 0.5;
		CRm = 0.95;
	} else {
		CRm = groupCR[index];
		Fm = groupF[index];
	}

	F = new double[NP];
	CR = new double[NP];

	if (param->Afactor > 0) {
		// define and initialize the archive
		archive = new Archive( (NP * param->Afactor), param->dimension);
		r1 = new unsigned[NP];
		r2 = new unsigned[NP+ archive->getCapacity()];
	} else {
		r1 = new unsigned[NP];
		r2 = new unsigned[NP];
	}

	/***************************** Population Initialization and evaluation **************************/
	PopulationT<double> parents(NP, ChromosomeT<double>(param->dimension));
	PopulationT<double> offsprings(NP, ChromosomeT<double>(param->dimension));

	parents.setMinimize();
	offsprings.setMinimize();

	for (unsigned i=0; i < parents.size(); i++){
		// run on individuals
		if (learnStageFlag == true)  {
			// For learning stage, groupInfo may merge internally, while pop always maintain one-dimensional group
			// jth dimension's group belongs to the same group as current optimized group
			for (unsigned j=0; j<parents[i][0].size(); j++){
				// run on dimensions
				if (j == index)  {
					parents[i][0][j] = (pop[index])[i][0][0];
				} else {
					parents[i][0][j] = (*bestCand)[0][j];
				}
			}
			preValInBestCand = (*bestCand)[0][index];
		} else {
			// For optimization stage, use goupInfo to enumerate all elements in the same group
			for (unsigned j=0; j<parents[i][0].size(); j++){
				// not belongs to the current optimized group
				if (lookUpGroup[j] != index){
					parents[i][0][j] = (*bestCand)[0][j];
				}

				// deal with the current optimized dimensions
				for (unsigned k=0; k< groupInfo[index].size(); k++){
					unsigned I = groupInfo[index][k];
					parents[i][0][I] = (pop[index])[i][0][k];
				}
			}
		}
	}


	for (unsigned i=0; i<parents.size(); i++){
		parents[i].setFitness( fp->compute(parents[i][0]) );
	}

	/* Print the struture of the entire population */

	//	if (cycle==5 && index ==0){
	//		printf ( "Indices for optimization in this phase\n" );
	//		printVector(vecIndex);
	//
	//		printf ( "The amount of subpopulation = %d\n", (int)pop.size() );
	//
	//		printf("The whole population\n");
	//		print2Dvector(pop);
	//
	//		printf("Best Candidate at the beginning of each phase\n");
	//		printPopulation((*bestCand));
	//
	//		printf("Parents Population\n");
	//		printPopulation(parents);
	//		printf ( "Fitness of Parents\n" );
	//		printFitness(parents);
	//	}

	unsigned bestIndex = parents.bestIndex();

	//	printf("Best Index = %d, Value =  %f\n", parents.bestIndex(), parents.best().fitnessValue());

	preBestVal = parents.best().fitnessValue();
	if (preBestVal < bestFit){
		// improve on bestFit
		bestFit = preBestVal;
	}

	g = 1;
	fes = fes + NP;
	sampleInfo(preBestVal);


	/***************************** Iterations **************************/
	while ( g<=G && fes < MaxFitEval){
		vector<double> goodCR, goodF, f_rec; 

		goodCR.clear();
		goodF.clear();
		f_rec.clear();

		offsprings = parents;
		PopulationT<double> popAll(parents);

		//		printf("Parents Population\n");
		//		printPopulation(parents);

		// Generate F according to a cauchy distribution with location parameter Fm & scale parameter 0.1
		//		printf ( "\nrandFCR, CRm = %f, Fm = %f\n", CRm, Fm );

		randFCR(NP, CRm, 0.1, Fm, 0.1, F, CR);

		//		printf ( "F array:\n" );
		//		printArray(F, NP);
		//		
		//		printf ( "CR array:\n" );
		//		printArray(CR, NP);

		//		printf ( "\ngnR1R2\n" );
		if (param->Afactor == 0){
			// without archive (it is actually a special case of JADE with archive)
			gnR1R2(NP, NP, r1, r2);
		} else {
			// with archive
			//	for (unsigned i=0; i<archive->getNP(); i++){
			//		popAll.append((*archive->getPop())[i]);
			//	}

			if (archive->getNP()>0){
				popAll.insert(popAll.size(), (*archive->getPop()));
			}

			gnR1R2(NP, (int)(NP + archive->getNP()), r1, r2);
		}

		//		printf ( "r1\n" );
		//		printArray(r1, NP);
		//		printf ( "r2\n" );
		//		printArray(r2, (int)(NP + archive->getNP()));
		//
		//		printf("offsprings\n");
		//		printPopulation(offsprings);
		//
		//		printf ( "Fitness of offsprings\n" );
		//		printFitness(offsprings);
		//		
		//		printf("popAll size = %d, dimension = %d\n", popAll.size(), popAll[0][0].size());
		//		printPopulation(popAll);
		//
		//		// Find the p-best solutions
		//		printf("Find the p-best solutions\n");
		unsigned pNP = max(round(p*NP), (double)2);
		unsigned* randIndex = new unsigned[NP];

		// the indices of individual in current generation with decreasing order
		unsigned* indBest = new unsigned[pNP]; 

		PopulationT<double> pBestIndiv(NP, ChromosomeT<double>(param->dimension));
		findPbestIndex(offsprings, pNP, indBest);
		for (unsigned i=0; i<NP; i++) {
			randIndex[i] = floor(Rng::uni()*pNP);
			pBestIndiv[i] = offsprings[indBest[randIndex[i]]];
		}

		//		printf("randIndex \n");
		//		printArray(randIndex, NP);
		//
		//		printf("indBest \n");
		//		printArray(indBest, pNP);
		//
		//		printf("pBestIndiv \n");
		//		printPopulation(pBestIndiv);


		//		printf("\nBegin Mutation\n");
		//************************************ Mutation ************************************//
		PopulationT<double> vi(NP, ChromosomeT<double>(param->dimension));
		vi = offsprings;
		for (unsigned i=0; i<NP; i++){
			// for each individual
			//			printf("random individual: r1[i] = %d, r2[i] = %d\n", r1[i], r2[i]);
			for (unsigned j=0; j<vecIndex.size(); j++){
				//				printf("vector index = %d\n", vecIndex[j]);
				// for each dimension to be optimized
				vi[i][0][vecIndex[j]] = offsprings[i][0][vecIndex[j]] + F[i] * ( pBestIndiv[i][0][vecIndex[j]] - offsprings[i][0][vecIndex[j]] + offsprings[r1[i]][0][vecIndex[j]] - popAll[r2[i]][0][vecIndex[j]] );
			}
		}

		//		printf("vi, before bound constrain \n");
		//		printPopulation(vi);

		boundConstrain(vi, offsprings, LB, UB, vecIndex);

		//		printf("vi, after bound constrain \n");
		//		printPopulation(vi);

		//************************************ Crossover ************************************//
		//		printf("\nCrossover\n");
		PopulationT<double> ui(offsprings);
		ui.setMinimize();
		for (unsigned i=0; i<NP; i++){
			unsigned randIndexCrossover = floor(Rng::uni()*vecIndex.size());
			for (unsigned j=0; j<vecIndex.size(); j++){
				//				printf ( "NP %d, D %d, Inherit mutation: ", i, j );
				if ( j == randIndexCrossover || Rng::uni() < CR[i]){
					ui[i][0][vecIndex[j]] = vi[i][0][vecIndex[j]];
				}
			}
		}

		//		printf("\nSelection\n");
		//
		//		// Parents update
		//		printf("Parents population before selection");
		//		printPopulation(parents);
		//		printf ( "Fitness of Parents\n" );
		//		printFitness(parents);
		//
		//		// archive update
		//		printf("archive population BEFORE selection");
		//		printPopulation(*(archive->getPop()));
		//		printf ( "Fitness of Parents\n" );
		//		printFitness(*(archive->getPop()));

		//************************************ Selection ************************************//
		for (unsigned i=0; i<NP; i++){
			ui[i].setFitness(fp->compute(ui[i][0]));
			double fitImprv = offsprings[i].fitnessValue() - ui[i].fitnessValue();

			if ( fitImprv > 0){
				// save CR & F & f_rec
				goodCR.push_back(CR[i]);
				goodF.push_back(F[i]);
				f_rec.push_back(fitImprv);

				// improved mutation is saved to offspring, archive the failed solution
				//				printf ( "NP %d, Improved and Updated\n", i );
				archive->addToArchive(parents[i]);
				parents.replace(i, ui[i]);
			}
		}

		//		printf("Parents population after selection");
		//		printPopulation(parents);
		//		printf ( "Fitness of Parents\n" );
		//		printFitness(parents);

		if (parents.best().getFitness() < bestFit){
			// improve on bestFit
			bestFit = parents.best().getFitness();
		}

		fes += NP;
		sampleInfo(parents.best().getFitness());

		//		printf("Update Archive...\n");
		//		printf("Remove Duplicate Elememt\n");
		archive->removeDuplicateElem();

		//		printf("Truncate Archive\n");
		archive->truncateArchive();

		//		printf("archive population AFTER selection");
		//		printPopulation(*(archive->getPop()));
		//		printf ( "Fitness of Parents\n" );
		//		printFitness(*(archive->getPop()));

		// update CRm and Fm
		//printf("update CRm and Fm, goodCR size = %d, goodF size = %d\n", (int)goodCR.size(), (int)goodF.size());

		//		/* Adaptation in original JADE */
		////		printf("CRm = %f, Fm = %f\n", CRm, Fm);
		//		if (goodCR.size()>0 && sum(goodF)>0){
		//			CRm = (1-c)*CRm + c*sum(goodCR)/goodCR.size();
		//			Fm = (1-c)*Fm + c*sum(dotMultiply(goodF, goodF))/sum(goodF);
		////			printf ( "after adaptation, CRm = %f, Fm = %f\n", CRm, Fm );
		//		}

		//		printf("rJADE adaptation: CRm = %f, Fm = %f\n", CRm, Fm);
		if (goodCR.size()>0 && sum(goodF)>0){
			/* Original Adaptation in JADE by [Zhang, Anderson. TEC 2009] */
			CRm = (1-c)*CRm + c*sum(goodCR)/goodCR.size();
			// CRm = (1-c)*CRm + c*sum(dotMultiply(f_rec, goodCR))/sum(f_rec);
			Fm = (1-c)*Fm + c*sum(dotMultiply(goodF, goodF))/sum(goodF);
			//			printf ( "after adaptation, CRm = %f, Fm = %f\n", CRm, Fm );
		}

		//	if (cycle==5 && index ==0){
		//		printf("Population after Mutation, population size of it = %d \n", ui.size());
		//		printPopulation(ui);
		//
		//		printf("Fitness of ui, bestIndex = %d\n", parents.bestIndex());
		//		printFitness(ui);
		//	}

		if (g % (20000/vecIndex.size()) == 0){
			if ( learnStageFlag == true ) {
				printf("LearningStage, GroupNum = %d, D = %d, groupIndex = %d, Cycle = %d, G = %d, fes = %ld, Best Fitness = %.8e\n", (int)groupInfo.size(), (int)vecIndex.size(), index, cycle, g, fes, parents.best().fitnessValue() );
			} else {
				printf("OptimizationStage, GroupNum = %d, D = %d, groupIndex = %d, Cycle = %d, G = %d, fes = %ld, Best Fitness = %.8e\n", (int)groupInfo.size(), (int)vecIndex.size(), index, cycle, g, fes, parents.best().fitnessValue());
			}
		}

		delete[] randIndex;
		delete[] indBest;

		bestIndex = parents.bestIndex();

		g++;
	}// iterations

	/***************************** After Iterations Process **************************/
	//	printf("goodCR size = %d, goodF size = %d\n", (int)goodCR.size(), (int)goodF.size());

	//	// update the optimized optimized dimensions to bestCand and pop simultaneously
	//	if (learnStageFlag == true) {
	//		for (unsigned j=0; j<offsprings[0][0].size(); j++){
	//			// update bestCand
	//			if (j == index) {
	//				(*bestCand)[0][j] = offsprings[bestIndex][0][j];
	//			}
	//			// update pop, the entire external population
	//			for(unsigned i=0; i<NP; i++){
	//				(pop[index])[i][0][0] = offsprings[i][0][index];
	//			}
	//		}
	//	} else {
	//		for (unsigned j=0; j< groupInfo[index].size(); j++){
	//			unsigned I = groupInfo[index][j];
	//			(*bestCand)[0][I] = offsprings[bestIndex][0][I];
	//			for (unsigned i=0; i<NP; i++){
	//				(pop[index])[i][0][j] = offsprings[i][0][I];
	//			}
	//		}
	//		groupCR[index] = CRm ;
	//		groupF[index] = Fm;
	//	}
	//
	//	//	printf("Best Candidate after update\n");
	//	//	printPopulation(*bestCand);
	//	bestCand->setFitness(offsprings.best().getFitness());


	// update the optimized optimized dimensions to bestCand and pop simultaneously
	if (learnStageFlag == true) {
		// update bestCand
		(*bestCand)[0][index] = parents[bestIndex][0][index];
		// update pop, the entire external population
		for(unsigned i=0; i<NP; i++){
			(pop[index])[i][0][0] = parents[i][0][index];
		}
	} else {
		for (unsigned j=0; j< groupInfo[index].size(); j++){
			unsigned I = groupInfo[index][j];
			(*bestCand)[0][I] = parents[bestIndex][0][I];
			for (unsigned i=0; i<NP; i++){
				(pop[index])[i][0][j] = parents[i][0][I];
			}
		}

		//		groupCR[index] = mean(CR, NP);
		//		groupF[index] = mean(F, NP);

		groupCR[index] = CRm;
		groupF[index] = Fm;
	}

	//	printf ( "================After Interactions Process================\n" );
	//
	//	printf("Parents Population\n");
	//	printPopulation(parents);
	//
	//	printf ( "Fitness of Parents, bestIndex = %d\n", parents.bestIndex() );
	//	printFitness(parents);
	//
	//	printf("The whole population\n");
	//	print2Dvector(pop);
	//
	//	printf("Best Candidate after update\n");
	//	printPopulation(*bestCand);

	bestCand->setFitness(parents.best().getFitness());

	delete[] F;
	delete[] CR;
	delete[] r1;
	delete[] r2;

	if (param->Afactor > 0){
		delete archive;
	}

	if (learnStageFlag == true){
		//	elimintation for the fitness value at the end of running 
		eliminateError(parents, index);
		bestCandIndex = parents.bestIndex();
	}


	if (learnStageFlag == false){
		// For optimization stage
		if (preBestVal - parents.best().getFitness()>0){
			//			printf ( "Inner Fitness Improved!\n" );
			return 0;
		}else{
			return 1;
		}
	}else{
		// For learning stage
		if (preValInBestCand != (*bestCand)[0][index]){
			//			printf ( "Best Candidate Update!\n" );
			return 1;
		}else{
			return 0;
		}

	}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  popSizeVary
 *  Description:  Change the structure of global population 'pop' according to the factor
 * =====================================================================================
 */
	void
CCVIL::popSizeVary ( double factor )
{
	vector<unsigned> newPopSize;

	//	compute the new size of subpopualtion size and store them in vector
	for (unsigned i=0; i<pop.size(); i++){
		newPopSize.push_back(round(factor*pop[i].size()));
	}

	pop.clear();

	//	population generation in optimization stage
	for (unsigned i=0; i<groupInfo.size(); i++){
		PopulationT<double> tempPop(newPopSize[i], ChromosomeT<double>(groupInfo[i].size()));
		tempPop.setMinimize();
		pop.push_back(tempPop);
	}
}		/* -----  end of function popSizeVary  ----- */

void CCVIL::printFitness(PopulationT<double> printPop){
	for (unsigned i=0; i < printPop.size(); i++){
		printf("Individual %d, Fitness = %f\n", i, printPop[i].fitnessValue());
	}
}

// if the boundary constraint is violated, set the value to be the middle of
// the previous value and the bound
void CCVIL::boundConstrain(PopulationT<double> &vi, PopulationT<double> offsprings, int LB, int UB, vector<unsigned> vecIndex){
	for (unsigned i=0; i<vi.size(); i++){
		for (unsigned j=0; j<vecIndex.size(); j++){
			//			printf("vi = %f, LB = %f, UB = %f\n", vi[i][0][vecIndex[j]], (double)LB, (double)UB);
			if ( vi[i][0][vecIndex[j]] < LB){
				//				printf("Smaller than lower bound\n");
				vi[i][0][vecIndex[j]] = (offsprings[i][0][vecIndex[j]] + (double)LB)/2;
			} else if (vi[i][0][vecIndex[j]]> UB){
				//				printf("Larger than upper bound\n");
				vi[i][0][vecIndex[j]] = (offsprings[i][0][vecIndex[j]] + (double)UB)/2;
			}
		}
	}
}

// Find and return the indexes of the P% best individuals
void CCVIL::findPbestIndex(PopulationT<double> inPop, unsigned pNP, unsigned* indBest){
	IndividualT<double> temp;
	unsigned iMin, i, tempI, *index = new unsigned[inPop.size()];

	for (i=0; i<inPop.size(); i++){
		index[i] = i;
	}

	for (unsigned iPos=0; iPos<pNP; iPos++){
		iMin = iPos;
		for (i=iPos+1; i<inPop.size(); i++){
			if (inPop[i].getFitness()<inPop[iMin].getFitness()){
				iMin = i;
			}
		}

		if (iMin != iPos){
			// swap individuals in population
			temp = inPop[iPos];
			inPop[iPos] = inPop[iMin];
			inPop[iMin] = temp;

			// swap index
			tempI = index[iPos];
			index[iPos] = index[iMin];
			index[iMin] = tempI;
		}
	}

	for (i=0; i<pNP; i++){
		indBest[i] = index[i];
	}

	if (indBest[0]==indBest[1]){
		cerr<<"ERROR: p best index should be distinct"<<endl;
		exit(-1);
	}

	delete[] index;
}

// combine two populations for mutation process
PopulationT<double> CCVIL::combinePopulation(PopulationT<double> p1, PopulationT<double> p2){
	unsigned D1 = p1[0][0].size(), NP1 = p1.size();
	unsigned NP2 = p2.size();

	// one-on-one copying method
	PopulationT<double> popAll(NP1+NP2, ChromosomeT<double>(D1));
	for (unsigned i=0; i<popAll.size(); i++) {
		for (unsigned j=0; j<D1; j++){
			if (i < NP1){
				popAll[i][0][j] = p1[i][0][j];
			}
			else{
				popAll[i][0][j] = p2[i-NP1][0][j];
			}
		}
	}

	/*
	// try using internal interface "append"
	PopulationT<double> popAll(0, ChromosomeT<double>(D1));
	for (unsigned i=0; i<p1.size() + p2.size(); i++){
	cout<<"size of popAll = "<<popAll.size()<<endl;
	if (i < NP1){
	popAll.append(p1[i]);
	} else {
	popAll.append(p2[i-NP1]);
	}
	}
	*/
	return popAll;
}

void CCVIL::gnR1R2(unsigned NP1, unsigned NP2, unsigned *r1, unsigned *r2){
	unsigned* r0 = new unsigned[NP1];
	bool r2ReGenerateFlag;

	//	printf("NP1 = %d, NP2 = %d\n", NP1, NP2);

	// initialize r0 array
	for (unsigned i=0; i<NP1; i++){
		r0[i] = i;
	}

	for (unsigned i=0; i<NP1; i++){
		r1[i] = floor(Rng::uni()*NP1);
		for (unsigned j=0; j<INT_MAX; j++){
			if (j > 100){
				cerr<<"Can not genrate r1 in 100 iterations"<<endl;
				exit(EXIT_FAILURE);
			}

			if (r1[i]!=r0[i]){
				// distinct, continue to genenrate next random number
				break;
			}else{
				//	printf("Re-generate NP1\n");
				r1[i] = floor(Rng::uni()*NP1);
			}
		}
		//		printf("r1[%d] = %d\t", i, r1[i]);
	}
	// eliminate the duplication

	for (unsigned i=0; i<NP2; i++){
		r2[i] = floor(Rng::uni()*NP2);
		for (unsigned j=0; j<INT_MAX ; j++){
			r2ReGenerateFlag = (r2[i]==r0[i] || r2[i]==r1[i]);
			if (j > 1000){
				cerr<<"Can not genrate r2 in 1000 iterations"<<endl;
				exit(-1);
			} 

			//	printf("Begin checking ending criterion, r2ReGenerateFlag = %d\n", r2ReGenerateFlag);
			if (r2ReGenerateFlag == false){
				// distinct
				break;
			}

			if (r2ReGenerateFlag == true){
				r2[i] = floor(Rng::uni()*NP2);
			}
		}
		//		printf("r2[%d] = %d\t", i, r2[i]);
	}
	delete[] r0;
}

double CCVIL::sum(vector<double> vec){
	double s = 0;
	for (unsigned i=0; i<vec.size(); i++){
		s += vec[i];
	}
	return s;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  sum
 *  Description:  
 * =====================================================================================
 */
	unsigned
CCVIL::sum ( unsigned* arr, unsigned N )
{
	unsigned totalSum = 0;
	for (unsigned i=0; i<N; i++){
		totalSum += arr[i];
	}
	return totalSum;
}		/* -----  end of function sum  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  mean
 *  Description:  
 * =====================================================================================
 */
	double
CCVIL::mean ( double* arr, unsigned size )
{
	double sum = 0;
	for (unsigned i=0; i<size; i++){
		sum += arr[i];
	}
	return sum/((double)size);
}		/* -----  end of function mean  ----- */

vector<double> CCVIL::dotMultiply(vector<double> v1, vector<double> v2){
	vector<double> vec;
	if (v1.size()!=v2.size()){
		printf("ERROR: Dimensions of vector are not match!");
		exit(-1);
	}else{
		for(unsigned i=0; i < v1.size(); i++){
			vec.push_back(v1[i] * v2[i]);
		}
	}
	return vec;
}

void CCVIL::printArray(double* a, unsigned D){
	for(unsigned i = 0; i<D; i++){
		printf("%f\n",a[i]);
	}
}

void CCVIL::printArray(unsigned* a, unsigned D){
	for(unsigned i = 0; i<D; i++){
		printf("%d\n",a[i]);
	}
}

void CCVIL::printPopulation(PopulationT<double> printPop){
	//	printf ( "Inside the print POP\n" );
	if (printPop.size()>0){
		printf("Dimension of printed Population = %d\n", (printPop)[0][0].size());
		for (unsigned j=0; j<printPop.size(); j++){
			for (unsigned k=0; k<printPop[j][0].size(); k++){
				printf("%f\t",(printPop)[j][0][k]);
			}
			printf("\n=============================================\n");
		}
		cout<<endl;
	}
	else
		printf("size of population is zero!");
}

void CCVIL::printPopulation(IndividualT<double> printIndiv){
	printf("Dimension of printed individual = %d\n", (printIndiv)[0].size());
	for (unsigned j=0; j<(printIndiv)[0].size(); j++){
		printf("%f\t",(printIndiv)[0][j]);
	}
	cout<<endl;
}

void CCVIL::print2Dvector(vector< PopulationT<double> > vector2D){
	for (unsigned i=0; i<vector2D.size(); i++){
		printf("********************\n");
		for (unsigned j=0; j<(vector2D[i]).size(); j++){
			for (unsigned k=0; k<(vector2D[i])[j][0].size(); k++){
				printf("%f\t",(vector2D[i])[j][0][k]);
			}
			printf("\n");
		}
	}
}

void CCVIL::print2Dvector(vector< vector<unsigned> > vector2D){
	for (unsigned i=0; i<vector2D.size(); i++){
		printf("********************\n");
		for (unsigned k=0; k<(vector2D[i]).size(); k++){
			printf("%d\t",(vector2D[i])[k]);
		}
		printf("\n");
	}
}

void CCVIL::printVector(vector<unsigned> v){
	for (unsigned i=0; i<v.size(); i++){
		printf("%d\t", v[i]);
	}
	cout<<endl;
}

void CCVIL::printVector(vector<double> v){
	for (unsigned i=0; i<v.size(); i++){
		printf("%f\t", v[i]);
	}
	cout<<endl;
}

void CCVIL::printVector(vector<bool> v){
	for (unsigned i=0; i<v.size(); i++){
		printf("%d\t", (int)v[i]);
	}
	cout<<endl;
}

void CCVIL::randFCR(unsigned NP, double CRm, double CRsigma, double Fm, double Fsigma, double* &F, double* &CR){
	Normal norRnd;
	for (unsigned i=0; i<NP; i++){
		CR[i] = CRm + CRsigma * (norRnd)(); 
		CR[i] = min( 1.0 , max( 0.0, CR[i] ) ); // truncated to [0 1]
	}

	for (unsigned i=0; i<NP; i++){
		F[i] = cauchyRnd(Fm, Fsigma);
		F[i] = min(1.0, F[i]); // truncated to [-INF, 1]
	}

	// instead of truncation for dealing with the lower bound -1, we regenerate invalid random numbers
	vector<unsigned> pos1, pos2;
	pos1.clear();
	pos2.clear();

	for (unsigned i=0; i<NP; i++){
		if (F[i]<=0){
			pos1.push_back(i);
		}
		else if (F[i]>=1){
			pos2.push_back(i);
		}
	}

	while ( !pos1.empty() || !pos2.empty() ) {

		for (unsigned i=0; i<pos1.size(); i++){
			//			F[pos1[i]] = cauchyRnd(Fm, Fsigma);
			//			F[pos1[i]] = min(1.0, F[pos1[i]]); // truncated to [-INF, 1]
			F[pos1[i]] = min(1.0, cauchyRnd(Fm, Fsigma)); // truncated to [-INF, 1]
		}

		for (unsigned i=0; i<pos2.size(); i++){
			//			F[pos2[i]] = cauchyRnd(Fm, Fsigma);
			//			F[pos2[i]] = min(1.0, F[pos2[i]]); // truncated to [-INF, 1]
			F[pos2[i]] = min(1.0, cauchyRnd(Fm, Fsigma)); // truncated to [-INF, 1]
		}

		//		printf ( "size of pos1 = %d, pos2 = %d, before clear\n", size1, size2);
		pos1.clear();
		pos2.clear();

		for (unsigned i=0; i<NP; i++){
			if (F[i]<=0){
				pos1.push_back(i);
			}
			else if (F[i]>=1){
				pos2.push_back(i);
			}
		}
	}

	for (unsigned i=0; i<NP; i++){
		if (F[i]>=1){
			cerr<<"ERROR: F = 1"<<endl;
			exit(-1);
		}
	}
}

double CCVIL::cauchyRnd(double mu, double delta){
	Uniform uniRnd;
	return (mu + delta * tan(PI * (uniRnd()-0.5)));
	//	return (mu + delta * tan(PI * (uniRnd())));
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  popGenerate
 *  Description:  generate the basic structure of the entire population, the pre-process
 *  for population initialization
 * =====================================================================================
 */
void CCVIL::popGenerate(bool learnStageFlag){

	pop.clear();

	if (learnStageFlag == true){
		for (unsigned i=0; i<param->dimension; i++){
			PopulationT<double> tempPop(param->NP, ChromosomeT<double>(param->initialGroupSize));
			tempPop.setMinimize();
			pop.push_back(tempPop);
		}
	}else{
		//		population generation in optimization stage
		for (unsigned i=0; i<groupInfo.size(); i++){
			unsigned NP = groupInfo[i].size()+10;
			//TODO: Change the population back to +10
			//			unsigned NP = groupInfo[i].size()+2;
			PopulationT<double> tempPop(NP, ChromosomeT<double>(groupInfo[i].size()));
			tempPop.setMinimize();
			pop.push_back(tempPop);
		}
	}
}

/*
 * Initialize the population
 */
void CCVIL::popInit(){
	for (unsigned i=0; i<pop.size(); i++){
		for(unsigned j=0; j<pop[i].size(); j++){
			(pop[i])[j][0].initialize(fp->getMinX(), fp->getMaxX());
		}
	}
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  popInitZeros
 *  Description:  Initial all genes with 0
 *   Dependency:  popGenerate
 * =====================================================================================
 */
void CCVIL::popInitZeros ()
{
	for (unsigned i=0; i<pop.size(); i++){
		for(unsigned j=0; j<pop[i].size(); j++){
			(pop[i])[j][0].initialize(fp->getMinX(), fp->getMaxX());

			for (unsigned k=0; k<pop[i][0].size(); k++){
				(pop[i])[j][0][k] = 0;
			}
		}
	}

	//	printf ( "Print the entire population\n" );
	//	print2Dvector(pop);


	for (unsigned i = 0; i < (*bestCand)[0].size(); i++){
		(*bestCand)[0][i] = 0;
	}

	bestCand->setFitness(fp->compute((*bestCand)[0]));

	printf ( "Fitness of bestCand = %f\n", bestCand->getFitness());

	//	printf ( "Print the best Candidate\n" );
	//	printPopulation(*bestCand);

}		/* -----  end of function popInitZeros  ----- */

void CCVIL::captureInter(unsigned curDim, unsigned lastDim){

	unsigned NP = pop[lastDim].size(), randi, counter=0;
	IndividualT<double> randIndiv = (*bestCand);
	//	printf("Capture Interaction Between %d & %d\n", curDim, lastDim);
	//	printf("The population  size of last optimized dimension = %d\n", NP);

	while (true) {
		randi = floor(rand() / ((double) RAND_MAX) * NP);
		if ( pop[lastDim][randi][0][0] != (*bestCand)[0][lastDim] ){
			break;
		} else if ( counter > 1000 ) {
			printf("ERROR: Fail to generate random index within 1000 trying\n");	
		}
	}

	randIndiv[0][lastDim] = pop[lastDim][randi][0][0];
	randIndiv.setFitness(fp->compute(randIndiv[0]));
	fes++;

	//		printf ( "bestCand\n" );
	//		printPopulation(*bestCand);
	//
	//		printf ( "randIndiv\n" );
	//		printPopulation(randIndiv);
	//
	//		printf ( "Fitness of randindiv == %f\n", randIndiv.getFitness());
	//		printf ( "Fitness of bestCand == %f\n", (*bestCand).getFitness());

	// if there is any interaction detected, combine the groupInfo
	if (randIndiv.getFitness()<bestCand->getFitness()){
		//			printf("Interaction Detected between %d & %d\n", curDim, lastDim);
		unsigned group1, group2;

		//		vector<unsigned> rmVec;
		//		printf("============= Before Merge =============\n");
		//		printf("Look Up Group Table:\n");
		//		printArray(lookUpGroup, param->dimension);

		group1 = lookUpGroup[curDim];
		group2 = lookUpGroup[lastDim];

		//		print2Dvector(groupInfo);
		//		printf("groupInfo:\n");
		//		printf("group1 = %d, group2 =%d\ncurrent D = %d, last D = %d\n", group1, group2, curDim, lastDim)
		// through comparison, always join the latter one into the previous
		if(group1 < group2){
			combineGroup(group1, group2);

			//				// join group2 into group1
			//				for (unsigned i=0; i<groupInfo[group2].size(); i++){
			//					// instead of push at the back of vector
			//					groupInfo[group1].push_back(groupInfo[group2][i]);
			//				}
			//				//			groupInfo[group2].clear();
			//
			//				for (unsigned i = 0; i<groupInfo[group2].size(); i++){
			//					lookUpGroup[groupInfo[group2][i]] = group1;
			//				}
			//
			//				groupInfo.erase(groupInfo.begin() + group2);
			//				for (unsigned i=0; i<param->dimension; i++){
			//					if (lookUpGroup[i]>group2){
			//						lookUpGroup[i]--;
			//					}
			//				}
		}else{// group2 < group1
			combineGroup(group2, group1);

			//				// join group1 into group2
			//				// rmVec = groupInfo[group1];
			//				for (unsigned i=0; i<groupInfo[group1].size(); i++){
			//					groupInfo[group2].push_back(groupInfo[group1][i]);
			//					//		groupInfo[group2].push_back(rmVec[i]);
			//				}
			//				//	groupInfo[group1].clear();
			//
			//				for (unsigned i = 0; i<groupInfo[group1].size(); i++){
			//					lookUpGroup[groupInfo[group1][i]] = group2;
			//				}
			//
			//				groupInfo.erase(groupInfo.begin() + group1);
			//
			//				for (unsigned i=0; i<param->dimension; i++){
			//					if (lookUpGroup[i]>group1){
			//						lookUpGroup[i]--;
			//					}
			//				}
		}

		//		printf("============= After Merge =============\n");
		//		printf("Look Up Group Table:\n");
		//		printArray(lookUpGroup, param->dimension);
		//		print2Dvector(groupInfo); 
		//		printf("groupInfo:\n");
		//		printf("group1 = %d, group2 =%d\ncurrent D = %d, last D = %d\n", lookUpGroup[curDim], lookUpGroup[lastDim], curDim, lastDim);

	}
	//		else if (randIndiv.getFitness() == bestCand->getFitness()){
	//			printf ( "bestCand\n" );
	//			printPopulation(*bestCand);
	//
	//			printf ( "randIndiv\n" );
	//			printPopulation(randIndiv);
	//
	//			printf ( "Fitness of randindiv == %f\n", randIndiv.getFitness());
	//			printf ( "Fitness of bestCand == %f\n", (*bestCand).getFitness());
	//
	//			printf ( "Error: Generating Random Individual, it is the same with the bestCand\n" );
	//			exit(EXIT_SUCCESS);
	//		}
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CCVIL::combineGroup
 *  Description:  
 * =====================================================================================
 */
void CCVIL::combineGroup ( unsigned group1, unsigned group2 )
{
	for (unsigned i=0; i<groupInfo[group2].size(); i++){
		// instead of push at the back of vector
		groupInfo[group1].push_back(groupInfo[group2][i]);
	}

	for (unsigned i = 0; i<groupInfo[group2].size(); i++){
		lookUpGroup[groupInfo[group2][i]] = group1;
	}

	groupInfo.erase(groupInfo.begin() + group2);
	for (unsigned i=0; i<param->dimension; i++){
		if (lookUpGroup[i]>group2){
			lookUpGroup[i]--;
		}
	}
	groupRec.push_back(groupInfo.size());
	groupFesRec.push_back(fes);
}		/* -----  end of function CCVIL::combineGroup  ----- */

/*
 * check whether the two variable belong to the same group or not
 */
bool CCVIL::sameGroup(unsigned v1, unsigned v2){
	if (lookUpGroup[v1] == lookUpGroup[v2]){
		return true;
	}else{
		return false;
	}
}

/*
 * generate an random permutation with the length N
 */
unsigned* CCVIL::randPerm(unsigned N)
{
	unsigned* p = new unsigned[N];
	for (unsigned i = 0; i < N; ++i) {
		int j = rand() % (i + 1);
		p[i] = p[j];
		p[j] = i;
	}
	return p;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  itos
 *  Description:  
 * =====================================================================================
 */
	string
CCVIL::itos ( int i )
{
	stringstream s;
	s << i;
	return s.str();
}		/* -----  end of function itos  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  sampleInfo
 *  Description:  
 * =====================================================================================
 */
	void
CCVIL::sampleInfo ( double curFit )
{
	if (!samplingPoints.empty()&& fes>=samplingPoints.front()){
		samplingPoints.erase(samplingPoints.begin());
		fesRec.push_back(fes);
		valRec.push_back(curFit);
	}
}
/* -----  end of function sampleInfo  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  initBestCand
 *  Description:  initialize the best candidate according to global population
 * =====================================================================================
 */
	void
CCVIL::initBestCand (bool learnStageFlag )
{
	unsigned curDim=0, randi = 0;
	// assign bestCand from pop in a group by group fashion
	if (learnStageFlag == false){
		for (unsigned i=0; i<groupInfo.size(); i++){
			for (unsigned j=0; j<groupInfo[i].size(); j++){
				curDim = groupInfo[i][j];
				randi = Rng::uni()*pop[i].size();
				(*bestCand)[0][curDim] = pop[i][randi][0][j] ;
			}
		}
	}else{
		for (unsigned i=0; i<pop.size(); i++){
			randi = floor(rand() / ((double) RAND_MAX) * pop[i].size());
			(*bestCand)[0][i] = pop[i][randi][0][0];
		}
	}
}		/* -----  end of function initBestCand  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  eliminateError
 *  Description:  
 * =====================================================================================
 */
void
CCVIL::eliminateError ( PopulationT<double> population, unsigned index){
	double bestFitVal = population.best().getFitness();
	for (unsigned i=0; i<population.size(); i++){
		if (bestFitVal == population[i].getFitness() && i!=population.bestIndex()){
			impreciseGroup[index] = true;
			break;
		}
	}
}		/* -----  end of function eliminateError  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  sortGroupInfo
 *  Description:  sort the 2D vector groupInfo
 * =====================================================================================
 */
	void
CCVIL::sortGroupInfo (  )
{
	for (unsigned i=0; i<groupInfo.size(); i++){
		if (groupInfo[i].size()>1){
			sort(groupInfo[i].begin(), groupInfo[i].end());
		}
	}
}		/* -----  end of function sortGroupInfo  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  countFailGroupNum
 *  Description:  
 * =====================================================================================
 */
	unsigned
CCVIL::countFailGroupNum (  )
{
	unsigned sum = 0;
	for (unsigned i=0; i<groupInfo.size(); i++){
		if(failCounter[i]>param->failThreshold){
			sum += 1;
		}
	}
	return sum;
}		/* -----  end of function countFailGroupNum  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  CCVIL::getPriorInterStage
 *  Description:  this stage 1) generate the groupInfo; 2) initialize bestCand 
 * =====================================================================================
 */
	void
CCVIL::getPriorInterStage ()
{
	unsigned arrSize = param->dimension*(param->dimension-1)/2;
	/**********************************
	 * 1) generate groupInfo
	 **********************************/

	// generate random permutation
	unsigned* randPermInter = randPerm(arrSize);
	//			for (unsigned i=0; i<arrSize; i++){
	//				printf ( "%d\t", (int)randPermInter[i] );
	//			}
	//			printf( "\n" );

	// firstly suppose that there is no prior information
	bool* interPartArray = new bool[arrSize];
	for (unsigned i=0; i<arrSize; i++){
		if (param->learnStrategy == 1){
			interPartArray[i] = false;
		}else if (param->learnStrategy == 2){
			interPartArray[i] = true;
		}else if (param->learnStrategy == 3){
			double randN = rand()/(double)RAND_MAX;
			if (randN<0.5){
				//						printf ( "F\n" );
				interPartArray[i] = false;
			}else {
				//						printf ( "T\n" );
				interPartArray[i] = true;
			}
		}else{
			printf ( "Learning Strategy not found\n" );
			exit(EXIT_FAILURE);
		}
	}	


	// truncate the known ProirInterStage according to given percentage of prior interaction information
	unsigned knownArrSize = floor(arrSize* param->knownGroupPercent[0]);
	printf ( "Known Array Size = %d, Interaction Array Size = %d\n", knownArrSize, fp->getInterArray().size() );

	for (unsigned i=0; i< knownArrSize; i++){
		interPartArray[randPermInter[i]] = (fp->getInterArray())[randPermInter[i]];
		//				interPartArray[i] = (fp->getInterArray())[i];
	}


	//			printf ( "Interaction Complete Information\n" );
	//			for (unsigned i=0; i<arrSize; i++){
	//				printf ( "%d\t", (int)fp->getInterArray()[i] );
	//			}
	//			printf( "\n" );
	//
	//			printf ( "Interaction Partial Information\n" );
	//			for (unsigned i=0; i<arrSize; i++){
	//				printf ( "%d\t", interPartArray[i] );
	//			}
	//			printf( "\n" );

	// generate groupInfo according to interPartArray
	groupInfo.clear();
	// initialize the groupInfo
	for (unsigned j = 0; j<param->dimension; j++){
		vector<unsigned> tempVec;
		tempVec.push_back(j);
		lookUpGroup[j] = j;
		groupInfo.push_back(tempVec);
	}

	// combine interactive groups as if the information is caputred by learning
	unsigned group1=-1, group2=-1;
	for (unsigned i=0; i<arrSize; i++){
		if (interPartArray[i]==true ){
			unsigned I1=-1, I2=-1;
			fp->MatToArr(I1, I2, i);
			group1 = lookUpGroup[I1];
			group2 = lookUpGroup[I2];
			if (group1 != group2){
				//						printf ( "Combine Index: %d & %d\n", I1, I2 );
				combineGroup(lookUpGroup[I1], lookUpGroup[I2]);

				//						print2Dvector(groupInfo);
			}
		}
	}

	sortGroupInfo();
	//			printf ( "After sorting, Group info\n" );
	//			print2Dvector(groupInfo);
	printf ( "The amount of groups = %d\n", groupInfo.size() );

	delete interPartArray;
	delete randPermInter;
}		/* -----  end of function CCVIL::getPriorInterStage  ----- */
