#include "Archive.h"

Archive::Archive(unsigned inNP, unsigned inD){
	MAX_NP = inNP; 
	dimension = inD;
	pop = new PopulationT<double>(0, ChromosomeT<double>(dimension));
}

Archive::~Archive(){
	delete pop;
}

unsigned Archive::getNP(){
	return (*pop).size();
}

PopulationT<double> *Archive::getPop(){
	return pop;
}

/*
 * Update the archive with input solutions
    Step 1: Add new solution to the archive
    Step 2: Remove duplicate elements
    Step 3: If necessary, randomly remove some solutions to maintain the archive size
 */

void Archive::addToArchive(IndividualT<double> failIndiv){
	// Step 1: Add new solution to the archive
	pop->append(failIndiv); 
}

void Archive::removeDuplicateElem(){
	// Step 2: Remove duplicate elements
	if (pop->size() > 1){
		unsigned D = (*pop)[0][0].size();
		bool duplicateIndiv; 

		for (unsigned i=0; i<pop->size(); i++){
			for(unsigned j=pop->size()-1; j>=i+1; j--){

				duplicateIndiv = true;

				for (unsigned k=0; k<D; k++){
					if ((*pop)[j][0][k] != (*pop)[i][0][k]){
						duplicateIndiv = false;
						break;
					}
				}

				if (duplicateIndiv == true){
					pop->remove(j);
				}

			}
		}
	}
}

void Archive::truncateArchive(){
	// Step 3: If necessary, randomly remove some solutions to maintain the archive size
	if (pop->size()>MAX_NP){
		unsigned numOfElemToRm =  pop->size() - MAX_NP;
//		printf("Number of Element to remove = %d\n", numOfElemToRm);
		for (unsigned i=0; i< numOfElemToRm; i++){
			//			printf(", popSize = %d\n", pop->size());
			pop->remove( floor(Rng::uni()*(pop->size())) );
			//			printf("After Remove, popSize = %d\n", pop->size());
		}
	}
}

/*
PopulationT<double> Archive::unique(PopulationT<double> popAll){
	if (popAll.size() > 0){
		unsigned D = popAll[0][0].size();
		PopulationT<double> uniPopAll(0, ChromosomeT<double>(D));
		bool alreadyExist = false;

		for (unsigned i=0; i<popAll.size(); i++){
			// check whether popAll[i] already has been in uniPopAll
			alreadyExist = true;

			for(unsigned j=0; j<uniPopAll.size(); j++){
				for (unsigned k=0; k<D; k++){
					if (uniPopAll[j][0][k] != popAll[i][0][k]){
						alreadyExist = false;
						break;
					}
				}
			}

			if (alreadyExist == false || uniPopAll.size()==0){
				uniPopAll.append(popAll[i]);
			}
		}

		//	printf("Size of popAll = %d, size of uniPopAll = %d\n", popAll.size(), uniPopAll.size());
		return uniPopAll;
	}else{
		return popAll;
	}
}
*/

/* 
 * ===  FUNCTION  ======================================================================
 *
 *         Name:  Archive::getCapacity()
 *  Description:  return the maximun size of archive
 * =====================================================================================
 */
	unsigned	
Archive::getCapacity() 
{
	return MAX_NP;
}		/* -----  end of function Archive::getCapacity()  ----- */
