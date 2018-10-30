#ifndef _ARCHIVE_H
#define _ARCHIVE_H

#include <EALib/PopulationT.h>
#include <Rng/Uniform.h>

class Archive{
protected:
	unsigned MAX_NP;
	unsigned dimension;
	PopulationT<double>* pop;

public:
	Archive(unsigned inNP, unsigned inD);
	~Archive();

	unsigned getNP();
	unsigned getCapacity();
	void updateArchive(PopulationT<double> popAll);
	void addToArchive(IndividualT<double> failIndiv);
	void removeDuplicateElem();
	void truncateArchive();
	//	PopulationT<double> unique(PopulationT<double> popAll);
	PopulationT<double> *getPop();
};

#endif
