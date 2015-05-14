/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#include "ccfit.h"
using namespace std;

void ccfit(int argc, char *argv[], FitData *fD, int chaseTimeIdx) {
	ProgramOptionParser *opt = new ProgramOptionParser("ccfit");
	opt->parseArgs(argc,argv);
	opt->setMinTime(chaseTimeIdx);

	// Get random seed
	base_generator_type *generator = new base_generator_type;
	initializeRandomNumberGenerator(opt->randomSeed, generator);
	int randomInt = randUniformInt(0,INT_MAX-opt->numSim,generator);
	delete generator; generator=NULL;

	// initialize storage structures
	GeometryProfile gg;
	opt->getGeometryProfile(&gg);
	int totalNumberCells = gg.totalNumberOfCells(0);
	Array2D<int> *position = new Array2D<int>;
	Array2D<float> *dna = new Array2D<float>;
	Array2D<int> *edu = new Array2D<int>;
	Array2D<int> *mPhase = new Array2D<int>;
	position->resize(opt->numSim,totalNumberCells);
	dna->resize(opt->numSim,totalNumberCells);
	edu->resize(opt->numSim,totalNumberCells);
	mPhase->resize(opt->numSim,totalNumberCells);

	Array2D<float> *G1idx = new Array2D<float>;
	Array2D<float> *Sidx =  new Array2D<float>;
	Array2D<float> *G2idx = new Array2D<float>;
	Array2D<float> *Midx =  new Array2D<float>;
	G1idx->resize(opt->getNumCellRows(),opt->numSim);
	Sidx->resize(opt->getNumCellRows(),opt->numSim);
	G2idx->resize(opt->getNumCellRows(),opt->numSim);
	Midx->resize(opt->getNumCellRows(),opt->numSim);

	// run simulations
	#if defined USE_LIBDISPATCH
	dispatch_apply((size_t) opt->numSim, dispatch_get_global_queue(0, 0), ^(size_t ii) {
		int i=(int)ii;
	#else
	for (int i=0;i<opt->numSim;i++) {
	#endif
		// Read in program options (again) to ensure thread safety
		ProgramOptionParser *opt_thread = new ProgramOptionParser("ccfit");
		opt_thread->parseArgs(argc,argv);
		opt_thread->setMinTime(chaseTimeIdx);

		// run germline simulation
		GermLine *gL = new GermLine(opt_thread,randomInt+i);
		gL->run("filled",200);
		int counter=0;
		for (unsigned int j=0;j<gL->germCellVector_->size();j++) {
			if ((*gL->germCellVector_)[j]->mitoticOrMeioticOrLaid()==0) {
				(*dna)(i,counter) = (*gL->germCellVector_)[j]->getDna();
				(*edu)(i,counter) = (*gL->germCellVector_)[j]->edu();
				(*mPhase)(i,counter) = (*gL->germCellVector_)[j]->isMphase();
				int spatialNoise = 0;
				(*position)(i,counter) = spatialNoise + (*gL->germCellVector_)[j]->getX();
				counter++;
			}
		}
		// get cell cycle phase indices at the end of the simulation
		Array1D<float> G1,S,G2,M;
		gL->getSpatialCellCyclePhaseIndices(&G1,&S,&G2,&M);
		for (int j=0;j<opt->getNumCellRows();j++) {
			(*G1idx)(j,i) = G1(j);
			(*Sidx)(j,i) = S(j);
			(*G2idx)(j,i) = G2(j);
			(*Midx)(j,i) = M(j);
		}

		delete gL; gL=NULL;
		delete opt_thread; opt_thread=NULL;
	#if defined USE_LIBDISPATCH
	});
	#else
	}
	#endif

	Array1D<float> dnaBin, chaseTimes;
	opt->getDnaBin(&dnaBin);
	opt->getChaseTimes(&chaseTimes);
	fD->initializeStorage(opt->getNumCellRows(),&dnaBin,&chaseTimes);
	fD->parseData(chaseTimeIdx,position,dna,edu,mPhase);

	if (opt->verbose) {
		printCellPhaseIndices(G1idx,Sidx,G2idx,Midx);
	}

	delete G1idx; G1idx=NULL;
	delete Sidx;  Sidx=NULL;
	delete G2idx; G2idx=NULL;
	delete Midx;  Midx=NULL;
	delete position; position=NULL;
	delete dna; dna=NULL;
	delete edu; edu=NULL;
	delete mPhase; mPhase=NULL;
	delete opt; opt=NULL;
}
