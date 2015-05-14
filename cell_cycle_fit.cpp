/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#include "ccfit.h"

int main(int argc, char *argv[]) {

	// read arguments
	ProgramOptionParser *opt = new ProgramOptionParser("ccfit");
	opt->parseArgs(argc,argv);

	// initialize structure to store simulated data
	FitData *fD = new FitData;

	// run simulations for all chase times
	for (int i=0;i<opt->getNumChaseTimes();i++) {
		ccfit(argc,argv,fD,i);
	}

	fD->serializeOutput();
	delete fD; fD=NULL;
	delete opt; opt=NULL;
	return 0;
}
