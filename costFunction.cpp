/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#include "costFunction.h"
#include <pthread.h>

#ifdef USE_PTHREADS
#define MAX_NUM_THREADS 64
#endif
static int randomNumberGeneratorSeed = 0;

typedef struct thread_data_t {
    _Atomic(int)* taskCounter;
    PedigreeDepthStorage *pdStorage;
    ProgramOptionParser *opt;
    int argc;
    char **argv;
    int randomInt;
    int i;
} thread_data_t;

typedef enum memory_order {
    memory_order_relaxed, memory_order_consume, memory_order_acquire,
    memory_order_release, memory_order_acq_rel, memory_order_seq_cst
} memory_order;

void * run_thread(void *arg) {
    thread_data_t *data = (thread_data_t *) arg;
    
    PedigreeDepthStorage *pdStorage = data -> pdStorage;
    int randomInt = data -> randomInt;
    //int i = data -> i;
    
    int task_id;
    
    while ( (task_id = __c11_atomic_fetch_add(data-> taskCounter, 1, memory_order_seq_cst)) < data -> opt -> numSim) {
    
    // Read in program options (again) to ensure thread safety
    ProgramOptionParser *opt_thread = new ProgramOptionParser("pd");
    opt_thread -> parseArgs(data -> argc, (char **) data -> argv);
    
    // run a single germline simulation
    GermLine *gL = new GermLine(opt_thread,randomInt+task_id);
	if (opt_thread->filledGermline) { // default value is false
		gL->run("filled",0);
	} else {
		gL->run("oneCell",0);
	}

    
    // calculate fitness metrics
    float pdNCellsEnterMeioticRegion, pdSperm, pdOocyte, pdSpermAndOocyte, pdSpermAndOocyteSum, pdMitoticRegion, pdMeioticRegion, pdWholeGermline, pdGermlineAtEggLaying,timeFirstCellExitMeioticRegion,timeNthCellExitsMitoticRegion;
    int numCellsEnterMeioticRegion, numCellDivisions;
    gL->calculateOutput(pdNCellsEnterMeioticRegion, pdSperm,pdOocyte,pdSpermAndOocyte,pdSpermAndOocyteSum, pdMitoticRegion, pdMeioticRegion, pdWholeGermline, pdGermlineAtEggLaying,numCellsEnterMeioticRegion,numCellDivisions,timeFirstCellExitMeioticRegion,timeNthCellExitsMitoticRegion);
    int idx = int(task_id);
    (*pdStorage->pd_sperm)(idx) = pdSperm;
    (*pdStorage->pd_oocyte)(idx) = pdOocyte;
    (*pdStorage->pd_spermOocyte)(idx) = pdSpermAndOocyte;
    (*pdStorage->pd_mitoticRegion)(idx) = pdMitoticRegion;
    (*pdStorage->pd_meioticRegion)(idx) = pdMeioticRegion;
    (*pdStorage->pd_wholeGermline)(idx) = pdWholeGermline;
    (*pdStorage->pd_spermOocyteSum)(idx) = pdSpermAndOocyteSum;
    (*pdStorage->pd_GermlineAtEggLaying)(idx) = pdGermlineAtEggLaying;
    (*pdStorage->pd_AllCellsEnterMeioticRegion)(idx) = pdNCellsEnterMeioticRegion;
    (*pdStorage->num_meioticCellsProduced)(idx) = (float)numCellsEnterMeioticRegion;
    (*pdStorage->num_cellDivisions)(idx) = (float)numCellDivisions;
    (*pdStorage->time_firstMeioticExit)(idx) = timeFirstCellExitMeioticRegion;
    (*pdStorage->time_nthMitoticExit)(idx) = timeNthCellExitsMitoticRegion;
    
    // calculate schedules
    gL->calculateSchedules(idx,pdStorage);
    // clean up
    delete gL; gL=NULL;
    delete opt_thread; opt_thread=NULL;
        
    }

    
    pthread_exit(NULL);
}

/*
 * \brief Pedigree depth simulations
 * This function runs pedigree depth simulations based on parameters specified by user.
 * Return is the value of the fitness for user-specified parameters.  The lower the fitness metric, the better
 * Return value is -1 if constraints are not satisfied.
 */
extern "C"
float costFunction_pedigreeDepth(int argc, char *argv[], float resultStats[]) {

	// read user-specified program options
	ProgramOptionParser *opt = new ProgramOptionParser("pd");
	opt->parseArgs(argc,argv);

	// Get random seed
	base_generator_type *generator = new base_generator_type;
	if (opt->randomSeed <= 0 && randomNumberGeneratorSeed > 0) {
		//We have already generated a random seed for this executable using the
		//system entropy source; reuse it
		initializeRandomNumberGenerator(randomNumberGeneratorSeed, generator);
	} else {
		initializeRandomNumberGenerator(opt->randomSeed, generator);
	}
	randomNumberGeneratorSeed = randUniformInt(1, INT_MAX, generator);
	int randomInt = randUniformInt(0,INT_MAX-opt->numSim,generator);
	delete generator; generator=NULL;

	// initialize structures to store average pedigree depth for sperm, oocytes, MR, etc.
	PedigreeDepthStorage *pdStorage = new PedigreeDepthStorage(opt);

	if (mitoticRegionTooBig(opt)) {
		if (opt->verbose) {printf("Mitotic region too big\n");}
		if (opt->verbose_tedPedigreeDepth) {
			printSimulationResultsForTed_pedigreeDepth_badMR(pdStorage, opt);
		}
		delete opt; opt=NULL;
		delete pdStorage; pdStorage=NULL;
		return -1;
	}
    
    #if defined USE_PTHREADS
    
    int num_threads = std::min((int) (size_t)opt->numSim, MAX_NUM_THREADS);
    _Atomic(int) taskCounter = 0;
    
    pthread_t *thr = (pthread_t *) malloc (sizeof (pthread_t) * (unsigned long) num_threads);
    thread_data_t *thr_data = (thread_data_t *) malloc (sizeof (thread_data_t) * (unsigned long) num_threads);
    
    for (int i = 0; i < num_threads ; i++) {
        thr_data[i].pdStorage = pdStorage;
        thr_data[i].opt = opt;
        thr_data[i].argc = argc;
        thr_data[i].argv = argv;
        thr_data[i].randomInt = randomInt;
        thr_data[i].i=i;
        thr_data[i].taskCounter = &taskCounter;
        int rc;
        
        if ((rc = pthread_create(&thr[i], NULL, run_thread, &thr_data[i]))) {
            fprintf(stderr, "error: pthread_create, rc: %d\n", rc);
            return EXIT_FAILURE;
        }
    }
    
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(thr[i], NULL);
    }

    free(thr); thr = NULL;
    free(thr_data); thr_data = NULL;
    
    #else

	// run simulations
	#if defined USE_LIBDISPATCH
	dispatch_apply((size_t)opt->numSim, dispatch_get_global_queue(0, 0), ^(size_t ii) {
	int i=(int)ii;
	#else
	for (int i=0;i<opt->numSim;i++) {
	#endif
		// Read in program options (again) to ensure thread safety
		ProgramOptionParser *opt_thread = new ProgramOptionParser("pd");
		opt_thread->parseArgs(argc,argv);

		// run a single germline simulation
		GermLine *gL = new GermLine(opt_thread,randomInt+i);
		if (opt_thread->filledGermline) { // default value is false
			gL->run("filled",0);
		} else {
			gL->run("oneCell",0);
		}

		// calculate fitness metrics
		float pdNCellsEnterMeioticRegion, pdSperm, pdOocyte, pdSpermAndOocyte, pdSpermAndOocyteSum, pdMitoticRegion, pdMeioticRegion, pdWholeGermline, pdGermlineAtEggLaying,timeFirstCellExitMeioticRegion,timeNthCellExitsMitoticRegion;
		int numCellsEnterMeioticRegion, numCellDivisions;
		gL->calculateOutput(pdNCellsEnterMeioticRegion, pdSperm,pdOocyte,pdSpermAndOocyte,pdSpermAndOocyteSum, pdMitoticRegion, pdMeioticRegion, pdWholeGermline, pdGermlineAtEggLaying,numCellsEnterMeioticRegion,numCellDivisions,timeFirstCellExitMeioticRegion,timeNthCellExitsMitoticRegion);
		int idx = int(i);
		(*pdStorage->pd_sperm)(idx) = pdSperm;
		(*pdStorage->pd_oocyte)(idx) = pdOocyte;
		(*pdStorage->pd_spermOocyte)(idx) = pdSpermAndOocyte;
		(*pdStorage->pd_mitoticRegion)(idx) = pdMitoticRegion;
		(*pdStorage->pd_meioticRegion)(idx) = pdMeioticRegion;
		(*pdStorage->pd_wholeGermline)(idx) = pdWholeGermline;
		(*pdStorage->pd_spermOocyteSum)(idx) = pdSpermAndOocyteSum;
		(*pdStorage->pd_GermlineAtEggLaying)(idx) = pdGermlineAtEggLaying;
		(*pdStorage->pd_AllCellsEnterMeioticRegion)(idx) = pdNCellsEnterMeioticRegion;
		(*pdStorage->num_meioticCellsProduced)(idx) = (float)numCellsEnterMeioticRegion;
		(*pdStorage->num_cellDivisions)(idx) = (float)numCellDivisions;
		(*pdStorage->time_firstMeioticExit)(idx) = timeFirstCellExitMeioticRegion;
		(*pdStorage->time_nthMitoticExit)(idx) = timeNthCellExitsMitoticRegion;

		// calculate schedules
		gL->calculateSchedules(i,pdStorage);
		// clean up
		delete gL; gL=NULL;
		delete opt_thread; opt_thread=NULL;
	#if defined USE_LIBDISPATCH
	});
	#else
	}
	#endif
                   
    #endif

	// get return value
	resultStats[0] = calculateFitness(opt->fitnessMetric.c_str(),pdStorage,opt);

	// get whether constraints are met
	bool cellProductionConstraintsMet =    checkCellProductionConstraints(pdStorage,opt);
	bool fecundityScheduleConstraintsMet = checkFecundityConstraints(pdStorage,opt);
	bool allConstraintsMet = cellProductionConstraintsMet && fecundityScheduleConstraintsMet;
	if (opt->fitnessMetric.find("tedNIPS:")!=std::string::npos) {
		int t11, t13,v11,v13;
		float t12,v12;
		Array1D<int> t2, v2, geometryDiff;
		bool constraint1, constraint2;
		tedConstraints(t11, t12, t13, &t2, v11, v12, v13, &v2, &geometryDiff, constraint1, constraint2, opt);
		if (constraint1 && constraint2) {
			allConstraintsMet=true;
		} else {
			allConstraintsMet=false;
		}
	}


	// display results
	if (opt->verbose) {
		printSimulationResults(pdStorage, opt->numSperm, cellProductionConstraintsMet, fecundityScheduleConstraintsMet, allConstraintsMet,opt);
	}
	if (opt->verbose_tedCellProduction) {
		printSimulationResultsForTed_cellProduction(pdStorage,opt);
	}
	if (opt->verbose_tedPedigreeDepth) {
		printSimulationResultsForTed_pedigreeDepth(pdStorage, cellProductionConstraintsMet, fecundityScheduleConstraintsMet, allConstraintsMet, opt);
	}

	// return weighted pedigree depth if constraints satisfied, otherwise return -1
	float rn=0;
	if (allConstraintsMet) {rn=resultStats[0];}
	else {rn=-1;}
	if (opt->verbose) {printf("Return value = %f\n",rn);}

	// clean up
	delete opt; opt=NULL;
	delete pdStorage; pdStorage=NULL;
	// return
	return rn;
}
