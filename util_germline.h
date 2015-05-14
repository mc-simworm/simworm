/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#include <vector>
#include "Array1D.h"
#include "Array2D.h"
#include "util.h"
#include "ProgramOptionParser.h"
#include <Dense>
#include <Eigenvalues>
#include "GermLine.h"

#ifndef PD_STORAGE_
#define PD_STORAGE_
struct PedigreeDepthStorage {
	Array1D<float> *pd_sperm;			// average pedigree depth of sperm
	Array1D<float> *pd_oocyte;			// average pedigree depth of oocytes
	Array1D<float> *pd_spermOocyte;	// average pedigree depth of progeny (i.e., sperm+oocyte/2)
	Array1D<float> *pd_mitoticRegion;  // average pedigree depth of the mitotic region at the end of the simulation
	Array1D<float> *pd_meioticRegion;  // average pedigree depth of the meiotic region at the end of the simulation
	Array1D<float> *pd_wholeGermline;  // average pedigree depth of the entire germline (mitotic+meiotic region) at the end of the simulation
	Array1D<float> *pd_spermOocyteSum; // average pedigree depth of all gametes.  This is pd_spermOocyte generalized to the mated case (i.e., number of sperm != number of oocytes)
	Array1D<float> *pd_GermlineAtEggLaying;    // average pedigree depth at time of egg laying events
	Array1D<float> *pd_AllCellsEnterMeioticRegion;    // average pedigree depth at time of egg laying events
	Array1D<float> *num_meioticCellsProduced;  // total number of cells that enter the meiotic region
	Array1D<float> *num_cellDivisions;         // total number of cell divisions by the end of the simulation
	Array1D<float> *time_firstMeioticExit;		// time (hours) when the first cell exits the meiotic region and either undergoes oogenesis/apoptosi
	Array1D<float> *time_nthMitoticExit;		// time (hours) when the nth cell exists the meiotic region as specified by --fitnessMetric meanPedigreeDepthMeioticCells:N

	Array2D<float> *cellProductionConstraints;
	Array2D<float> *simulatedCellProductionSchedule;
	Array2D<float> *fecundityConstraints;
	Array2D<float> *simulatedFecunditySchedule;
	Array2D<float> *simulatedCellProductionSchedule_atFecundityBins;
	Array2D<float> *simulatedMeioticExitSchedule;

	PedigreeDepthStorage(ProgramOptionParser *opt) {
		pd_sperm = new Array1D<float>;
		pd_oocyte = new Array1D<float>;
		pd_spermOocyte = new Array1D<float>;
		pd_mitoticRegion = new Array1D<float>;
		pd_meioticRegion = new Array1D<float>;
		pd_wholeGermline = new Array1D<float>;
		pd_spermOocyteSum = new Array1D<float>;
		pd_GermlineAtEggLaying = new Array1D<float>;
		pd_AllCellsEnterMeioticRegion = new Array1D<float>;
		num_meioticCellsProduced = new Array1D<float>;
		num_cellDivisions = new Array1D<float>;
		time_firstMeioticExit = new Array1D<float>;
		time_nthMitoticExit = new Array1D<float>;
		cellProductionConstraints = new Array2D<float>;

		pd_sperm->resize(opt->numSim);
		pd_oocyte->resize(opt->numSim);
		pd_spermOocyte->resize(opt->numSim);
		pd_mitoticRegion->resize(opt->numSim);
		pd_meioticRegion->resize(opt->numSim);
		pd_wholeGermline->resize(opt->numSim);
		pd_spermOocyteSum->resize(opt->numSim);
		pd_GermlineAtEggLaying->resize(opt->numSim);
		pd_AllCellsEnterMeioticRegion->resize(opt->numSim);
		num_meioticCellsProduced->resize(opt->numSim);
		num_cellDivisions->resize(opt->numSim);
		time_firstMeioticExit->resize(opt->numSim);
		time_nthMitoticExit->resize(opt->numSim);

		opt->getCellProductionConstraints(cellProductionConstraints);
		Array2D<float>::array_index dimx,dimy;
		cellProductionConstraints->getDimensions(dimx,dimy);
		simulatedCellProductionSchedule = new Array2D<float>;
		simulatedCellProductionSchedule->resize(opt->numSim,dimy);
		fecundityConstraints = new Array2D<float>;
		opt->getFecundityConstraints(fecundityConstraints);
		fecundityConstraints->getDimensions(dimx,dimy);
		simulatedFecunditySchedule = new Array2D<float>;
		simulatedFecunditySchedule->resize(opt->numSim,dimy);
		simulatedCellProductionSchedule_atFecundityBins = new Array2D<float>;
		simulatedCellProductionSchedule_atFecundityBins->resize(opt->numSim,dimy);
		simulatedMeioticExitSchedule = new Array2D<float>;
		simulatedMeioticExitSchedule->resize(opt->numSim,dimy);
    }

	~PedigreeDepthStorage() {
		delete pd_sperm; pd_sperm=NULL;
		delete pd_oocyte; pd_oocyte=NULL;
		delete pd_spermOocyte; pd_spermOocyte=NULL;
		delete pd_mitoticRegion; pd_mitoticRegion=NULL;
		delete pd_meioticRegion; pd_meioticRegion=NULL;
		delete pd_wholeGermline; pd_wholeGermline=NULL;
		delete pd_spermOocyteSum; pd_spermOocyteSum=NULL;
		delete pd_GermlineAtEggLaying; pd_GermlineAtEggLaying=NULL;
		delete pd_AllCellsEnterMeioticRegion; pd_AllCellsEnterMeioticRegion=NULL;
		delete num_meioticCellsProduced; num_meioticCellsProduced=NULL;
		delete num_cellDivisions; num_cellDivisions=NULL;
		delete cellProductionConstraints; cellProductionConstraints=NULL;
		delete simulatedCellProductionSchedule;  simulatedCellProductionSchedule=NULL;
		delete fecundityConstraints; fecundityConstraints=NULL;
		delete simulatedFecunditySchedule; simulatedFecunditySchedule=NULL;
		delete simulatedCellProductionSchedule_atFecundityBins; simulatedCellProductionSchedule_atFecundityBins=NULL;
		delete simulatedMeioticExitSchedule; simulatedMeioticExitSchedule=NULL;
		delete time_firstMeioticExit; time_firstMeioticExit=NULL;
		delete time_nthMitoticExit; time_nthMitoticExit=NULL;
	}
};
#endif

bool mitoticRegionTooBig(ProgramOptionParser *opt);
void sort2(Array1D<float> *v1, Array1D<float> *v2);
void printSimulationResults(PedigreeDepthStorage *pdStorage, int numSperm, bool cellProductionConstraintsMet, bool fecundityScheduleConstraintsMet, bool allConstraintsMet, ProgramOptionParser *opt);
void printSimulationResultsForTed_cellProduction(PedigreeDepthStorage *pdStorage, ProgramOptionParser *opt);
bool checkCellProductionConstraints(PedigreeDepthStorage *pdStorage, ProgramOptionParser *opt);
bool checkFecundityConstraints(PedigreeDepthStorage *pdStorage, ProgramOptionParser *opt);
float calculateFitness(const char *option, PedigreeDepthStorage *pdStorage, ProgramOptionParser *opt);
int getClosestIdx(Array1D<float> *vec, float val);
float calculateGenerationRate(Array2D<float> *simulatedFecunditySchedule,Array2D<float> *fecundityConstraints);
float calculateGenerationFitness(Array2D<float> *simulatedFecunditySchedule,Array2D<float> *fecundityConstraints, float pd, float deleteriousMutationRate, float maximalGenerationRate, float maximalPedigreeDepth);
float calculateGenerationFitness_selectionCoefficient(Array2D<float> *simulatedFecunditySchedule,Array2D<float> *fecundityConstraints, float pd, float deleteriousMutationRate, float maximalGenerationRate, float maximalPedigreeDepth);
float calculateDoublingTime(float r);
void printCellPhaseIndices(Array2D<float> *G1idx, Array2D<float> *Sidx, Array2D<float> *G2idx, Array2D<float> *Midx);
void tedConstraints(int &t11, float &t12, int &t13, Array1D<int> *t2, int &v11, float &v12, int &v13, Array1D<int> *v2, Array1D<int> *geometryDiff, bool &constraint1, bool &constraint2, ProgramOptionParser *opt);
void printSimulationResultsForTed_pedigreeDepth(PedigreeDepthStorage *pdStorage, bool cellProductionConstraintsMet, bool fecundityScheduleConstraintsMet, bool allConstraintsMet, ProgramOptionParser *opt);
void printSimulationResultsForTed_pedigreeDepth_badMR(PedigreeDepthStorage *pdStorage, ProgramOptionParser *opt);
