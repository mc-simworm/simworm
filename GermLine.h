/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#ifndef GERMLINE_H_
#define GERMLINE_H_

#include "util.h"
#include <stdio.h>
#include "GeometryProfile.h"
#include "CellCycleProfile.h"
#include "Array2D.h"
#include "GermCell.h"
#include "ProgramOptionParser.h"
#include "PhaseIndexProfile.h"
#include "util_germline.h"
#include "MeioticFractionProfile.h"
#include <algorithm>

/*!
 * This class implements a germline simulation
 */
class GermLine {
private:

	// Below are schedules that determine germline properties over time
	GeometryProfile *geometry_;				// geometry profile of the mitotic region over time
	CellCycleProfile *cellCycle_;			// cell cycle profile in the mitotic region over time
	PhaseIndexProfile *phaseIndexProfile_;  // phase indices profile in the mitotic region over time
	MeioticFractionProfile *meioticFractionProfile_; // fraction of pre-meiotic cells in the mitotic region over time
	Array2D<int> *meioticRegionSize_;		// geometry (i.e., number of cells) in the meiotic region over time
	Array2D<float> *oocyteProbability_;		// probability of oogenesis (vs. apoptosis) over time

	// position of mitotic region cells
	std::vector<std::vector<int>*> *mapPosition_;	// _mapPosition[x][y] gives the index of the GermCell at position (x,y).  -1 if there is no cell

	// Simulation parameters
	float dnaNoise_;
	float cellCycleNoise_;
	bool verbose_;
	float movementProbabilityForward_;
	std::string oocyteProbabilityMethod_;
	int pedigreeDepthIncrementMethod_;

	// timers
	float tElapse_;
	int cellDivisionsElapse_;
	int countSperm_;
	int countOocyteLaid_;
	int cellsPushedIntoMR_;
	int cellsExitMeioticRegion_;
	int countCellsEnterMeioticRegion_;

	// stopping conditions
	int minNumOocytesProduced_;
	float minTimeElapsed_;
	int minCellDivisions_;
	int numSpermToProduce_;
	int minCellsEnterMeioticRegion_;

	// storage vectors for simulation results
	std::vector<int> *meioticGermCellVector_;		// vector storing indices of GermCells in the meiotic region
	std::vector<float> *cellProductionTimes_;  	    // (*cellProductionTimes_)[i-1] gives the time (hours) of the ith cell division
	std::vector<int> *exitMeioticRegionTime_;	    // (*exitMeioticRegionTime_)[i-1] gives the time (hours) of the ith cell exiting the meiotic region (i.e., undergoing oogenesis/apoptosis)
	std::vector<float> *pdGermlineAtEggLaying;		// (*pdGermlineAtEggLaying)[i-1] gives the average pedigree depth of entire germline at the ith egg laying event

	// random number generator
	base_generator_type *generator;

	// record time that first cell exits meiotic region

public:

	std::vector<GermCell *> *germCellVector_;		// vector storing every GermCell in the simulation

	GermLine(ProgramOptionParser *args, int seed);
	virtual ~GermLine();
	void initialize_oneCell();
	void initialize_filledGermline();
	void run(std::string initialCondition, int preRunSteps);
	void calculateNextEvent(unsigned int &cellIndex, float &timeToNextEvent, int &event);
	void incrementAge(float dt);
	void display();
	void pushCellFromMR(unsigned int index);
	void moveCells(unsigned int pushingCellIdx, int direction_UpDown, int direction_Forward);
	void restructureGeometry();
	void rescaleCellCycleLengths();
	void reproductiveDiapause();
	void calculateOutput(float &pdAllCellsEnterMeioticRegion, float &pdSperm, float &pdOocyte, float &pdSpermAndOocyte, float &pdSpermAndOocyteSum, float &pdMitoticRegion, float &pdMeioticRegion, float &pdWholeGermline, float &pdGermlineAtEggLaying, int &numCellsEnterMeioticRegion, int &numCellDivisions, float &timeFirstCellExitMeioticRegion, float &timeNthCellExitsMitoticRegion);
	void calculateRawOutput(Array1D<float> *allCellsEnterMeioticRegionPedigreeDepth, Array1D<float> *oocyteProductionTimes, Array1D<float> *spermProductionTimes, Array1D<float> *oocytePedigreeDepth, Array1D<float> *spermPedigreeDepth, Array1D<float> *oocyteProductionTimes_cellDivision, Array1D<float> *spermProductionTimes_cellDivision, Array1D<float> *mitoticRegionPedigreeDepths, Array1D<float> *meioticRegionPedigreeDepths, Array1D<float> *wholeGermlinePedigreeDepths, Array1D<float> *cellProductionTimes, Array1D<int> *cellExitTimes, Array1D<float> *pedigreeDepthAtEggLaying);
	void calculateSchedules(int i, PedigreeDepthStorage *pdStorage);
	void setPremeioticCells();
	void displayPreMeioticState();
	void displayPreMeioticFraction();
	void getSpatialCellCyclePhaseIndices(Array1D<float> *G1idx,Array1D<float> *Sidx,Array1D<float> *G2idx,Array1D<float> *Midx);
	float getFractionMphase(int x);
	void displayPhase();
	float calculatePercentPremeiotic();
	int step(float &time_nextEvent);
	void pulseEdU();
	void moveCellsWrapper(unsigned int pushingCellIdx);
	void displayPedigreeDepth();
	//void forcePremeioticEntry();
	//int calculateOutput_numExitMeioticRegion(float t1, float t2);
	//int calculateOutput_numDivisions(float t1, float t2);
};

#endif /* GERMLINE_H_ */
