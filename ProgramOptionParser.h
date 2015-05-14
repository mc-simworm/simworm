/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#ifndef PROGRAMOPTIONPARSER_H_
#define PROGRAMOPTIONPARSER_H_

#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include "CellCycleProfile.h"
#include "GeometryProfile.h"
#include "PhaseIndexProfile.h"
#include "MeioticFractionProfile.h"
namespace PO = boost::program_options;

/*!
 * \class ProgramOptionParser
 * \brief This class provides functionality to parse user-defined program options
 */
class ProgramOptionParser {
public:
	std::string pedigreeDepthOrCellCycleFit;	// Parse program options for pedigree depth simulations (pd) or cell cycle fits (ccfit)
	std::string cellCycleProfile;			// Serialized cell cycle profile
	std::string mitoticRegionGeometry;		// Serialized mitotic region geometry
	std::string meioticRegionGeometry;		// Serialized meiotic region geometry
	std::string truncatedMitoticIndex;		// fraction of cells in mitotic region that are pre-meiotic
	bool debug;						    // Write debugging messages to the standard output
	bool verbose;						// Write simulation results to the standard output
	bool verbose_tedCellProduction;					// Write simulation results to the standard output in a form parsable by Ted's optimizer
	bool verbose_tedPedigreeDepth;
	int numSim;							// Number of germline simulations to run
	float noise;						// Noise in cell cycle length for individual cells.
	int randomSeed;						// Random number seed.  -1 to choose a random seed.
	float deleteriousMutationRate;      // Number of deleterious mutations per genome per generation.  Used to calculate fitness
	float maximalGenerationRate;		// Maximal generation rate in hours^-1.  Used to calculate fitness
	float maximalPedigreeDepth;
	int pedigreeDepthIncrementMethod;	// Method of incrementing pedigree depth counter

	// parameters specific to pedigree depth simulations
	std::string cellProductionConstraints;
	std::string oocyteProbability;
	std::string fecunditySchedule;
	int numSperm;
	int numOocyte;
	int numDivisions;
	std::string fitnessMetric;
	std::string constraintMethod;
	float forwardMovementProbability;
	int maxMitoticRegionSize;
	bool filledGermline;

	// parameters specific to ccfit
	float minTime;
	std::string minTimeStr;
	std::string mitoticPhaseIndices;
	std::string dnaBin;
	std::string fitData;
	float dnaNoise;

	ProgramOptionParser(std::string mode);
	virtual ~ProgramOptionParser();
	void parseArgs (int argc, char *argv[]);
	void display();
	void getCellCycleProfile(CellCycleProfile *cc);
	void getGeometryProfile(GeometryProfile *gg);
	bool isPedigreeDepthSimulation();
	int getNumCellRows();
	void setMinTime(int i);
	void getDnaBin(Array1D<float> *dnaBins);
	void getChaseTimes(Array1D<float> *chaseTimes);
	int getNumChaseTimes();
	void getPhaseIndicesProfile(PhaseIndexProfile *pp);
	void getMeioticFractionProfile(MeioticFractionProfile *mfp);
	void getCellProductionConstraints(Array2D<float> *cellProductionConstraintsMatrix);
	void getMeioticRegionSize(Array2D<int> *meioticRegionSz);
	void getOocyteProbability(Array2D<float> *probability);
	std::string oocyteProbabilityInterpolationMethod();
	void getFecundityConstraints(Array2D<float> *pedigreeDepthWeights);
};
#endif /* PROGRAMOPTIONPARSER_H_ */
