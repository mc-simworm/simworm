/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#ifndef GERMCELL_H_
#define GERMCELL_H_

#include "util.h"
#include "Array2D.h"
#include "PhaseIndexProfile.h"
#include "CellCycleProfile.h"
#include <boost/math/special_functions/round.hpp>
#include "MeioticFractionProfile.h"

/*!
 * This class records cellular data for a single cell in the simulation
 */
class GermCell {
private:
	int idx_;								// number used to identify cell in mapPosition_
	int x_;									// cell row position in mitotic region.
	int y_;									// y-position (i.e., position in circular ring) in the mitotic region
	float pedigreeDepth_;						// number of cell divisions removed from founder cell
	float age_;								// age of cell
	float meanCellCycleLength_;				// mean cell cycle length of cells in x_
	float cellCycleNoise_;					// % noise in cell cycle
	float cellCycleLength_;					// actual length of cell cycle for cell
	PhaseIndexProfile *phaseIndexProfile_;  // spatial profile of cell cycle phase indices
	float dnaNoise_;						// amount of noise in DNA content measurements
	std::vector<std::vector<int>*> *mapPosition_;		// communication interface between germ cells
	std::vector<GermCell *> *germCellVector_;			// communication interface between germ cells
	CellCycleProfile *cellCycleProfile_;				// communication interface between germ cells
	int *cellDivisionsElapsed_;							// communication interfact between germ cells
	MeioticFractionProfile *meioticFractionProfile_;	// probability of cell entering pre-meiotic G2 instead of dividing
	int mitoticOrMeioticOrLaid_; 			// 0 if germ cell is in mitotic region, 1 if germ cell is in meiotic region, 2 if germ cell has undergone gametogenesis, 3 if germ cell is dead
	bool isOocyte_;							// germ cell has undergone oogenesis
	bool isSpermatocyte_;					// germ cell has undergone spermatogenesis
	float timeGametogenesis_hours;			// time (in hours) when the germ cell underwent gametogenesis.  -1 if the cell never underwent gametogenesis
	float timeGametogenesis_cellDivisions; 	// time (in cell divisions) when the germ cell underwent gametogenesis.  -1 if the cell never underwent gametogenesis
	int pedigreeDepthIncrementMethod_;		// 0 if simple increment, 1 if increment by 1/ccLength, 2 if increment by 1/(ccLength/minCcLength)
	base_generator_type *generator_;
public:
	int preMeiotic_;						// -1 if undecided, 0 if mitotic, 1 if preMeiotic
	bool edu_;								// EdU content of germ cell
	float timeLeftMitoticRegion_;			// time (in hours) that the germ cell left the mitotic region.  -1 if cell never left mitotic region
	int meioticEntryOrder_; // time (in cell divisions) that the germ cell left the mitotic region.  -1 if the cell never left the mitotic region

	GermCell(int idx, int x, int y, float percentAge, float meanCellCycleLength, float cellCycleNoise, float pedigreeDepth, std::vector<std::vector<int> *> *mapPosition, std::vector<GermCell *> *germCellVector, bool edu, PhaseIndexProfile *phaseIndexProfile, float dnaNoise, int pedigreeDepthIncrementMethod, CellCycleProfile *cellCycleProfile, int *cellDivisionsElapsed, MeioticFractionProfile *mfp,base_generator_type *generator);
	virtual ~GermCell();
	float getTimeToNextDivision();
	void incrementAge(float dt);
	int mitoticOrMeioticOrLaid();
	float getCellCycleLength();
	int getX();
	int getY();
	float getAge();
	void setNoisyCellCycleLength();
	int divide();
	float getPercentAge();
	bool edu();
	bool isGamete();
	bool isSperm();
	bool isOocyte();
	void die();
	bool isMphase();
	int getCellCyclePhase();
	float getDna();
	void setEdU(bool EdU);
	float getPedigreeDepth();
	int setPosition(int x, int y, float meanCellCycleLengthAtX);
	void rescaleCellCycleLength(float newMeanCellCycleLength);
	void pushCellFromMitoticRegion(float time_left_MR,int time_left_MR_cellDivisions);
	void spermatogenesis(float time_left_meioticRegion, float time_left_meioticRegion_cellDivision);
	void oogenesis(float time_left_meioticRegion, float time_left_meioticRegion_cellDivision);
	float getTimeLeftMitoticRegion();
	int getTimeLeftMitoticRegion_cellDivisions();
	float getTimeGametogenesis_hours();
	float getTimeGametogenesis_cellDivisions();
	float getTimeToNextG2toM();
};

#endif /* GERMCELL_H_ */
