/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/
#include "GermCell.h"
using namespace std;

/*!
 * Constructor
 */
GermCell::GermCell(int idx, int x, int y, float percentAge, float meanCellCycleLength, float cellCycleNoise, float pedigreeDepth, std::vector<std::vector<int> *> *mapPosition, std::vector<GermCell *> *germCellVector, bool edu, PhaseIndexProfile *phaseIndexProfile, float dnaNoise, int pedigreeDepthIncrementMethod, CellCycleProfile *cellCycleProfile, int *cellDivisionsElapsed, MeioticFractionProfile *mfp,base_generator_type *generator) {
	generator_ = generator;
	idx_ = idx;
	x_=x;
	y_=y;
	pedigreeDepth_=pedigreeDepth;
	meanCellCycleLength_=meanCellCycleLength;
	cellCycleNoise_=cellCycleNoise;
	mapPosition_=mapPosition;
	edu_ = edu;
	phaseIndexProfile_ = phaseIndexProfile;
	mitoticOrMeioticOrLaid_=0;
	isOocyte_=false;
	isSpermatocyte_=false;
	timeLeftMitoticRegion_=-1;
	meioticEntryOrder_=-1;
	timeGametogenesis_hours=-1;
	timeGametogenesis_cellDivisions = -1;
	dnaNoise_ = dnaNoise;
	setNoisyCellCycleLength();
	age_=percentAge*cellCycleLength_;
	if ((*(*mapPosition_)[(unsigned int)x])[(unsigned int)y]==-1) {
		(*(*mapPosition_)[(unsigned int)x])[(unsigned int)y]=idx_;
	}
	germCellVector_=germCellVector;
	germCellVector_->push_back(this);
	pedigreeDepthIncrementMethod_=pedigreeDepthIncrementMethod;
	cellCycleProfile_=cellCycleProfile;
	cellDivisionsElapsed_=cellDivisionsElapsed;
	preMeiotic_=-1;
	meioticFractionProfile_ = mfp;
}

/*!
 * Destructor
 */
GermCell::~GermCell() {

}

/*!
 * Return time (hours) remaining until cell divides
 */
float GermCell::getTimeToNextDivision() {
	return cellCycleLength_-age_;
}

/*!
 * Return time (hours) remaining until cell enters G2. Return can be negative if cell is in M-phase
 */
float GermCell::getTimeToNextG2toM() {
	float percentAgeTransition = phaseIndexProfile_->getPercentAgeG2toM(x_,*cellDivisionsElapsed_);
	float ageTransition = cellCycleLength_*percentAgeTransition;
	float timeToG2M = ageTransition-age_;
	return timeToG2M;
}

/*!
 * Increment age of cell by dt
 */
void GermCell::incrementAge(float dt) {
	age_ += dt;
	if ((age_ > cellCycleLength_) & (preMeiotic_==1)) {  // pre-meiotic cells can age without dividing, reset age in this case
		age_ = cellCycleLength_;
	}
}

/*!
 * Return 0 if germ cell is in mitotic region, 1 if germ cell is in meiotic region, 2 if germ cell has left meiotic region
 */
int GermCell::mitoticOrMeioticOrLaid() {
	return mitoticOrMeioticOrLaid_;
}

/*!
 * Set cell cycle length with a bit of noise
 */
void GermCell::setNoisyCellCycleLength() {
	float lower =std::max(meanCellCycleLength_-cellCycleNoise_*meanCellCycleLength_,(float)0);
	float upper =meanCellCycleLength_+cellCycleNoise_*meanCellCycleLength_;
	cellCycleLength_= randUniformFloat(lower,upper,generator_);
}

/*!
 * Set spatial position (x,y position in mitotic region) of cell
 * Also rescale cell cycle length to account for spatial profile of cell cycle lengths along x-axis
 */
int GermCell::setPosition(int x, int y, float meanCellCycleLengthAtX) {
	int cellBeingPushed =(*(*mapPosition_)[(unsigned int)x])[(unsigned int)y];
	x_=x;
	y_=y;
	(*(*mapPosition_)[(unsigned int)x])[(unsigned int)y]=idx_;
	rescaleCellCycleLength(meanCellCycleLengthAtX);
	return cellBeingPushed;
}

/*!
 * Make the cell divide.
 * Returns 1 if cell did not divide due to being set to pre-meiotic
 */
int GermCell::divide() {
	// Sanity check: pre-meiotic cells should not divide
	if (preMeiotic_==1) {
		cerr<<"Pre-meiotic cells cannot divide"<<endl;
		exit(1);
	}

	// otherwise, if cell is not pre-meiotic, increment pedigree depth
	float daughterPedigreeDepth;
	if (pedigreeDepthIncrementMethod_==0) {
		pedigreeDepth_++;
		daughterPedigreeDepth=pedigreeDepth_;
	} else if (pedigreeDepthIncrementMethod_==1) {
		pedigreeDepth_+=1/meanCellCycleLength_;
		daughterPedigreeDepth=pedigreeDepth_;
	} else if (pedigreeDepthIncrementMethod_==2) {
		float minCellCycleLength = cellCycleProfile_->getCurrentProfileMinimum(*cellDivisionsElapsed_);
		pedigreeDepth_+=1/(meanCellCycleLength_/minCellCycleLength);
		daughterPedigreeDepth=pedigreeDepth_;
	} else if (pedigreeDepthIncrementMethod_==3) {	// golden strand
		if (x_==0) {
			daughterPedigreeDepth=pedigreeDepth_+2;
		} else {
			pedigreeDepth_++;
			daughterPedigreeDepth=pedigreeDepth_;
		}
	} else {
		cerr<<"Invalid option for pedigreeDepthIncrementMethod_ (should be 0,1,2)" <<endl;
		exit(1);
	}

	// get new cell cycle length
	setNoisyCellCycleLength();

	// reset age to zero
	age_=0;

	// reset premeiotic state
	preMeiotic_=-1;

	// create daughter cell
	int daughterIdx=(int)germCellVector_->size();
	new GermCell(daughterIdx,x_,y_,age_,meanCellCycleLength_,cellCycleNoise_,daughterPedigreeDepth,mapPosition_,germCellVector_,edu_,phaseIndexProfile_,dnaNoise_,pedigreeDepthIncrementMethod_,cellCycleProfile_,cellDivisionsElapsed_,meioticFractionProfile_,generator_);
	return daughterIdx;
}

/*!
 * Rescale cell cycle length
 */
void GermCell::rescaleCellCycleLength(float newMeanCellCycleLength) {
	if (abs(meanCellCycleLength_-newMeanCellCycleLength)>EPS) {
		float rescaleFactor = newMeanCellCycleLength/meanCellCycleLength_;
		meanCellCycleLength_ = newMeanCellCycleLength;
		cellCycleLength_ = rescaleFactor*cellCycleLength_;
		age_ = rescaleFactor*age_;
	}
}

/*!
 * Push a cell from the mitotic region into the meiotic region
 */
void GermCell::pushCellFromMitoticRegion(float time_left_MR,int meioticEntryOrder){
	// Sanity check
	if (mitoticOrMeioticOrLaid_!=0) {
		cerr<<"Error: Pushing cell from mitotic region when cell is not in mitotic region"<<endl;
		exit(1);
	}
	mitoticOrMeioticOrLaid_=1;
	timeLeftMitoticRegion_=time_left_MR;
	meioticEntryOrder_=meioticEntryOrder;
}


/*!
 * Undergo spermatogenesis
 */
void GermCell::spermatogenesis(float timeSpermatogenesis_hours, float timeSpermatogenesis_cellDivisions){
	// Sanity check
	if (mitoticOrMeioticOrLaid_!=1) {
		cerr<<"Error: Pushing cell from meiotic region when cell is not in meiotic region"<<endl;
		exit(1);
	}
	mitoticOrMeioticOrLaid_=2;
	timeGametogenesis_hours=timeSpermatogenesis_hours;
	timeGametogenesis_cellDivisions = timeSpermatogenesis_cellDivisions;
	isOocyte_ = false;
	isSpermatocyte_ = true;
}

/*!
 *  Undergo oogenesis
 */
void GermCell::oogenesis(float timeOogenesis_hours, float timeOogenesis_cellDivisions){
	// Sanity check
	if (mitoticOrMeioticOrLaid_!=1) {
		cerr<<"Error: Pushing cell from meiotic region when cell is not in meiotic region"<<endl;
		exit(1);
	}
	mitoticOrMeioticOrLaid_=2;
	timeGametogenesis_hours=timeOogenesis_hours;
	timeGametogenesis_cellDivisions = timeOogenesis_cellDivisions;
	isOocyte_ = true;
	isSpermatocyte_ = false;
}

/*!
 * Return cell cycle length (hours)
 */
float GermCell::getCellCycleLength() {
	return cellCycleLength_;
}

/*!
 * Return x-position of cell (distal-proximal axis)
 */
int GermCell::getX() {
	return x_;
}

/*!
 * Return y-position of cell (location in cell row)
 */
int GermCell::getY() {
	return y_;
}

/*!
 * Return the age of the cell (hours)
 */
float GermCell::getAge() {
	return age_;
}

/*!
 * Return percent completion of the cell cycle, exists in [0,1]
 */
float GermCell::getPercentAge() {
	return age_/cellCycleLength_;
}

/*!
 * Return EdU content of cell
 */
bool GermCell::edu() {
	return edu_;
}

/*!
 * Return true if the cell is a gamete (oocyte or spermatocyte)
 */
bool GermCell::isGamete() {
	if (isOocyte_ || isSpermatocyte_) {
		return true;
	}
	return false;
}

/*!
 * Return true if the cell is sperm
 */
bool GermCell::isSperm() {
	if (isSpermatocyte_) {
		return true;
	}
	return false;
}

/*!
 * Return true if the cell is oocyte
 */
bool GermCell::isOocyte() {
	if (isOocyte_) {
		return true;
	}
	return false;
}

/*!
 * Undergo apoptosis
 */
void GermCell::die() {
	isOocyte_=false;
	isSpermatocyte_=false;
	mitoticOrMeioticOrLaid_=2;
}

/*!
 * Return true if cell is in M-phase
 */
bool GermCell::isMphase() {
	return phaseIndexProfile_->isMphase(age_/cellCycleLength_,x_,*cellDivisionsElapsed_);
}

/*!
 * Return cell cycle phase (1 for G1, 2 for S, 3 for G2, 4 for M)
 */
int GermCell::getCellCyclePhase() {
	int rn =  phaseIndexProfile_->getPhase(age_/cellCycleLength_,x_,*cellDivisionsElapsed_);
	if (preMeiotic_==1) { // pre-meiotic cells must be in G2
		rn = 3;
	}
	return rn;
}

/*!
 * Return DNA content of cell
 */
float GermCell::getDna() {
	float dnaContent = phaseIndexProfile_->getDna(age_/cellCycleLength_,x_,*cellDivisionsElapsed_);
	float dnaContent_withNoise = dnaContent + randNormal(0,dnaNoise_,generator_);
	return dnaContent_withNoise;
}

/*!
 * Set EdU content of the cell
 */
void GermCell::setEdU(bool EdU) {
	edu_ = EdU;
}

/*!
 * Return pedigree depth of the cell
 */
float GermCell::getPedigreeDepth() {
	return pedigreeDepth_;
}

/*!
 * Return time (hours) cell left the mitotic region
 */
float GermCell::getTimeLeftMitoticRegion() {
	return timeLeftMitoticRegion_;
}

/*!
 * Return time (cell divisions) cell left the mitotic region
 */
int GermCell::getTimeLeftMitoticRegion_cellDivisions() {
	return meioticEntryOrder_;
}

/*!
 * Return time (hours) of gametogenesis
 */
float GermCell::getTimeGametogenesis_hours() {
	if (mitoticOrMeioticOrLaid_!=2 || (!isSpermatocyte_ & !isOocyte_)) {
		cerr<<"Error: Cell is not a gamete"<<endl;
		exit(1);
	}
	return timeGametogenesis_hours;
}

/*!
 * Return time (cell divisions) cell left the meiotic region
 */
float GermCell::getTimeGametogenesis_cellDivisions() {
	if (mitoticOrMeioticOrLaid_!=2 || (!isSpermatocyte_ & !isOocyte_)) {
		cerr<<"Error: Cell is not a gamete"<<endl;
		exit(1);
	}
	return timeGametogenesis_cellDivisions;
}
