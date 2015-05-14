/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#include "PhaseIndexProfile.h"
using namespace std;
using namespace boost;

/*!
 * Constructor
 */
PhaseIndexProfile::PhaseIndexProfile() {
	G1_ = new Array2D<float>();
	S_  = new Array2D<float>();
	G2_ = new Array2D<float>();
	M_  = new Array2D<float>();
}

/*!
 * Destructor
 */
PhaseIndexProfile::~PhaseIndexProfile() {
	delete G1_; G1_=NULL;
	delete S_;  S_=NULL;
	delete G2_; G2_=NULL;
	delete M_;  M_=NULL;
}

/*!
 * Return the mitotic index in cell row x at time d (cell divisions)
 */
float PhaseIndexProfile::getMphaseIndex(int x, int d) {
	// return (*M_)(x);
	return M_->calculateInterpolatedValue(x,d);
}

/*!
 * Write cell cycle phase indices to standard output
 */
void PhaseIndexProfile::display() {
	cout << "G1 indices: " << endl;
	G1_->display();
	cout << "S indices: " << endl;
	S_->display();
	cout << "G2 indices: " << endl;
	G2_->display();
	cout << "M indices: " << endl;
	M_->display();
}

/*!
 * Read formatted string
 */
void PhaseIndexProfile::parseFormattedParameters(std::string serializedParameters, std::string mitoticRegionGeometry) {

	// input string is of format 0.1*0.6*0.25*0.05
	if (serializedParameters.find(":")==std::string::npos) {
		GeometryProfile gg;
		gg.parseFormattedParameters(mitoticRegionGeometry);
		int dimx = gg.numCellRows();
		G1_->resize(dimx,1);
		S_->resize(dimx,1);
		G2_->resize(dimx,1);
		M_->resize(dimx,1);
		vector<string> splitStr;
		boost::split(splitStr,serializedParameters,boost::is_any_of("*"));
		G1_->fill(boost::lexical_cast<float>(splitStr[0].c_str()));
		S_->fill(boost::lexical_cast<float>(splitStr[1].c_str()));
		G2_->fill(boost::lexical_cast<float>(splitStr[2].c_str()));
		M_->fill(boost::lexical_cast<float>(splitStr[3].c_str()));
	} else {
		vector<string> splitStr;
		boost::split(splitStr,serializedParameters,boost::is_any_of("*"));
		G1_->parseString(splitStr[0]);
		S_->parseString(splitStr[1]);
		G2_->parseString(splitStr[2]);
		M_->parseString(splitStr[3]);

		// make sure phase indices are normalized
		int dimx,dimy;
		G1_->getDimensions(dimx,dimy);
		for (int x=0;x<dimx;x++) {
			for (int y=0;y<dimy;y++) {
				float G1 = (*G1_)(x,y);
				float S  = (*S_)(x,y);
				float G2 = (*G2_)(x,y);
				float M  = (*M_)(x,y);
				float sum = G1+S+G2+M;
				(*G1_)(x,y) = G1/sum;
				(*S_)(x,y) = S/sum;
				(*G2_)(x,y) = G2/sum;
				(*M_)(x,y) = M/sum;
			}
		}
	}
}

/*!
 * return 1 if cell is in G1, 2 if S, 3 if G2, 4 if M
 * d is time in cell divisions
 */
int PhaseIndexProfile::getPhase(float percentAge, int x, int d) {
	// calculate transition times from phase indices assuming expontential age distribution
	/*
	float G1 = (*G1_)(x);
	float S  = (*S_)(x);
	float G2 = (*G2_)(x);
	*/
	float G1 = G1_->calculateInterpolatedValue(x,d);
	float S  = S_->calculateInterpolatedValue(x,d);
	float G2 = G2_->calculateInterpolatedValue(x,d);
	float transition_G1toS = (logf(2)-logf(2-G1))/logf(2);
	float transition_StoG2 = (logf(2)-logf(2-(G1+S)))/logf(2);
	float transition_G2toM = (logf(2)-logf(2-(G1+S+G2)))/logf(2);

	if (percentAge<transition_G1toS) {
		return 1;
	} else if (percentAge>=transition_G1toS && percentAge<transition_StoG2) {
		return 2;
	} else if (percentAge>=transition_StoG2 && percentAge<transition_G2toM) {
		return 3;
	} else {
		return 4;
	}
}

float PhaseIndexProfile::getPercentAgeG2toM(int x, int d) {
	// calculate transition times from phase indices assuming expontential age distribution
	float G1 = G1_->calculateInterpolatedValue(x,d);
	float S  = S_->calculateInterpolatedValue(x,d);
	float G2 = G2_->calculateInterpolatedValue(x,d);
	return (logf(2)-logf(2-(G1+S+G2)))/logf(2);
}


/*!
 * Calculates DNA content of cell in row r assuming exponential age distribution
 */
float PhaseIndexProfile::getDna(float percentAge, int x, int d) {
	// calculate transition times from phase indices assuming exponential age distribution
	float G1 = G1_->calculateInterpolatedValue(x,d);
	float S  = S_->calculateInterpolatedValue(x,d);
	float transition_G1toS = (logf(2)-logf(2-G1))/logf(2);
	float transition_StoG2 = (logf(2)-logf(2-(G1+S)))/logf(2);

	// DNA content at G1 is 1
	if (percentAge<transition_G1toS) {
		return 1;
	}
	// DNA content at G2, M is 2
	else if (percentAge>=transition_StoG2) {
		return 2;
	}
	// DNA content at S is such that dna is 1 at beginning of S, 2 at end of S
	else {
		float percentS_complete = (percentAge-transition_G1toS)/(transition_StoG2-transition_G1toS);
		return 1+percentS_complete;
	}

}

/*!
 * Determine whether a cell in row r is in M-phase assuming exponential age distribution
 */
bool PhaseIndexProfile::isMphase(float percentAge, int x, int d) {
	float G1 = G1_->calculateInterpolatedValue(x,d);
	float S  = S_->calculateInterpolatedValue(x,d);
	float G2 = G2_->calculateInterpolatedValue(x,d);
	float transition_G2toM = (logf(2)-logf(2-(G1+S+G2)))/logf(2);
	if (percentAge>=transition_G2toM) {
		return true;
	} else {
		return false;
	}
}
