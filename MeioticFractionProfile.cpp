/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#include "MeioticFractionProfile.h"
using namespace boost;
using namespace std;

/*!
 * Constructor
 */
MeioticFractionProfile::MeioticFractionProfile() {
	meioticFractionProfile_ = new Array2D<float>();
}

/*!
 * Destructor
 */
MeioticFractionProfile::~MeioticFractionProfile() {
	delete meioticFractionProfile_; meioticFractionProfile_=NULL;
}

/*!
 * Interpolate to find the spatial cell cycle profile at cell at time d (cell divisions)
 */
void MeioticFractionProfile::getCurrentProfile(int d, Array1D<float> *currentProfile) {
	int dimx, dimy;
	meioticFractionProfile_->getDimensions(dimx,dimy);
	currentProfile->resize(dimx);
	for (int x=0;x<dimx;x++) {
		(*currentProfile)(x) = calculateMeioticFraction(x,d);
	}
}

/*!
 * Write cell cycle profile to standard output
 */
void MeioticFractionProfile::display() {
	cout << "Spatiotemporal meiotic fraction profile:" << endl;
	meioticFractionProfile_->display();
}

/*!
 * Calculate meiotic fraction at row r at cell division d
 */
float MeioticFractionProfile::calculateMeioticFraction(int r, int d) {
	return meioticFractionProfile_->calculateInterpolatedValue(r,d);
}

/*!
 * Calculate the M-phase threshold at row r at cell division d
 * If the simulated M-phase index is above the M-phase threshold, then G2-M transitions are pre-meiotic
 * If the simulated M-phase index is belwo the M-phase threshold, then G2-M transitions are mitotic
 */
float MeioticFractionProfile::calculateMphaseThreshold(int r, int d) {
	return meioticFractionProfile_->calculateInterpolatedValue(r,d);
}

/*
 * Parse formatted parameters.
 */
void MeioticFractionProfile::parseFormattedParameters(string serializedParameters, string mitoticRegionGeometry) {
	if (serializedParameters.find(":")==std::string::npos) {
		GeometryProfile gg;
		gg.parseFormattedParameters(mitoticRegionGeometry);
		int dimx = gg.numCellRows();
		meioticFractionProfile_->resize(dimx,1);
		meioticFractionProfile_->fill(boost::lexical_cast<float>(serializedParameters.c_str()));
	} else {
		meioticFractionProfile_->parseString(serializedParameters);
	}
}

/*
 * Code below is old way of formatting meioticFractionProfile parameters (ex. 0:1,11,19)
 *
void MeioticFractionProfile::parseFormattedParameters(string serializedParameters, int numCellRows) {
	// get spatiotemporal control points
	Array2D<float> ctrlPts;
	ctrlPts.parseString(serializedParameters);
	Array2D<float>::array_index dimx, dimy;
	ctrlPts.getDimensions(dimx,dimy);
	// shift from indexing at 1 to indexing at 0
	for (int x=0;x<dimx;x++) {
		for (int y=0;y<dimy;y++) {
			ctrlPts(x,y)-=1;
		}
	}

	// Sanity check
	if (dimx!=3) {
		cerr<<"Should have 3 spatial control points for meiotic fraction profile"<<endl;
		exit(1);
	}

	// initialize meioticFractionProfile_
	meioticFractionProfile_->resize(numCellRows,dimy);

	// Get explicit meiotic fraction for each cell row
	for (Array2D<float>::array_index y=0;y<dimy;y++) {
		(*meioticFractionProfile_)(y)=ctrlPts(y);
		for (int r=0;r<numCellRows;r++) {
			if (r<=ctrlPts(1,y)) {
				(*meioticFractionProfile_)(r,y)=0;
			} else if (r>=ctrlPts(2,y)) {
				(*meioticFractionProfile_)(r,y)=1;
			} else {
				float m = 1/(ctrlPts(2,y)-ctrlPts(1,y));
				float dr = (float)r - ctrlPts(1,y);
				float dy = m*dr;
				(*meioticFractionProfile_)(r,y)=dy;
			}
		}
	}
}
*/

