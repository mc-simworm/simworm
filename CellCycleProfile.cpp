/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#include "CellCycleProfile.h"
using namespace boost;
using namespace std;

/*!
 * Constructor
 */
CellCycleProfile::CellCycleProfile() {
	cellCycleProfile_ = new Array2D<float>();
}

/*!
 * Destructor
 */
CellCycleProfile::~CellCycleProfile() {
	delete cellCycleProfile_; cellCycleProfile_=NULL;
}

/*!
 * Interpolate to find the spatial cell cycle profile at cell at time d (cell divisions)
 */
void CellCycleProfile::getCurrentProfile(int d, Array1D<float> *currentProfile) {
	int dimx, dimy;
	cellCycleProfile_->getDimensions(dimx,dimy);
	currentProfile->resize(dimx);
	for (int x=0;x<dimx;x++) {
		(*currentProfile)(x) = calculateCellCycleLength(x,d);
	}
}

/*!
 * Get minimum cell cycle length across all cell rows at time d (cell divisions)
 */
float CellCycleProfile::getCurrentProfileMinimum(int d) {
	int dimx,dimy;
	cellCycleProfile_->getDimensions(dimx,dimy);
	Array1D<float> currentProfile;
	currentProfile.resize(dimx);
	for (int x=0;x<dimx;x++) {
		currentProfile(x) = calculateCellCycleLength(x,d);
	}
	return currentProfile.min();
}

/*!
 * Return number of temporal control points
 */
int CellCycleProfile::getNumberControlPoints_t() {
	Array2D<float>::array_size dimx, dimy;
	cellCycleProfile_->getDimensions(dimx,dimy);
	return (int)dimy;
}

/*!
 * Return true if geometry is the same from cell division d1 to d2.
 * Return false if geometry changes between d1, d2
 */
bool CellCycleProfile::profileConstant(int d1, int d2) {
	Array1D<float> profile1,profile2;
	getCurrentProfile(d1,&profile1);
	getCurrentProfile(d2,&profile2);
	for (int i=0;i<profile1.size();i++) {
		if (abs(profile1(i)-profile2(i))>EPS) {
			return false;
		}
	}
	return true;
}

/*!
 * Write cell cycle profile to standard output
 */
void CellCycleProfile::display() {
	cout << "Spatiotemporal cell cycle profile:" << endl;
	cellCycleProfile_->display();
}

/*!
 * Calculate cell cycle length at row r at cell division d
 */
float CellCycleProfile::calculateCellCycleLength(int r, int d) {
	return cellCycleProfile_->calculateInterpolatedValue(r,d);
}

/*
 * Parse formatted cell cycle parameters.
 */
void CellCycleProfile::parseFormattedParameters(string serializedParameters, int numCellRows) {
	// split input string
	vector<string> split1;
	split(split1,serializedParameters,is_any_of("*"));

	// get spatiotemporal control points
	Array2D<float> ccCtrlPts;
	boost::erase_all(split1[1], "{");
	boost::erase_all(split1[1], "}");
	ccCtrlPts.parseString(split1[1]);
	Array2D<float>::array_index dimx, dimy;
	ccCtrlPts.getDimensions(dimx,dimy);
	// shift from indexing at 1 to indexing at 0
	for (int x=0;x<dimx;x++) {
		for (int y=0;y<dimy;y++) {
			ccCtrlPts(x,y)-=1;
		}
	}

	// get values of cell cycle profile at control points
	Array1D<float> ccVal_1D;
	boost::erase_all(split1[0], "{");
	boost::erase_all(split1[0], "}");
	ccVal_1D.parseString(split1[0]);
	// Sanity check
	for (int i=0;i<ccVal_1D.size();i++) {
		if (abs(ccVal_1D(i))<EPS) {
			cerr<<"Cannot have 0 as cell cycle profile value"<<endl;
			exit(1);
		}
	}

	// reshape
	Array2D<float> ccVal;
	ccVal_1D.reshape(dimx,dimy,&ccVal);

	// get interpolation method
	string interpolationMethod=split1[2];

	// initialize cellCycleProfile_
	cellCycleProfile_->resize(numCellRows,dimy);

	// Get cell cycle length for each cell row
	for (Array2D<float>::array_index y=0;y<dimy;y++) {
		(*cellCycleProfile_)(y)=ccCtrlPts(y);
		for (int r=0;r<numCellRows;r++) {
			if (r<=ccCtrlPts(0,y)) {
				(*cellCycleProfile_)(r,y)=ccVal(0,y);
			} else if (r>=ccCtrlPts(dimx-1,y)) {
				(*cellCycleProfile_)(r,y)=ccVal(dimx-1,y);
			} else {
				for (Array2D<float>::array_index x=0;x<dimx-1;x++) {
					if (r>=ccCtrlPts(x,y) && r<ccCtrlPts(x+1,y)) {
						if (interpolationMethod.compare("step")==0) {
							(*cellCycleProfile_)(r,y)=ccVal(x,y);
						// For linear gradient, interpolate to get full matrix of cell cycle lengths
						} else if (interpolationMethod.compare("linear")==0) {
							float ccVal_interpolated = interpolation_linear((float)r,ccCtrlPts(x,y),ccCtrlPts(x+1,y),ccVal(x,y),ccVal(x+1,y));
							(*cellCycleProfile_)(r,y)=ccVal_interpolated;
						} else {
							printf("Invalid choice of interpolation %s\n",interpolationMethod.c_str());
							exit(1);
						}
					}
				}
			}
		}
	}
}


