/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#include "GeometryProfile.h"
using namespace boost;
using namespace std;

/*!
 * Constructor
 */
GeometryProfile::GeometryProfile() {
	geometryProfile_ = new Array2D<int>();
}

/*!
 * Destructor
 */
GeometryProfile::~GeometryProfile() {
	delete geometryProfile_; geometryProfile_=NULL;
}

/*!
 *
 */
bool GeometryProfile::isMonotonicIncreasing() {
	int dimx,dimy;
	geometryProfile_->getDimensions(dimx,dimy);
	bool rn=true;
	for (int y=0;y<dimy;y++) {
		for (int x=0;x<dimx-1;x++) {
			if ((*geometryProfile_)(x,y)>(*geometryProfile_)(x+1,y)) {
				rn = false;
			}
		}
	}
	return rn;
}

/*!
 * Write profile to standard output
 */
void GeometryProfile::display() {
	cout << "Temporal geometry profile: " << endl;
	geometryProfile_->display();
}

/*!
 * Parse formatted string
 */
void GeometryProfile::parseFormattedParameters(std::string formattedParameters) {
	geometryProfile_->parseString(formattedParameters);
}

/*!
 * Return the maximum cell row size (over all times)
 */
int GeometryProfile::maxSizeOfCellRows() {
	Array2D<int>::array_index dimx,dimy;
	geometryProfile_->getDimensions(dimx,dimy);
	int maxVal=-INF;
	for (Array2D<int>::array_index x=0;x<dimx;x++) {
		for (Array2D<int>::array_index y=0;y<dimy;y++) {
			if ((*geometryProfile_)(x,y)>maxVal) {
				maxVal=(*geometryProfile_)(x,y);
			}
		}
	}
	return maxVal;
}

/*!
 * Return the peak number of cells (over all times
 */
int GeometryProfile::maxSizeOfMR() {
	Array2D<int>::array_index dimx,dimy;
	geometryProfile_->getDimensions(dimx,dimy);
	int maxVal=-INF;
	for (Array2D<int>::array_index y=0;y<dimy;y++) {
		int sum = 0;
		for (Array2D<int>::array_index x=0;x<dimx;x++) {
			sum+=(*geometryProfile_)(x,y);
		}
		if (sum>maxVal) {maxVal = sum;}
	}
	return maxVal;
}

/*!
 * Get the number of cell rows, including cell rows with zero cells
 */
int GeometryProfile::numCellRows() {
	Array2D<int>::array_size dimx,dimy;
	geometryProfile_->getDimensions(dimx,dimy);
	return (int)dimx;
}

/*!
 * Returns the total number of cells in the mitotic region at time "cellDivision"
 */
int GeometryProfile::totalNumberOfCells(int cellDivision) {
	Array1D<int> geometry;
	getCurrentGeometry(cellDivision,&geometry);
	int sum=0;
	for (int i=0;i<geometry.size();i++) {
		sum+=geometry(i);
	}
	return sum;
}

/*!
 * Calculate number of cells in row r at time d (cell divisions)
 */
int GeometryProfile::calculateGeometry(int r, int d) {
	return geometryProfile_->calculateInterpolatedValue(r,d);
}

/*!
 * Calculate geometry profile at time d (cell divisions)
 */
void GeometryProfile::getCurrentGeometry(int d, Array1D<int> *currentGeometry) {
	Array1D<int>::array_index dimx,dimy;
	geometryProfile_->getDimensions(dimx,dimy);
	currentGeometry->resize(dimx);
	for (Array1D<int>::array_index x=0;x<dimx;x++) {
		(*currentGeometry)(x) = calculateGeometry((int)x,d);
	}
}

/*!
 * Return true if geometry is the same from cell division d1 to d2.
 * Return false if geometry changes between d1, d2
 */
bool GeometryProfile::geometryConstant(int d1, int d2) {
	Array1D<int> geometry1,geometry2;
	getCurrentGeometry(d1,&geometry1);
	getCurrentGeometry(d2,&geometry2);
	for (int i=0;i<geometry1.size();i++) {
		if (geometry1(i)!=geometry2(i)) {
			return false;
		}
	}
	return true;
}
