
/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#ifndef GEOMETRYPROFILE_H_
#define GEOMETRYPROFILE_H_
#include <vector>
#include <iostream>
#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <boost/algorithm/string.hpp>
#include "Array1D.h"
#include "Array2D.h"
#include "definitions.h"

/*!
 * \class GeometryProfile
 * \brief Provides functionality to work with mitotic region geometry profiles
 */
class GeometryProfile {
private:
	Array2D<int> *geometryProfile_;

public:
	GeometryProfile();
	virtual ~GeometryProfile();
	void parseFormattedParameters(std::string formattedParameters);
	void display();
	int maxSizeOfCellRows();
	int maxSizeOfMR();
	int numCellRows();
	int totalNumberOfCells(int cellDivision);
	int calculateGeometry(int r, int d);
	void getCurrentGeometry(int d, Array1D<int> *currentGeometry);
	bool geometryConstant(int d1, int d2);
	bool isMonotonicIncreasing();
};

#endif /* GEOMETRYPROFILE_H_ */
