/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/


/*!
 * \class CellCycleProfile
 * \brief Provides functionality to work with spatiotemporal cell cycle profiles
 */
#ifndef CELLCYCLELENGTH_H_
#define CELLCYCLELENGTH_H_
#include <vector>
#include <iostream>
#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <boost/algorithm/string.hpp>
#include "Array1D.h"
#include "Array2D.h"
#include <boost/math/special_functions/round.hpp>

class CellCycleProfile {
private:
	Array2D<float> *cellCycleProfile_;

public:
	CellCycleProfile();
	virtual ~CellCycleProfile();
	void parseFormattedParameters(std::string formattedParameters, int numCellRows);
	void display();
	float calculateCellCycleLength(int r, int d);
	void getCurrentProfile(int d, Array1D<float> *currentProfile);
	float getCurrentProfileMinimum(int d);
	int getNumberControlPoints_t();
	bool profileConstant(int d1, int d2);
};

#endif /* CELLCYCLELENGTH_H_ */
