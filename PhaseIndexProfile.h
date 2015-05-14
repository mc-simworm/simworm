/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

/*!
 * \class PhaseIndexProfile
 * \brief Provides functionality to work with phase index profiles
 */
#ifndef PHASEINDEXPROFILE_H_
#define PHASEINDEXPROFILE_H_
#include <vector>
#include <iostream>
#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <boost/algorithm/string.hpp>
//#include "Array1D.h"
#include "Array2D.h"
#include "definitions.h"
#include "GeometryProfile.h"

class PhaseIndexProfile {
private:
	Array2D<float> *G1_;	// spatial profile of G1 indices
	Array2D<float> *S_;		// spatial profile of S indices
	Array2D<float> *G2_;	// spatial profile of G2 indices
	Array2D<float> *M_;		// spatial profile of M indices

public:
	PhaseIndexProfile();
	virtual ~PhaseIndexProfile();
	float getMphaseIndex(int x, int d);
	void parseFormattedParameters(std::string formattedParameters, std::string mitoticRegionGeometry);
	void display();
	int getPhase(float percentAge, int x, int d);
	float getDna(float percentAge, int x, int d);
	bool isMphase(float percentAge,int x, int d);
	float getPercentAgeG2toM(int x, int d);
};

#endif /* PHASEINDEXPROFILE_H_ */
