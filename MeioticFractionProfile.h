/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/


/*!
 * \class MeioticFractionProfile
 * \brief Provides functionality to work with spatiotemporal meiotic fraction profiles
 */
#ifndef MEIOTICREGIONPROFILE_H_
#define MEIOTICREGIONPROFILE_H_
#include <vector>
#include <iostream>
#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <boost/algorithm/string.hpp>
#include "Array1D.h"
#include "Array2D.h"
#include <boost/math/special_functions/round.hpp>
#include "GeometryProfile.h"

class MeioticFractionProfile {
private:
	Array2D<float> *meioticFractionProfile_;

public:
	MeioticFractionProfile();
	virtual ~MeioticFractionProfile();
	void getCurrentProfile(int d, Array1D<float> *currentProfile);
	void display();
	float calculateMeioticFraction(int r, int d);
	void parseFormattedParameters(std::string formattedParameters, std::string mitoticRegionGeometry);
	float calculateMphaseThreshold(int r, int d);
};
#endif /* MEIOTICREGIONPROFILE_H_ */
