/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

/*!
 * Class: FitData
 * Provides cell cycle fit functionality
 */

#ifndef FITDATA_H_
#define FITDATA_H_

#include "Array1D.h"
#include "Array2D.h"
#include "util_germline.h"
#include <boost/multi_array.hpp>
typedef boost::multi_array<int, 4> array4D_int;
typedef boost::multi_array<float, 3> array3D_float;


class FitData {

private:

	Array1D<float> *chaseTimes_;
	Array1D<float> *dnaBin_;
	array4D_int *EMM_data_;
	array3D_float *FLM_data_;
	long numRows_;

public:

	FitData();
	virtual ~FitData();
	void initializeStorage(int numRows, Array1D<float> *dnaBin, Array1D<float> *chaseTimes);
	void parseData(int chaseTimeIdx, Array2D<int> *position, Array2D<float> *dna, Array2D<int> *edu, Array2D<int> *mPhase);
	void serializeOutput();
	void parseString(std::string str);

};

#endif /* FITDATA_H_ */
