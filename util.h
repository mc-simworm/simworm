/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/numeric/conversion/bounds.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/linear_congruential.hpp>
#include <sys/time.h>
#include <fstream>
#include <vector>
#include "definitions.h"
#include <errno.h>
typedef boost::minstd_rand base_generator_type;

float interpolation_linear(float x, float x0, float x1, float y0, float y1);
void initializeRandomNumberGenerator(int randomSeed, base_generator_type *generator);
float randUniformFloat (float lower, float upper, base_generator_type *generator);
int minCircularDistance(int y1, int y2, int size,base_generator_type *generator);
float percentAgeFromExponentialDistribution(base_generator_type *generator);
int randUniformInt(int lower, int upper,base_generator_type *generator);
float randNormal(float mean, float std,base_generator_type *generator);
int calculateSmallestIntegerGreaterThan(float f);
int calculateLargestIntegerLessEqualThan(float f);
base_generator_type& myGenerator();
std::vector<unsigned int> permutedVector(size_t sz,base_generator_type *generator);

