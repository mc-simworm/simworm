/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#include "util.h"
#include <fcntl.h>
using namespace std;

//static base_generator_type *generator=NULL;

/*!
 * Return a permuted vector
 */
vector<unsigned int> permutedVector(size_t sz,base_generator_type *generator) {
	vector<unsigned int> v;
	for (unsigned int i=0; i<sz; i++) {
		v.push_back(i); // 0 1 2 3 4 5 6 7 8 9
	}

	vector<unsigned int> rn;
	for (unsigned int i=0;i<sz;i++) {
		int idx = randUniformInt(0,(int)v.size()-1,generator);
		rn.push_back(v[(unsigned int)idx]);
		v.erase(v.begin()+idx);
	}
	return rn;
}


/*!
 * Return a random number drawn from gaussian distribution with specified and mean, standard deviation
 */
float randNormal(float mean, float std,base_generator_type *generator) {
	if (abs(std)<EPS) {
		return mean;
	}
	boost::normal_distribution<> dist(mean, std);
	boost::variate_generator<base_generator_type&, boost::normal_distribution<> > gen(*generator,dist);
	return (float)gen();
}

/*!
 * Return random int drawn from uniform distribution with lower, upper bound
 */
int randUniformInt(int lower, int upper,base_generator_type *generator) {
	boost::uniform_int<> dist(lower,upper);
	boost::variate_generator<base_generator_type&, boost::uniform_int<> > gen(*generator,dist);
	// Sanity check
	if ((int)gen()<lower) {
		cerr<<"Overflow in random number generator"<<endl;
		exit(1);
	}
	return (int)gen();
}

/*!
 * Return a random float drawn from uniform distribution
 */
float randUniformFloat (float lower, float upper, base_generator_type *generator) {
	if (abs(lower-upper)<EPS) {
		return lower;
	}
	boost::uniform_real<> dist(lower,upper);
	boost::variate_generator<base_generator_type&, boost::uniform_real<> > gen(*generator,dist);
	return (float)gen();
}

/*!
 * Linear interpolation
 */
float interpolation_linear(float x, float x0, float x1, float y0, float y1) {
	return y0 + (x-x0)*(y1-y0)/(x1-x0);
}

/*!
 * Initialize random number generator
 */
void initializeRandomNumberGenerator(int randomSeed, base_generator_type *generator) {
	unsigned int seed = (unsigned int) randomSeed;
	if (randomSeed<0) {
		int randomDevice=open("/dev/urandom", O_RDONLY);
		if (randomDevice == -1){
			std::cerr<<"Could not open random device";
			exit(1);
		}
		ssize_t result=0;
		do {
			result=read(randomDevice, &seed, sizeof (unsigned int));
			if ((unsigned long) result < sizeof(unsigned int) && errno!=EINTR && errno!= EAGAIN){
				std::cerr<<"Problem reading from random device; error "<<errno<<std::endl;
				exit(1);
			}
		} while (result!=sizeof (unsigned int));
		close(randomDevice);
	}
	generator->seed(seed);
}


/*!
 * minimum circular distance from y1 to y2 assuming ring of size sz
 */
int minCircularDistance(int y1, int y2, int size,base_generator_type *generator){
	// get plus, minus circular distance from y1 to y2
	int shortest_distance, distance_plusDirection,distance_minusDirection;
	if (y2>y1) {
		distance_plusDirection = y2-y1;
		distance_minusDirection = size-y2+y1;
	} else if (y1>y2) {
		distance_plusDirection = size-y1+y2;
		distance_minusDirection = y1-y2;
	} else {
		distance_plusDirection = 0;
		distance_minusDirection= 0;
	}

	// get shortest circular distance
	if (distance_plusDirection<distance_minusDirection) {
		shortest_distance = distance_plusDirection;
	} else if (distance_minusDirection<distance_plusDirection) {
		shortest_distance = -distance_minusDirection;
	} else {
		if (randUniformFloat(0,1,generator)>0.5) {
			shortest_distance = distance_plusDirection;
		} else {
			shortest_distance = -distance_minusDirection;
		}
	}
	return shortest_distance;
}

/*!
 * Get percent age (between [0,1]) sampled from exponential distribution
 */
float percentAgeFromExponentialDistribution(base_generator_type *generator) {
	float percentAge=INF;
	while(percentAge>1) {
		float u = randUniformFloat(0,1,generator);
		float L = logf(2);
		percentAge = logf(1-u)/(-L);
	}
	return percentAge;
}


/*!
 * Calculate the smallest integer > f given f
 */
int calculateSmallestIntegerGreaterThan(float f) {
	float f2 = ceilf(f);
	if (f2>f) {
		return (int)f2;
	} else {
		return (int)f2+1;
	}
}

/*!
 * Calculate the smallest integer >= f given f
 */
int calculateLargestIntegerLessEqualThan(float f) {
	float f2 = floorf(f);
	return (int)f2;
}
