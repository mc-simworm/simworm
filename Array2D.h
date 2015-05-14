/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

/*!
 * \class Array2D
 * \brief Provides functionality to work with 2-dimensional arrays that store schedules.
 *
 * A common task in simworm is to store "schedules".  For example, you might need to store
 * the cell cycle profile schedule (i.e., the spatial cell cycle profile over time), or the
 * apoptosis schedule (the probability of apoptosis over time).  These schedules are stored
 * as a 1D array combined with an 2D array.
 *
 * For example, suppose the cell cycle profile is 5,4,3,3,3 at 10 cell divisions and 6,5,4,3,3
 * at 20 cell divisions.  This schedule is stored as:
 *
 * v1 = [10,20]
 * v2 = [5,4,3,3,3;
 *       6,5,4,3,3]
 *
 * This class provides functionality to access and mutate v1,v2 as well as calculate useful values
 * like the cell cycle profile at intermediate times through interpolation.
 *
 */
#ifndef ARRAY2D_H_
#define ARRAY2D_H_

#include <boost/multi_array.hpp>
#include "definitions.h"
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/lexical_cast.hpp>
#include <stdio.h>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include "util.h"

template <class T>
class Array2D {
public:
	typedef boost::multi_array<T, 2> array_type;
	typedef typename array_type::index array_index;
	typedef typename array_type::size_type array_size;

private:
	array_type* __restrict__ v2;
	boost::multi_array<int, 1>* __restrict__ v1;

public:
	/*!
	 * Constructor
	 */
	Array2D() {
		v2 = new array_type(boost::extents[0][0]);
		v1 = new boost::multi_array<int, 1>(boost::extents[0]);
	}

	/*
	 * Copy contructor
	 */
	Array2D(const Array2D& other) {
		// deep copy of boost matrix
		array_size dimx,dimy;
		other.getDimensions(dimx,dimy);
		v1 = new boost::multi_array<int,1>(boost::extents[dimy]);
		v2 = new array_type(boost::extents[dimx][dimy]);
		(*v1) = *(other.v1);
		(*v2) = *(other.v2);
	}

	/*!
	 * Destructor
	 */
	~Array2D() {
		delete v1; v1=NULL;
		delete v2; v2=NULL;
	}

	/*!
	 * Write matrix values to standard output
	 */
	void display() {
		array_index dimx,dimy;
		getDimensions(dimx,dimy);
		std::cout<<"Dimensions: dimx="<<dimx<<", dimy="<<dimy<<std::endl;
		for (array_index y=0;y<dimy;y++) {
			std::cout<<(*v1)[y]<<": ";
			for (array_index x=0;x<dimx;x++) {
				std::cout<<(*v2)[x][y]<<",";
			}
			std::cout<<std::endl;
		}
	}

	/*!
	 * Write matrix values to standard output in vector format (i.e., all elements in one row)
	 */
	void displayVector() {
		array_index dimx,dimy;
		getDimensions(dimx,dimy);
		for (array_index y=0;y<dimy;y++) {
			for (array_index x=0;x<dimx;x++) {
				std::cout<<(*v2)[x][y];
				if (y!=dimy-1 || x!=dimx-1) {
					std::cout<<",";
				}
			}
		}
	}

	/*!
	 * Calculate mean of a given matrix row y
	 */
	T mean(array_index y) {
		array_index dimx,dimy;
		getDimensions(dimx,dimy);
		T sum=0;
		for (array_index x=0;x<dimx;x++) {
			sum+=(*v2)[x][y];
		}
		return sum/(T)dimx;
	}

	/*!
	 * Calculate mean of a given matrix row y
	 */
	T mean_x(array_index x) {
		array_index dimx,dimy;
		getDimensions(dimx,dimy);
		T sum=0;
		for (array_index y=0;y<dimy;y++) {
			sum+=(*v2)[x][y];
		}
		return sum/(T)dimy;
	}

	/*!
	 * Access/mutate matrix values
	 */
	T& operator() (array_index x,array_index y) {
		#if defined INCLUDESANITYCHECKS
		// Sanity check
		array_index dimx,dimy;
		getDimensions(dimx,dimy);
		if (x<0 || x>=dimx || y<0 || y>=dimy) {
			printf("Coordinates (%d,%d) outside of image range (%d,%d)",x,y,dimx,dimy);
			exit(1);
		}
		#endif
		return (*v2)[x][y];
	}

	/*!
	 * Access/mutate matrix value
	 */
	int& operator() (array_index y) {
		#if defined INCLUDESANITYCHECKS
		// Sanity check
		array_index dimx,dimy;
		getDimensions(dimx,dimy);
		if (y<0 || y>=dimy) {
			printf("Coordinate (%d) outside of image range (%d)",y,dimy);
			exit(1);
		}
		#endif
		return (*v1)[y];
	}

	/*!
	 * Deep copy
	 */
	void operator= (const Array2D& other) {
		// deep copy of boost matrix
		array_index dimx,dimy;
		other.getDimensions(dimx,dimy);
		v2->resize(boost::extents[dimx][dimy]);
		(*v2) = *(other.v2);
		v1->resize(boost::extents[dimy]);
		(*v1) = *(other.v1);
	}

	/*!
	 * Parse formatted string into matrix
	 */
	void parseString(std::string str) {
		if (str.find('*')==std::string::npos) {
			std::vector<std::string> split2;
			boost::split(split2,str,boost::is_any_of("/"));
			array_index dimy = (array_index) split2.size();
			v1->resize(boost::extents[dimy]);
			for (array_index y=0;y<dimy;y++) {
				std::vector<std::string> split3;
				boost::split(split3,split2[(array_size)y],boost::is_any_of(":"));
				(*v1)[y]=atoi(split3[0].c_str());
				std::vector<std::string> split4;
				boost::split(split4,split3[1],boost::is_any_of(","));
				array_index dimx = (array_index) split4.size();
				if (y==0) {
					v2->resize(boost::extents[dimx][dimy]);
				}
				for (array_index x=0;x<dimx;x++) {
					(*v2)[x][y] = (T)atof(split4[(array_size) x].c_str());

				}
			}
		} else {
			std::vector<std::string> split2;
			boost::split(split2,str,boost::is_any_of("/"));
			array_index dimy = (array_index) split2.size();
			v1->resize(boost::extents[dimy]);
			// initialize v2
			int dimx = -1;
			for (array_index y=0;y<dimy;y++) {
				std::vector<std::string> split3;
				boost::split(split3,split2[(array_size)y],boost::is_any_of(":"));
				(*v1)[y]=atoi(split3[0].c_str());
				std::vector<std::string> split4;
				boost::split(split4,split3[1],boost::is_any_of("*"));
				dimx = std::max(dimx,atoi(split4[0].c_str()));
			}
			v2->resize(boost::extents[dimx][dimy]);

			// fill up matrix
			for (array_index y=0;y<dimy;y++) {
				std::vector<std::string> split3;
				boost::split(split3,split2[(array_size)y],boost::is_any_of(":"));
				(*v1)[y]=atoi(split3[0].c_str());
				std::vector<std::string> split4;
				boost::split(split4,split3[1],boost::is_any_of("*"));
				array_index dim_x = atoi(split4[0].c_str());
				for (array_index x=0;x<dim_x;x++) {
					(*v2)[x][y] = (T)atof(split4[1].c_str());
				}
			}
		}
	}

	/*!
	 * Return x,y,z dimensions of the boost matrix.
	 */
	void getDimensions(array_size &dimx, array_size &dimy) const {
		const boost::multi_array_types::size_type* dim = v2->shape();
		dimx = dim[0];
		dimy = dim[1];
		#if defined INCLUDESANITYCHECKS
		// Sanity check
		const boost::multi_array_types::size_type* dim_sanity = v1->shape();
		array_size dimy_sanity=dim_sanity[0];
		if (dimy!=dimy_sanity) {
			printf("v1,v2 dimensions inconsistent (%d, %d)\n",dimy,dimy_sanity);
			exit(1);
		}
		#endif
	}

	/*!
	 * Return x,y,z dimensions of the boost matrix.
	 */
	void getDimensions(array_index &dimx, array_index &dimy) const {
		array_size x,y;
		getDimensions(x,y);
		dimx=(array_index) x;
		dimy=(array_index) y;
	}

	/*!
	 * Return x,y,z dimensions of the boost matrix.
	 */
	void getDimensions(int &dimx, int &dimy) const {
		array_size x,y;
		getDimensions(x,y);
		dimx=(int) x;
		dimy=(int) y;
	}

	/*!
	 * Resize matrix
	 */
	void resize_unsigned(array_size dimx, array_size dimy) {
		v2->resize(boost::extents[(array_index) dimx][(array_index) dimy]);
		v1->resize(boost::extents[(array_index) dimy]);
	}

	/*!
	 * Resize matrix
	 */
	void resize(array_index dimx, array_index dimy){
		resize_unsigned((array_size) dimx, (array_size) dimy);
	}

	/*!
	 * Calculate intermediate matrix values through linear interpolation
	 */
	T calculateInterpolatedValue(array_index x, int t) {
		array_index dimx,dimy;
		getDimensions(dimx,dimy);
		#if defined INCLUDESANITYCHECKS
		if (x>=dimx || x<0) {
			printf("x exceeds bounds (x=%d, dimx=%d",x,dimx);
			exit(1);
		}
		#endif
		if (t<=(*v1)[0]) {
			return (*v2)[x][0];
		} else if (t>=(*v1)[dimy-1]) {
			return (*v2)[x][dimy-1];
		} else {
			float rn=0;
			for (array_index y=0;y<dimy-1;y++) {
				if (t>=(*v1)[y] && t<(*v1)[y+1]) {
					rn = interpolation_linear((float)t,(float)(*v1)[y],(float)(*v1)[y+1],(float)(*v2)[x][y],(float)(*v2)[x][y+1]);
				}
			}
			return (T)rn;
		}
	}

	/*!
	 * Calculate intermediate matrix values through floor
	 */
	T calculateFloorValue(array_index x, int t) {
		array_index dimx,dimy;
		getDimensions(dimx,dimy);
		#if defined INCLUDESANITYCHECKS
		if (x>=dimx || x<0) {
			printf("x exceeds bounds (x=%d, dimx=%d",x,dimx);
			exit(1);
		}
		#endif
		for (array_index y=0;y<dimy-1;y++) {
			if (t>=(*v1)[y] && t<(*v1)[y+1]) {
				return (*v2)[x][y];
			}
		}
		if (t<=(*v1)[0]) {
			return (*v2)[x][0];
		} else if (t>=(*v1)[dimy-1]) {
			return (*v2)[x][dimy-1];
		} else {
			printf("calculateFloorValue not working correctly\n");
			exit(1);
		}
	}

	/*!
	 * Fill matrix with specified values
	 */
	void fill(T val) {
		array_index dimx,dimy;
		getDimensions(dimx,dimy);
		for (array_index x=0;x<dimx;x++) {
			for (array_index y=0;y<dimy;y++) {
				(*v2)[x][y]=val;
			}
		}
	}
};

#endif
