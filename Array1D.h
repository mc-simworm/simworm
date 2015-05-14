/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/


/*!
 * \class Array1D
 * \brief Provides functionality to work with 1-dimensional arrays
 */
#ifndef ARRAY1D_H_
#define ARRAY1D_H_

#include <boost/multi_array.hpp>
#include "definitions.h"
#include <boost/algorithm/string.hpp>
#include "Array2D.h"
#include <boost/bind.hpp>

template <class T>
class Array1D {

public:
	typedef boost::multi_array<T, 1> array_type;
	typedef typename array_type::index array_index;
	typedef typename array_type::size_type array_size;

private:
	array_type* __restrict__ vec;

public:

	/*!
	 * Constructor
	 */
	Array1D() {
		vec = new boost::multi_array<T, 1>(boost::extents[0]);
	}

	/*
	 * Copy contructor
	 */
	Array1D(const Array1D& other) {
		// deep copy of boost matrix
		array_index sz = other.size();
		vec = new array_type(boost::extents[sz]);
		(*vec) = *(other.vec);
	}

	/*!
	 * Destructor
	 */
	~Array1D() {
		delete vec;
		vec=NULL;
	}

	/*!
	 * Deep copy
	 */
	void operator= (const Array1D& other) {
		// deep copy of boost matrix
		array_index sz = other.size();
		vec->resize(boost::extents[sz]);
		(*vec) = *(other.vec);
	}

	/*!
	 * Copy vector into Array1D
	 */
	void copy(std::vector<T> *other) {
		array_size sz = other->size();
		vec->resize(boost::extents[sz]);
		for (array_index i=0;i<sz;i++) {
			(*vec)[i] = (*other)[i];
		}
	}

	/*!
	 * Reshape 1D array into a 2D matrix
	 */
	void reshape(array_size dimx, array_size dimy, Array2D<T> *v2) {
		#if defined INCLUDESANITYCHECKS
		array_size sz = size();
		// Sanity check
		if (sz!=dimx*dimy) {
			printf("Reshaping with inconsistent sizes (%d vs %d,%d)",sz,dimx,dimy);
			exit(1);
		}
		#endif
		v2->resize_unsigned(dimx,dimy);
		int counter=0;
		for (array_index y=0;y< (array_index) dimy;y++) {
			for (array_index x=0; x< (array_index) dimx; x++) {
				(*v2)(x,y)=(*vec)[counter];
				counter++;
			}
		}
	}

	/*!
	 * Reshape 1D array into a 2D matrix
	 */
	void reshape(array_index dimx, array_index dimy, Array2D<T> *v2) {
		#if defined INCLUDESANITYCHECKS
		array_size sz = size();
		// Sanity check
		if (sz!=dimx*dimy) {
			printf("Reshaping with inconsistent sizes (%d vs %d,%d)",sz,dimx,dimy);
			exit(1);
		}
		#endif
		v2->resize(dimx,dimy);
		int counter=0;
		for (array_index y=0;y<dimy;y++) {
			for (array_index x=0; x<dimx; x++) {
				(*v2)(x,y)=(*vec)[counter];
				counter++;
			}
		}
	}


	/*!
	 * Access/mutate element idx
	 */
	T& operator() (array_index idx) {
		#if defined INCLUDESANITYCHECKS
		// Sanity check
		array_size sz = size();
		if (idx>=sz) {
			printf("Index (%d) outside of array range (%d)",idx,sz);
			exit(1);
		}
		#endif
		return (*vec)[idx];
	}

	/*!
	 * Print array to standard output
	 */
	void display() {
		array_index sz=size();
		for (array_index i=0;i<sz;i++) {
			std::cout<<(*vec)[i];
			if (i!=sz-1) {
				std::cout<<",";
			}
		}
		std::cout << std::endl;
	}

	/*!
	 * Resize array
	 */
	void resize_unsigned(array_size sz) {
		vec->resize(boost::extents[(array_index) sz]);
	}

	/*!
	 * Resize array
	 */
	void resize(array_index sz){
		resize_unsigned((array_size) sz);
	}

	void resize(unsigned long sz) {
		resize_unsigned((array_size) sz);
	}

	void resize(int sz) {
		resize_unsigned((array_size) sz);
	}

	/*!
	 * Sort array in ascending order
	 */
	void sort() {
		std::sort(vec->begin(),vec->end());
	}

	/*!
	 * Sort vector in ascending order and returns sorted indices
	 */
	void sort(Array1D<int> *idx) {
		// resize
		idx->resize_unsigned(size_unsigned());

		// make vector of paired values,indices
		std::vector<std::pair<T,int> > vectorPair;
		array_index sz=size();
		for (array_index i=0;i<sz;i++) {
			vectorPair.push_back(std::make_pair((*vec)[i],i));
		}
		std::sort(vectorPair.begin(), vectorPair.end(),
		          boost::bind(&std::pair<T, int>::first, _1) < boost::bind(&std::pair<T, int>::first, _2));

		typename std::vector<std::pair<T,int> >::const_iterator it;
		array_type vec_sorted(boost::extents[size()]);
		int counter=0;
		for( it = vectorPair.begin(); it != vectorPair.end(); it++ ) {
			std::pair<T,int> p = *it;
			(*idx)(counter) = p.second;
			vec_sorted[counter] = p.first;
			counter++;
		}
		(*vec) = vec_sorted;
	}

	/*!
	 * Find minimum of vector
	 */
	T min() {
		T minimumValue = INF;
		for (array_index i=0;i<size();i++) {
			minimumValue=std::min(minimumValue,(*vec)[i]);
		}
		return minimumValue;
	}

	/*!
	 * Find maximum of vector
	 */
	T max() {
		T maximumValue = -INF;
		for (array_index i=0;i<size();i++) {
			maximumValue=std::max(maximumValue,(*vec)[i]);
		}
		return maximumValue;
	}

	/*!
	 * Find pth percentile of vector
	 */
	T perctile(double percentile) {
		// get index corresponding to given percentile
		int idxPerc = (int)(percentile*(double)(vec->size()-1));

		// get copy of vector
		array_index sz = size();
		array_type temp(boost::extents[sz]);
		temp = *vec;

		// sort copy and return appropriate percentile
		std::sort(temp.begin(),temp.end());
		return temp[idxPerc];
	}

	/*!
	 * Return number of elements in boost matrix
	 */
	array_size size_unsigned() const {
		return vec->size();
	}

	/*!
	 * Return number of elements in boost matrix
	 */
	array_index size() const {
		return (array_index) vec->size();
	}

	/*!
	 * Fill vector up with specified value
	 */
	void fill(T val) {
		array_index sz = size();
		for (array_index i=0;i<sz;i++) {
			(*vec)[i]=val;
		}
	}

	/*!
	 * Parse comma delimited std::string into array
	 */
	void parseString(std::string str) {
		// split std::string
		std::vector<std::string> split;
		boost::split(split,str,boost::is_any_of(","));

		// resize matrix
		array_size sz = split.size();
		resize_unsigned(sz);

		// Read string contents and put in vector
		for (unsigned int i=0;i<sz;i++) {
			(*vec)[i]=(T)(atof(split[i].c_str()));
		}

		/*
		element *vec_it=vec->origin();
		for (auto split_it=split.begin(); split_it!=split.end(); ++split_it, ++vec_it) {
			*vec_it = (T)(atof(split_it->c_str()));
		}
		*/
	}

	/*!
	 * Calculate mean value of array
	 */
	T mean() {
		T sum=0;
		array_index sz=size();
		for (array_index i=0;i<sz;i++) {
			sum+=(*vec)[i];
		}
		return sum/(T)size();
	}

	/*!
	 * Append value to array
	 */
    void push_back(T val) {
		array_index sz = size();
		vec->resize(boost::extents[sz+1]);
		(*vec)[sz] = val;
    }

    /*!
     * Return the last element of the array
     */
    T back() {
		array_size sz = size();
		// Sanity check
		if (sz==0) {
			std::cerr<<"Empty vector"<<std::endl;
			exit(1);
		}
		return (*vec)[sz-1];
    }
};

#endif
