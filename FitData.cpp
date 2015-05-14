/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#include "FitData.h"
using namespace std;

FitData::FitData() {
	EMM_data_ = new array4D_int(boost::extents[0][0][0][0]);
	FLM_data_ = new array3D_float(boost::extents[0][0][0]);
	chaseTimes_ = new Array1D<float>;
	dnaBin_ = new Array1D<float>;
}

FitData::~FitData() {
	delete EMM_data_; EMM_data_=NULL;
	delete FLM_data_; FLM_data_=NULL;
	delete chaseTimes_; chaseTimes_=NULL;
	delete dnaBin_; dnaBin_=NULL;
}

void FitData::parseString(std::string str) {
	vector<string> split2;
	boost::split(split2,str,boost::is_any_of("/"));
	string chaseTimes = split2[0];
	string dnaBin = split2[1];
	string emm = split2[2];
	string flm = split2[3];

	chaseTimes_->parseString(chaseTimes);
	dnaBin_->parseString(dnaBin);

	Array1D<int> emm_vec;
	emm_vec.parseString(emm);
	numRows_ = emm_vec.size()/(dnaBin_->size()*2*chaseTimes_->size());

	EMM_data_->resize(boost::extents[numRows_][dnaBin_->size()][2][chaseTimes_->size()]);
	int counter=0;
	for (int t=0;t<chaseTimes_->size();t++) {
		for (int e=0;e<2;e++) {
			for (int d=0;d<dnaBin_->size();d++) {
				for (int r=0;r<numRows_;r++) {
					(*EMM_data_)[r][d][e][t] = emm_vec(counter);
					counter++;
				}
			}
		}
	}

	Array1D<float> flm_vec;
	flm_vec.parseString(flm);
	FLM_data_->resize(boost::extents[numRows_][2][chaseTimes_->size()]);
	counter=0;
	for (int t=0;t<chaseTimes_->size();t++) {
		for (int e=0;e<2;e++) {
			for (int r=0;r<numRows_;r++) {
				(*FLM_data_)[r][e][t]=flm_vec(counter);
				counter++;
			}
		}
	}
}

void FitData::serializeOutput() {
	// print chaseTimes
	for (int i=0;i<chaseTimes_->size();i++) {
		cout << (*chaseTimes_)(i);
		if (i!=chaseTimes_->size()-1) {
			cout << ",";
		}
	}
	cout << "/";

	// print dnaBin
	for (int i=0;i<dnaBin_->size();i++) {
		cout << (*dnaBin_)(i);
		if (i!=dnaBin_->size()-1) {
			cout << ",";
		}
	}
	cout << "/";

	// print EMM
	for (int t=0;t<chaseTimes_->size();t++) {
		for (int e=0;e<2;e++) {
			for (int d=0;d<dnaBin_->size();d++) {
				for (int r=0;r<numRows_;r++) {
					cout << (*EMM_data_)[r][d][e][t];
					if (t!=chaseTimes_->size()-1 || e!=1 || d!=dnaBin_->size() || r!=numRows_-1) {
						cout << ",";
					}
				}
			}
		}
	}
	cout << "/";

	// print FLM
	for (int t=0;t<chaseTimes_->size();t++) {
		for (int e=0;e<2;e++) {
			for (int r=0;r<numRows_;r++) {
				cout << (*FLM_data_)[r][e][t];
				if (t!=chaseTimes_->size()-1 || e!=1 || r!=numRows_-1) {
					cout << ",";
				}
			}
		}
	}
}

void FitData::parseData(int chaseTimeIdx, Array2D<int> *position, Array2D<float> *dna, Array2D<int> *edu, Array2D<int> *mPhase) {
	Array2D<int>::array_index dimx,dimy;
	position->getDimensions(dimx,dimy);
	for (int x=0;x<dimx;x++) {
		for (int y=0;y<dimy;y++) {
			int r = (*position)(x,y);
			float d = (*dna)(x,y);
			int e = (*edu)(x,y);
			int m = (*mPhase)(x,y);

			if (r<0) {r=0;}
			if (r>(int)numRows_-1) {r=(int)numRows_-1;}

			// Increment EMM histogram
			int dnaBinIdx=getClosestIdx(dnaBin_,d);
			(*EMM_data_)[r][dnaBinIdx][e][chaseTimeIdx]++;

			// Increment FLM histogram
			if (m==1) {
				(*FLM_data_)[r][e][chaseTimeIdx]++;
			}
		}
	}
}

void FitData::initializeStorage(int numRows, Array1D<float> *dnaBin, Array1D<float> *chaseTimes) {
	numRows_=numRows;
	(*chaseTimes_)=(*chaseTimes);
	(*dnaBin_)=(*dnaBin);
	EMM_data_->resize(boost::extents[numRows_][dnaBin_->size()][2][chaseTimes_->size()]);
	FLM_data_->resize(boost::extents[numRows_][2][chaseTimes_->size()]);
}

