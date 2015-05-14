/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#include "GermLine.h"
using namespace std;

/*!
 * Constructor
 */
GermLine::GermLine(ProgramOptionParser *opt, int seed) {
	// get schedules
	geometry_ = new GeometryProfile;
	opt->getGeometryProfile(geometry_);
	cellCycle_ = new CellCycleProfile;
	opt->getCellCycleProfile(cellCycle_);
	phaseIndexProfile_ = new PhaseIndexProfile;
	opt->getPhaseIndicesProfile(phaseIndexProfile_);
	meioticFractionProfile_ = new MeioticFractionProfile;
	opt->getMeioticFractionProfile(meioticFractionProfile_);
	meioticRegionSize_ = new Array2D<int>;
	opt->getMeioticRegionSize(meioticRegionSize_);
	oocyteProbability_ = new Array2D<float>;
	opt->getOocyteProbability(oocyteProbability_);
	oocyteProbabilityMethod_ = opt->oocyteProbabilityInterpolationMethod();

	// initialize random number generator
	generator = new base_generator_type;
	initializeRandomNumberGenerator(seed, generator);

	// set simulation parameters
	cellCycleNoise_ = opt->noise;
	dnaNoise_ = opt->dnaNoise;
	movementProbabilityForward_ = opt->forwardMovementProbability;
	verbose_=opt->debug;
	pedigreeDepthIncrementMethod_=opt->pedigreeDepthIncrementMethod;

	// initialize timers
	countSperm_=0;
	countOocyteLaid_=0;
	tElapse_=0;
	cellDivisionsElapse_=0;
	cellsExitMeioticRegion_=0;
	countCellsEnterMeioticRegion_=0;

	// initialize storage vectors
	germCellVector_ = new vector<GermCell *>;
	meioticGermCellVector_ = new vector<int>;
	exitMeioticRegionTime_ = new vector<int>;
	cellProductionTimes_ = new vector<float>;
	pdGermlineAtEggLaying = new vector<float>;

	// set stopping conditions
	numSpermToProduce_=opt->numSperm;
	minCellDivisions_=opt->numDivisions;
	minNumOocytesProduced_=opt->numOocyte;
	minTimeElapsed_=opt->minTime;
	minCellsEnterMeioticRegion_ = 0;
	if (opt->fitnessMetric.find("meanPedigreeDepthMeioticCells:")!=string::npos || opt->fitnessMetric.find("tedNIPS:")!=string::npos) {
		vector<string> splitStr;
		boost::split(splitStr,opt->fitnessMetric,boost::is_any_of(":"));
		minCellsEnterMeioticRegion_ = atoi(splitStr[1].c_str());
	}

	// initialize mapPosition_ as matrix of -1 (no cells)
	mapPosition_ = new vector<vector<int>*>;
	Array1D<int> initialGeometry;
	geometry_->getCurrentGeometry(0,&initialGeometry);
	for (int i=0;i<initialGeometry.size();i++) {
		vector<int> *column = new vector<int>((unsigned int)initialGeometry(i));
		fill(column->begin(),column->end(),-1);
		mapPosition_->push_back(column);
	}
}

/*!
 * Destructor
 */
GermLine::~GermLine() {
	delete geometry_; geometry_=NULL;
	delete cellCycle_; cellCycle_=NULL;
	for (unsigned int i=0;i<mapPosition_->size();i++) {
		delete (*mapPosition_)[i]; (*mapPosition_)[i]=NULL;
	}
	delete mapPosition_; mapPosition_=NULL;
	for (unsigned int i=0;i<germCellVector_->size();i++) {
		delete (*germCellVector_)[i]; (*germCellVector_)[i]=NULL;
	}
	delete meioticGermCellVector_; meioticGermCellVector_=NULL;
	delete germCellVector_; germCellVector_=NULL;
	delete oocyteProbability_; oocyteProbability_=NULL;
	delete phaseIndexProfile_; phaseIndexProfile_=NULL;
	delete meioticFractionProfile_; meioticFractionProfile_=NULL;
	delete meioticRegionSize_; meioticRegionSize_=NULL;
	delete cellProductionTimes_; cellProductionTimes_=NULL;
	delete exitMeioticRegionTime_; exitMeioticRegionTime_=NULL;
	delete pdGermlineAtEggLaying; pdGermlineAtEggLaying=NULL;
	delete generator; generator=NULL;
}


/*!
 * Calculate cell production schedule, fecundity schedule, etc.
 */
void GermLine::calculateSchedules(int i, PedigreeDepthStorage *pdStorage) {
	Array1D<float> allCellsEnterMeioticRegionPedigreeDepth, spermProductionTimes, oocyteProductionTimes, cellProductionTimes, oocytePedigreeDepths, spermPedigreeDepths, spermProductionTimes_cellDiv, oocyteProductionTimes_cellDiv,pd_mitotic_region, pd_meiotic_region, pd_whole_germline,pd_atEggLaying;
	Array1D<int> cellExitTimes;
	calculateRawOutput(&allCellsEnterMeioticRegionPedigreeDepth,&oocyteProductionTimes,&spermProductionTimes,&oocytePedigreeDepths,&spermPedigreeDepths,&oocyteProductionTimes_cellDiv,&spermProductionTimes_cellDiv,&pd_mitotic_region,&pd_meiotic_region,&pd_whole_germline,&cellProductionTimes,&cellExitTimes,&pd_atEggLaying);

	// get cell production schedule
	Array2D<float>::array_index dummy,dimy;
	pdStorage->simulatedCellProductionSchedule->getDimensions(dummy,dimy);
	for (int j=0;j<dimy;j++) {
		float actualTime;
		int cellsProduced = (*pdStorage->cellProductionConstraints)(j);
		if (cellsProduced>cellProductionTimes.size()) {
			actualTime=INF;
		} else {
			actualTime = cellProductionTimes(cellsProduced-1);
		}
		(*pdStorage->simulatedCellProductionSchedule)(i,j)=actualTime;
	}

	// get fecundity schedule
	int spermProducedByConstraint = 0;
	for (int j=0;j<spermProductionTimes.size();j++) {
		float spermTimeConstraint = (*pdStorage->fecundityConstraints)(0,0);
		float t=spermProductionTimes(j);
		if (t>spermTimeConstraint) {
			break;
		}
		spermProducedByConstraint=j+1;
	}
	(*pdStorage->simulatedFecunditySchedule)(i,0)=(float)spermProducedByConstraint;
	Array2D<float>::array_index dimy_fecundity;
	pdStorage->fecundityConstraints->getDimensions(dummy,dimy_fecundity);
	for (int j=1; j<dimy_fecundity; j++) {
		float oocyteTimeConstraint = (*pdStorage->fecundityConstraints)(0,j);
		int numOocytesProduced=0;
		float lowerTimeConstraint = 0;
		if (j!=1) {
			lowerTimeConstraint = (*pdStorage->fecundityConstraints)(0,j-1);
		}
		for (int k=0;k<oocyteProductionTimes.size();k++) {
			float t = oocyteProductionTimes(k);
			if (t>lowerTimeConstraint && t<=oocyteTimeConstraint) {
				numOocytesProduced++;
			}
		}
		(*pdStorage->simulatedFecunditySchedule)(i,j)=(float)numOocytesProduced;
	}

	// get number of cell divisions and number of cells that hit end of meiotic region during time contraints
	for (int j=0; j<dimy_fecundity; j++) {
		float t2 = (*pdStorage->fecundityConstraints)(0,j);
		float t1;
		if (j==0) {
			t1 = 0;
		} else {
			t1 = (*pdStorage->fecundityConstraints)(0,j-1);
		}
		// get number of cell divisions
		int numDivisions = 0;
		for (int k=0; k<cellProductionTimes.size(); k++) {
			float divisionTime = cellProductionTimes(k);
			if (divisionTime>=t1 && divisionTime<t2) {
				numDivisions++;
			}
		}
		// get number of cells that exit the meiotic region (i.e., undergo apoptosis/oogenesis)
		int numHitEndMeioticRegion = 0;
		for (int k=0; k<cellExitTimes.size(); k++) {
			int numCellDivisions = cellExitTimes(k);
			float divisionTime = cellProductionTimes(numCellDivisions-1);
			if (divisionTime>=t1 && divisionTime<t2) {
				numHitEndMeioticRegion++;
			}
		}
		(*pdStorage->simulatedCellProductionSchedule_atFecundityBins)(i,j) = (float)numDivisions;
		(*pdStorage->simulatedMeioticExitSchedule)(i,j) = (float)numHitEndMeioticRegion;
	}
}


/*!
 * Calculates useful scalar quantities
 */
void GermLine::calculateOutput(float &pdNCellsEnterMeioticRegion, float &pdSperm, float &pdOocyte, float &pdSpermAndOocyte, float &pdSpermAndOocyteSum, float &pdMitoticRegion, float &pdMeioticRegion, float &pdWholeGermline, float &pdGermlineAtOogenesis, int &numCellsEnterMeioticRegion, int &numCellDivisions, float &timeFirstCellExitMeioticRegion, float &timeNthCellExitsMitoticRegion) {
	// get raw data from germline simulation
	Array1D<float> NCellsEnterMeioticRegionPedigreeDepth, spermProductionTimes, oocyteProductionTimes, cellProductionTimes, oocytePedigreeDepths, spermPedigreeDepths, spermProductionTimes_cellDiv, oocyteProductionTimes_cellDiv,pd_mitotic_region, pd_meiotic_region, pd_whole_germline,pd_atEggLaying;
	Array1D<int> cellExitTimes;
	calculateRawOutput(&NCellsEnterMeioticRegionPedigreeDepth,&oocyteProductionTimes,&spermProductionTimes,&oocytePedigreeDepths,&spermPedigreeDepths,&oocyteProductionTimes_cellDiv,&spermProductionTimes_cellDiv,&pd_mitotic_region,&pd_meiotic_region,&pd_whole_germline,&cellProductionTimes,&cellExitTimes,&pd_atEggLaying);

	// calculate average pedigree depth over first N cells that ever entered the meiotic region
	pdNCellsEnterMeioticRegion = NCellsEnterMeioticRegionPedigreeDepth.mean();

	// calculate average pedigree depth over all sperm produced
	pdSperm = spermPedigreeDepths.mean();

	// calculate average pedigree depth over all oocytes produced
	pdOocyte = oocytePedigreeDepths.mean();

	// calculate average pedigree depth over all embryos produced
	Array1D<float> spermOocytePedigreeDepths;
	spermOocytePedigreeDepths.resize(oocytePedigreeDepths.size());
	for (int i=0;i<oocytePedigreeDepths.size();i++) {
		Array1D<float>::array_index idxSperm=i/4; // divide by four to take into account each germ cell makes 4 sperm
		if (idxSperm>=spermPedigreeDepths.size()) {
			idxSperm = spermPedigreeDepths.size() - 1;
		}
		spermOocytePedigreeDepths(i) = (oocytePedigreeDepths(i) + spermPedigreeDepths(idxSperm))/2;
	}
	pdSpermAndOocyte = spermOocytePedigreeDepths.mean();

	// calculate sum of pedigree depths over sperm and oocytes
	pdSpermAndOocyteSum = 0;
	for (int i=0;i<spermPedigreeDepths.size();i++) {
		pdSpermAndOocyteSum += spermPedigreeDepths(i)*4;
	}
	for (int i=0;i<oocytePedigreeDepths.size();i++) {
		pdSpermAndOocyteSum += oocytePedigreeDepths(i);
	}
	pdSpermAndOocyteSum=pdSpermAndOocyteSum/(float)(spermPedigreeDepths.size()*4+oocytePedigreeDepths.size());

	// get average pedigree depth in mitotic region, meiotic region, and whole germline
	pdMitoticRegion = pd_mitotic_region.mean();
	pdMeioticRegion = pd_meiotic_region.mean();
	pdWholeGermline = pd_whole_germline.mean();
	pdGermlineAtOogenesis = pd_atEggLaying.mean();

	// get total number of cell divisions, total number of cells that enter meiotic region
	numCellDivisions=cellDivisionsElapse_;
	numCellsEnterMeioticRegion=countCellsEnterMeioticRegion_;

	// get time when first cell exits the meiotic region
	timeFirstCellExitMeioticRegion = cellProductionTimes(cellExitTimes(0));

	// get time when Nth cell exits the mitotic region
	timeNthCellExitsMitoticRegion = -1;
	for (unsigned int i=0;i<germCellVector_->size();i++) {
		if (minCellsEnterMeioticRegion_>0 && (*germCellVector_)[i]->meioticEntryOrder_==minCellsEnterMeioticRegion_) {
			timeNthCellExitsMitoticRegion = (*germCellVector_)[i]->timeLeftMitoticRegion_;
		}
	}
}


/*!
 * Extract useful data from simulations for calculating fitness, schedules, etc.
 */
void GermLine::calculateRawOutput(Array1D<float> *NCellsEnterMeioticRegionPedigreeDepth,Array1D<float> *oocyteProductionTimes, Array1D<float> *spermProductionTimes, Array1D<float> *oocytePedigreeDepth, Array1D<float> *spermPedigreeDepth, Array1D<float> *oocyteProductionTimes_cellDivision, Array1D<float> *spermProductionTimes_cellDivision, Array1D<float> *mitoticRegionPedigreeDepths, Array1D<float> *meioticRegionPedigreeDepths, Array1D<float> *wholeGermlinePedigreeDepths, Array1D<float> *cellProductionTimes, Array1D<int> *cellExitTimes, Array1D<float> *pedigreeDepthAtEggLaying) {
	for (unsigned int i=0;i<germCellVector_->size();i++) {
		// sperm statistics
		if ((*germCellVector_)[i]->isSperm()) {
			spermPedigreeDepth->push_back((*germCellVector_)[i]->getPedigreeDepth());
			spermProductionTimes->push_back((*germCellVector_)[i]->getTimeGametogenesis_hours());
			spermProductionTimes_cellDivision->push_back((*germCellVector_)[i]->getTimeGametogenesis_cellDivisions());
		}
		// oocyte statistics
		if ((*germCellVector_)[i]->isOocyte()) {
			oocytePedigreeDepth->push_back((*germCellVector_)[i]->getPedigreeDepth());
			oocyteProductionTimes->push_back((*germCellVector_)[i]->getTimeGametogenesis_hours());
			oocyteProductionTimes_cellDivision->push_back((*germCellVector_)[i]->getTimeGametogenesis_cellDivisions());
		}
		// pedigree depth of first N cells that enter mitotic region
		if ((*germCellVector_)[i]->mitoticOrMeioticOrLaid()>0 && (*germCellVector_)[i]->getTimeLeftMitoticRegion_cellDivisions()<=minCellsEnterMeioticRegion_) {
			NCellsEnterMeioticRegionPedigreeDepth->push_back((*germCellVector_)[i]->getPedigreeDepth());
		}
	}

	oocyteProductionTimes_cellDivision->sort();
	spermProductionTimes_cellDivision->sort();
    sort2(spermPedigreeDepth,spermProductionTimes);
    sort2(oocytePedigreeDepth,oocyteProductionTimes);

    // Get pedigree depths of cells in the mitotic region, meiotic region, and whole germline at the end of the simulation
	mitoticRegionPedigreeDepths->resize(0);
	meioticRegionPedigreeDepths->resize(0);
	for (unsigned int i=0;i<germCellVector_->size();i++) {
		int mitoticOrMeioticOrLaid = (*germCellVector_)[i]->mitoticOrMeioticOrLaid();
		float pedigreeDepth = (float)(*germCellVector_)[i]->getPedigreeDepth();
		if (mitoticOrMeioticOrLaid==0) {
			mitoticRegionPedigreeDepths->push_back(pedigreeDepth);
			wholeGermlinePedigreeDepths->push_back(pedigreeDepth);
		} else if (mitoticOrMeioticOrLaid==1) {
			meioticRegionPedigreeDepths->push_back(pedigreeDepth);
			wholeGermlinePedigreeDepths->push_back(pedigreeDepth);
		}
	}

	// Get cell production times (hours)
	cellProductionTimes->resize_unsigned(cellProductionTimes_->size());
	for (unsigned int i=0;i<cellProductionTimes->size();i++) {
		(*cellProductionTimes)(i) = (*cellProductionTimes_)[i];
	}

	// Get time (cell divisions) when cells exit the meiotic region (i.e., undergo oogenesis or apoptosis)
	cellExitTimes->resize_unsigned(exitMeioticRegionTime_->size());
	for (unsigned int i=0;i<cellExitTimes->size();i++) {
		(*cellExitTimes)(i) = (*exitMeioticRegionTime_)[i];
	}

	// Get pedigree depth of the entire germline at each egg laying event
	pedigreeDepthAtEggLaying->resize_unsigned(pdGermlineAtEggLaying->size());
	for (unsigned int i=0;i<pdGermlineAtEggLaying->size();i++) {
		(*pedigreeDepthAtEggLaying)(i) = (*pdGermlineAtEggLaying)[i];
	}

}

/*!
 * Initialize a single cell at the distal tip of the germline.
 * Starts at beginning of G1 (i.e., percentAge = 0)
 */
void GermLine::initialize_oneCell() {
	float meanCellCycleLength = cellCycle_->calculateCellCycleLength(0,0);
	new GermCell(0,0,0,1,meanCellCycleLength,cellCycleNoise_,0,mapPosition_,germCellVector_,false,phaseIndexProfile_,dnaNoise_,pedigreeDepthIncrementMethod_,cellCycle_,&cellDivisionsElapse_,meioticFractionProfile_,generator);
}

/*!
 * Initialize a filled germline with age distribution sampled from exponential distribution.
 * Set all S phase cells as EdU-positive
 */
void GermLine::initialize_filledGermline() {
	int idx=0;
	for (int x=0;x<(int)mapPosition_->size();x++) {
		for (int y=0;y<(int)(*mapPosition_)[(unsigned int)x]->size();y++) {
			float meanCellCycleLength = cellCycle_->calculateCellCycleLength(x,0);
			float percentAge = percentAgeFromExponentialDistribution(generator);
			new GermCell(idx,x,y,percentAge,meanCellCycleLength,cellCycleNoise_,0,mapPosition_,germCellVector_,false,phaseIndexProfile_,dnaNoise_,pedigreeDepthIncrementMethod_,cellCycle_,&cellDivisionsElapse_,meioticFractionProfile_,generator);
			idx++;
		}
	}
}

void GermLine::pulseEdU(){
	for (int x=0;x<(int)mapPosition_->size();x++) {
		for (int y=0;y<(int)(*mapPosition_)[(unsigned int)x]->size();y++) {
			int idx = (*(*mapPosition_)[(unsigned int)x])[(unsigned int)y];
			if (idx!=-1) {
				int phase = (*germCellVector_)[(unsigned int)idx]->getCellCyclePhase();
				if (phase==2) {
					(*germCellVector_)[(unsigned int)idx]->setEdU(true);
				}
			}
		}
	}
}

/*!
 * Calculate index of cell will divide next and time to next division excluding premeiotic cells.
 * Event = 1 if it is a division event, 2 if it is a G2-M transition
 */
void GermLine::calculateNextEvent(unsigned int &cellIndex, float &timeToNextEvent, int &event) {
	// calculate time to next division
	timeToNextEvent = INF;
	int count=0;
	cellIndex=0;
	for (unsigned int i=0; i<germCellVector_->size(); i++) {
		if (((*germCellVector_)[i]->mitoticOrMeioticOrLaid()==0) & ((*germCellVector_)[i]->preMeiotic_==-1 || (*germCellVector_)[i]->preMeiotic_==0)) {
			if ((*germCellVector_)[i]->getTimeToNextDivision() < timeToNextEvent) {
				timeToNextEvent = (*germCellVector_)[(unsigned long)i]->getTimeToNextDivision();
				cellIndex = i;
				event = 1;
				count++;
			}
			if ((*germCellVector_)[i]->getTimeToNextG2toM() < timeToNextEvent && (*germCellVector_)[i]->getTimeToNextG2toM()>0) {
				timeToNextEvent = (*germCellVector_)[(unsigned long)i]->getTimeToNextG2toM();
				cellIndex = i;
				event = 2;
				count++;
			}
		}
	}
	if (count==0) {
		cerr<<"No events"<<endl;
		exit(1);
	}
}

/*!
 * Increments age of all mitotic region cells by dt
 */
void GermLine::incrementAge(float dt) {
	for (unsigned int i=0; i<germCellVector_->size(); i++) {
		if ((*germCellVector_)[i]->mitoticOrMeioticOrLaid()==0) {
			(*germCellVector_)[i]->incrementAge(dt);
		}
	}
}

/*!
 * Write germline state to standard output
 */
void GermLine::display() {
	printf("*************\n");
	printf("division=%d, time=%f, Number of germ cells=%d\n",cellDivisionsElapse_,tElapse_,(int)germCellVector_->size());
	for (unsigned int x=0;x<mapPosition_->size();x++) {
		for (unsigned int y=0;y<(*mapPosition_)[x]->size();y++) {
			int idx=(*(*mapPosition_)[x])[y];
			if (idx!=-1) {
				printf("%3d (%.1f/%.1f) %d,",idx,(*germCellVector_)[(unsigned long)idx]->getAge(),(*germCellVector_)[(unsigned long)idx]->getCellCycleLength(),(*germCellVector_)[(unsigned long)idx]->preMeiotic_);
			} else {
				printf("%3d (%.1f/%.1f) %d,",idx,(float)0,(float)0,0);
			}
		}
		printf("\n");
	}
	printf("*************\n");
}

/*!
 * Write germline state to standard output
 */
void GermLine::displayPedigreeDepth() {
	printf("*************\n");
	printf("division=%d, time=%f, Number of germ cells=%d\n",cellDivisionsElapse_,tElapse_,(int)germCellVector_->size());
	for (unsigned int x=0;x<mapPosition_->size();x++) {
		for (unsigned int y=0;y<(*mapPosition_)[x]->size();y++) {
			int idx=(*(*mapPosition_)[x])[y];
			if (idx!=-1) {
				printf("%3d (%.1f),",idx,(*germCellVector_)[(unsigned long)idx]->getPedigreeDepth());
			} else {
				printf("%3d (  ),",idx);
			}
		}
		printf("\n");
	}
	printf("*************\n");
}

/*!
 * Write germline pre-meiotic state to standard output
 */
void GermLine::displayPreMeioticState() {
	printf("*************\n");
	printf("division=%d, time=%f, Number of germ cells=%d\n",cellDivisionsElapse_,tElapse_,(int)germCellVector_->size());
	for (unsigned int x=0;x<mapPosition_->size();x++) {
		for (unsigned int y=0;y<(*mapPosition_)[x]->size();y++) {
			int idx=(*(*mapPosition_)[x])[y];
			if (idx!=-1) {
				printf("%d, ",(*germCellVector_)[(unsigned long)idx]->preMeiotic_);
			} else {
				printf(" , ");
			}
		}
		printf("\n");
	}
	printf("*************\n");
}

/*!
 * Write cell cycle phase of cells in MR to standard output
 */
void GermLine::displayPhase() {
	printf("*************\n");
	printf("division=%d, time=%f, Number of germ cells=%d\n",cellDivisionsElapse_,tElapse_,(int)germCellVector_->size());
	for (unsigned int x=0;x<mapPosition_->size();x++) {
		for (unsigned int y=0;y<(*mapPosition_)[x]->size();y++) {
			int idx=(*(*mapPosition_)[x])[y];
			if (idx!=-1) {
				printf("%d(%d), ",(*germCellVector_)[(unsigned long)idx]->getCellCyclePhase(),(*germCellVector_)[(unsigned long)idx]->edu_);
			} else {
				printf(" , ");
			}
		}
		printf("\n");
	}
	printf("*************\n");
}

/*!
 * Write germline pre-meiotic fraction to standard output
 */
void GermLine::displayPreMeioticFraction() {
	printf("*************\n");
	printf("division=%d, time=%f, Number of germ cells=%d\n",cellDivisionsElapse_,tElapse_,(int)germCellVector_->size());
	for (unsigned int x=0;x<mapPosition_->size();x++) {
		float sum=0;
		for (unsigned int y=0;y<(*mapPosition_)[x]->size();y++) {
			int idx=(*(*mapPosition_)[x])[y];
			if (idx!=-1 && (*germCellVector_)[(unsigned long)idx]->preMeiotic_==1) {
				sum++;
			}
		}
		printf("%.2f,",sum/(float)(*mapPosition_)[x]->size());
	}
	printf("\n*************\n");
}

/*!
 * Get spatial cell cycle phase indices
 */
void GermLine::getSpatialCellCyclePhaseIndices(Array1D<float> *G1idx,Array1D<float> *Sidx,Array1D<float> *G2idx,Array1D<float> *Midx) {
	G1idx->resize(mapPosition_->size());
	G1idx->fill(0.0f);
	Sidx->resize(mapPosition_->size());
	Sidx->fill(0.0f);
	G2idx->resize(mapPosition_->size());
	G2idx->fill(0.0f);
	Midx->resize(mapPosition_->size());
	Midx->fill(0.0f);
	for (unsigned int x=0;x<mapPosition_->size();x++) {
		for (unsigned int y=0;y<(*mapPosition_)[x]->size();y++) {
			int idx=(*(*mapPosition_)[x])[y];
			if (idx!=-1) {
				int phase =  (*germCellVector_)[(unsigned long)idx]->getCellCyclePhase();
				if (phase==1) {
					(*G1idx)(x)=(*G1idx)(x)+1;
				} else if (phase==2) {
					(*Sidx)(x)=(*Sidx)(x)+1;
				} else if (phase==3) {
					(*G2idx)(x)=(*G2idx)(x)+1;
				} else if (phase==4) {
					(*Midx)(x)=(*Midx)(x)+1;
				}
			}
		}
		(*G1idx)(x) = (*G1idx)(x) / (float)((*mapPosition_)[x]->size());
		(*Sidx)(x) =  (*Sidx)(x) / (float)((*mapPosition_)[x]->size());
		(*G2idx)(x) = (*G2idx)(x) / (float)((*mapPosition_)[x]->size());
		(*Midx)(x) =  (*Midx)(x) / (float)((*mapPosition_)[x]->size());
	}
}

/*!
 * 	Get what simulated M-phase fraction is
 *
 */
float GermLine::getFractionMphase(int x) {
	unsigned int x0=(unsigned int)max(x-M_SEARCH_RANGE,0);
	unsigned int x1=(unsigned int)min(x+M_SEARCH_RANGE,(int)mapPosition_->size()-1);
	float total=0;
	float numMphase=0;
	for (unsigned int i=x0;i<=x1;i++) {
		for (unsigned int y=0;y<(*mapPosition_)[i]->size();y++) {
			int idx = (*(*mapPosition_)[i])[y];
			if (idx!=-1) {
				total++;
			}
			if (idx!=-1 && (*germCellVector_)[(unsigned int)idx]->getCellCyclePhase()==4) {
				numMphase++;
			}
		}
	}
	float simulatedFractionM = (numMphase-1)/total; // numMphase-1 because we don't count cell that just entered G2
	return simulatedFractionM;
}

/*!
 * Returns 1 if a cell division event, returns 2 if a G2-M transition
 */
int GermLine::step(float &time_nextEvent) {
	// calculate time and cell index of next dividing cell
	unsigned int idx_nextEvent=0;
	int event=-1;
	calculateNextEvent(idx_nextEvent,time_nextEvent,event);

	// increment cell ages
	incrementAge(time_nextEvent);

	// G2-M transition event
	if (event==2) {
		int x=(*germCellVector_)[idx_nextEvent]->getX();
		float theoreticalFractionM = meioticFractionProfile_->calculateMphaseThreshold(x,cellDivisionsElapse_);
		float simulatedFraction = getFractionMphase(x);
		if (simulatedFraction<theoreticalFractionM) {
			(*germCellVector_)[idx_nextEvent]->preMeiotic_=0;
		}
		else {
			(*germCellVector_)[idx_nextEvent]->preMeiotic_=1;
		}
	// Division event
	} else if (event==1) {
		// mother divides
		int daughter = (*germCellVector_)[idx_nextEvent]->divide();
		// push daughter
		moveCellsWrapper((unsigned int)daughter);
		// restructure geometry
		restructureGeometry();
		// rescale cell cycle lengths
		rescaleCellCycleLengths();
	} else {
		cerr<<"Invalid event"<<endl;
		exit(1);
	}
	return event;
}

/*!
 * Run germline simulation
 */
void GermLine::run(string initialCondition, int preRunSteps) {
	// initialize mitotic region
	if (initialCondition.compare("oneCell")==0) {
		initialize_oneCell();
	} else if (initialCondition.compare("filled")==0) {
		initialize_filledGermline();
	} else {
		cerr << "Invalid initial condition" << endl;
		exit(1);
	}

	// pre-run simulation for a set amount of cell divisions
	for (int i=0;i<preRunSteps;i++) {
		float time_nextEvent;
		step(time_nextEvent);
	}

	// pulse S phase cells with EdU
	pulseEdU();

	// Display debugging output
	if (verbose_) {
		cout << "********************" << endl;
		cellCycle_->display();
		cout << "********************" << endl;
		geometry_->display();
		cout << "********************" << endl;
		meioticFractionProfile_->display();
		cout << "********************" << endl;
		phaseIndexProfile_->display();
		cout << "********************" << endl;
	}

	while (1) {
		// Display debugging output
		if (verbose_) {
			//display();
			displayPedigreeDepth();
			//displayPreMeioticFraction();
			//displayPhase();
			//displayPreMeioticState();
			cin.ignore(1);
		}

		float time_nextEvent;
		int event = step(time_nextEvent);

		// increment timers
		if (event==1) {
			cellDivisionsElapse_++;
			cellProductionTimes_->push_back(tElapse_);
		}
		tElapse_+=time_nextEvent;


		// check stopping condition
		bool sufficientOocytesProduced        = countOocyteLaid_>=minNumOocytesProduced_;
		bool sufficientCellEnterMeioticRegion = countCellsEnterMeioticRegion_>= minCellsEnterMeioticRegion_;
		bool sufficientTimeElapsed            = tElapse_>=minTimeElapsed_;
		bool sufficientCellDivisions          = cellDivisionsElapse_>=minCellDivisions_;
		if (sufficientOocytesProduced && sufficientTimeElapsed && sufficientCellDivisions && sufficientCellEnterMeioticRegion) {
			break;
		}
	}
}

/*!
 * Calculate
 */
float GermLine::calculatePercentPremeiotic() {
	float preMeiotic=0;
	float mitotic=0;
	for (unsigned int i=0; i<germCellVector_->size(); i++) {
		if ((*germCellVector_)[i]->mitoticOrMeioticOrLaid()==0) {
			if ((*germCellVector_)[i]->preMeiotic_==1) {
				preMeiotic++;
			} else if ((*germCellVector_)[i]->preMeiotic_==0 || (*germCellVector_)[i]->preMeiotic_==-1) {
				mitotic++;
			}
		}
	}
	return preMeiotic/(mitotic+preMeiotic);
}


/*!
 * Rescale cell cycle lengths to account for cell cycle profile changing over time
 */
void GermLine::rescaleCellCycleLengths() {
	// only rescale every CCRESCALEINTERVAL divisions
	if (cellDivisionsElapse_ % CCRESCALEINTERVAL != 0) {
		return;
	}
	// don't need to rescale if there is only 1 temporal control point (i.e., cell cycle profile constant in time)
	if (cellCycle_->profileConstant(cellDivisionsElapse_,cellDivisionsElapse_+1)) {
		return;
	}
	for (unsigned int i=0;i<germCellVector_->size();i++) {
		if ((*germCellVector_)[i]->mitoticOrMeioticOrLaid()==0) {
			int xPosition = (*germCellVector_)[i]->getX();
			float newMeanCellCycleLength = cellCycle_->calculateCellCycleLength(xPosition,cellDivisionsElapse_+1);
			(*germCellVector_)[i]->rescaleCellCycleLength(newMeanCellCycleLength);
		}
	}
}

/*!
 * Push a cell from the mitotic region into the meiotic region
 * Simulate all subsequent events (spermatogenesis/oogenesis/apoptosis)
 */
void GermLine::pushCellFromMR(unsigned int index) {
	// Push cell out of mitotic region into meiotic region.
	countCellsEnterMeioticRegion_++;
	(*germCellVector_)[index]->pushCellFromMitoticRegion(tElapse_,countCellsEnterMeioticRegion_);
	meioticGermCellVector_->push_back((int)index);

	// cell becomes sperm as soon as it enters meiotic region if not enough sperm are produced
	if (countSperm_<numSpermToProduce_) {
		countSperm_++;
		(*germCellVector_)[index]->spermatogenesis(tElapse_,(float)cellDivisionsElapse_);
		meioticGermCellVector_->pop_back();
		return;
	}

	// Apoptosis/oogenesis until meiotic region meets size constraints
	int meioticSize = (int)meioticGermCellVector_->size();
	int meioticSize_constraint = meioticRegionSize_->calculateInterpolatedValue(0,cellDivisionsElapse_);
	int numCellsToRemove=max(meioticSize-meioticSize_constraint,0);
	for (int i=0;i<numCellsToRemove;i++) {
		// increment timer
		cellsExitMeioticRegion_++;
		// get index of oldest germ cell in the meiotic region
		int idx = (*meioticGermCellVector_)[0];

		// calculate probability of oogenesis (vs. apoptosis)
		float probability_oocyte;
		if (strcmp(oocyteProbabilityMethod_.c_str(),"linear")==0) {
			probability_oocyte = oocyteProbability_->calculateInterpolatedValue(0,cellsExitMeioticRegion_);
		} else if (strcmp(oocyteProbabilityMethod_.c_str(),"step")==0) {
			probability_oocyte = oocyteProbability_->calculateFloorValue(0,cellsExitMeioticRegion_);
		} else {
			cerr << "Invalid method choice for oocyteProbability" << endl;
			exit(1);
		}

		// simulate oogenesis/apoptosis
		float rnd=randUniformFloat(0,1,generator);
		if (rnd<=probability_oocyte && countOocyteLaid_<minNumOocytesProduced_) {  // oogenesis
			(*germCellVector_)[(unsigned long)idx]->oogenesis(tElapse_,(float)cellDivisionsElapse_);
			countOocyteLaid_++;
			// record average pedigree depth of germline at time of egg laying
			float count=0,sum=0;
			for (unsigned int j=0;j<germCellVector_->size();j++) {
				if ((*germCellVector_)[j]->mitoticOrMeioticOrLaid()==0 || (*germCellVector_)[j]->mitoticOrMeioticOrLaid()==1) {
					sum+=(float)(*germCellVector_)[j]->getPedigreeDepth();
					count++;
				}
			}
			pdGermlineAtEggLaying->push_back(sum/count);
		} else {	// apoptosis
			(*germCellVector_)[(unsigned long)idx]->die();
		}

		// remove cell from meiotic region
		meioticGermCellVector_->erase(meioticGermCellVector_->begin());
		exitMeioticRegionTime_->push_back(cellDivisionsElapse_);
	}
}

/*!
 * This is a wrapper function for moveCells to handle the "golden strand" method of moving cells (i.e., daughter cells stay in cell row 1)
 */
void GermLine::moveCellsWrapper(unsigned int pushingCellIdx) {
	if (pedigreeDepthIncrementMethod_ == 3 && (*germCellVector_)[pushingCellIdx]->getX()==0) {  	// golden strand simulations in row 1
		// check if row is filled
		bool rowFilled = true;
		unsigned long yDim = (*mapPosition_)[0]->size();
		for (unsigned int y=0;y<yDim;y++) {
			if ((*(*mapPosition_)[0])[y] == -1) {
				rowFilled = false;
			}
		}
		// move mother cells forward if row is filled
		if (rowFilled) {
			moveCells(pushingCellIdx,0,1);
		} else {
			moveCells(pushingCellIdx,0,0);
		}
	} else {   	// in non-golden strand simulations, cells in row 1 move like any other cells
		moveCells(pushingCellIdx,0,0);
	}
}

/*
 * Move cells so that no cells occupy the same space
 * Forward movement is only allowed when a cell row is completely filled
 * direction_UpDown is either negative (down within the same row) or positive (up within the same row) or zero (no movement within the same row, must move forward)
 */
void GermLine::moveCells(unsigned int pushingCellIdx, int direction_UpDown, int direction_Forward) {
	// Sanity check
	#if defined INCLUDESANITYCHECKS
	if (direction_Forward==1 && direction_UpDown!=0) {
		cerr << "Cannot move cells horizontally and vertically at the same time" << endl;
		exit(1);
	}
	#endif

	// useful quantities
	int current_X = (*germCellVector_)[pushingCellIdx]->getX();
	int current_Y = (*germCellVector_)[pushingCellIdx]->getY();
	int current_rowSize = (int)(*mapPosition_)[(unsigned int)current_X]->size();

	 // If cell is in empty space, then settle in
	if ((*(*mapPosition_)[(unsigned int)current_X])[(unsigned int)current_Y]==-1) {
		(*(*mapPosition_)[(unsigned int)current_X])[(unsigned int)current_Y]=(int)pushingCellIdx;
		return;
	}

	// simulate a cell moving forward
	if (direction_Forward==1) {
		// If cell is in last row then push cell out of MR
		int lastRow = (int)mapPosition_->size();
		if (current_X==lastRow-1) {
			pushCellFromMR(pushingCellIdx);
			return;
		}
		// If next row has 0 cells then push cell out of MR
		if ((*mapPosition_)[(unsigned int)current_X+1]->size()==0) {
			pushCellFromMR(pushingCellIdx);
			return;
		}
		// get position of nearest space in next row assuming a circular geometry
		int new_X = current_X+1;
		int new_Y;
		if (current_rowSize==1) {
			new_Y=0;
		} else {
			new_Y = current_Y*((int)(*mapPosition_)[(unsigned int)current_X+1]->size()-1)/(current_rowSize-1);
		}
		// update coordinates of cell being pushed
		int cellBeingPushedIdx = (*germCellVector_)[pushingCellIdx]->setPosition(new_X,new_Y,cellCycle_->calculateCellCycleLength(new_X,cellDivisionsElapse_));
		// move pushed cell
		if (cellBeingPushedIdx==-1) { 		// if pushed into an empty space, we're done
			return;
		}
		moveCells((unsigned int)cellBeingPushedIdx,0,0);
		return;
	}

	// push cell up or down in same row
	if (direction_UpDown!=0) {
		int new_X=current_X;
		int new_Y,new_direction_UpDown;
		// cell gets pushed down
		if (direction_UpDown<0) {
			new_direction_UpDown = direction_UpDown+1;
			new_Y = current_Y-1;
			if (new_Y<0) {
				new_Y += current_rowSize;
			}
		}
		// cell gets pushed up
		else {
			new_direction_UpDown = direction_UpDown-1;
			new_Y = current_Y+1;
			if (new_Y>=current_rowSize) {
				new_Y -= current_rowSize;
			}
		}
		#if defined INCLUDESANITYCHECKS
		// sanity check
		if (new_Y<0 || new_Y>=current_rowSize) {
			cerr << "Problem with moveCell, new coordinates not within y-range" << endl;
			exit(1);
		}
		#endif

		// update coordinates of cell being pushed
		int cellBeingPushedIdx = (*germCellVector_)[pushingCellIdx]->setPosition(new_X,new_Y,cellCycle_->calculateCellCycleLength(new_X,cellDivisionsElapse_));

		// move pushed cell
		if (cellBeingPushedIdx==-1) { 		// if pushed into an empty space, we're done
			return;
		}
		if (new_direction_UpDown==0) {
			moveCells((unsigned int)cellBeingPushedIdx,0,1);
		} else {
			moveCells((unsigned int)cellBeingPushedIdx,new_direction_UpDown,0);
		}
		return;
	}

	// If no movement direction is specified, then randomly whether to push cell forward or within same row
	if (direction_Forward==0 && direction_UpDown==0) {
		// count the number of empty spaces in the same row
		int countEmptyCells = 0;
		for (int y=0;y<current_rowSize;y++) {
			if ((*(*mapPosition_)[(unsigned int)current_X])[(unsigned int)y]==-1){
				countEmptyCells++;
			}
		}

		// if there are empty spaces in the same row, pick one at random and shift cells to fill in that direction
		if (countEmptyCells>0) {

			// y-position to push to excluding non-empty spaces
			int y_excludeNonempty =  (int)floorf(randUniformFloat(0,1,generator)*(float)countEmptyCells);

			// get y-position to push to including non-empty spaces
			int counter = 0;
			int y_includeNonempty=0;
			for (int y=0;y<current_rowSize;y++) {
				if ((*(*mapPosition_)[(unsigned int)current_X])[(unsigned int)y]==-1) {
					if (counter==y_excludeNonempty) {
						y_includeNonempty=y;
						break;
					} else {
						counter++;
					}
				}
			}

			// push cells toward chosen y-position
			int y1=current_Y;  // this is the y-position of the current cell
			int y2=y_includeNonempty;	// this is the y-position we are pushing toward
			int shortest_distance = minCircularDistance(y1,y2,current_rowSize,generator); // this is the shortest displacement from y1 to y2
			moveCells(pushingCellIdx,shortest_distance,0);
		}

		// if there are no empty cells in the same row, randomly shift within same row or forward
		else {

			// shift forward if only one cell in row, otherwise shift forward with random probability
			float rnd = randUniformFloat(0,1,generator);
			if (rnd<=movementProbabilityForward_ || current_rowSize==1) {
				moveCells(pushingCellIdx,0,1);
			}
			// shift to random space in same row with random probability
			else {
				// pick random cell in same row, making sure not to pick same cell
				int y1 = current_Y;
				int y2=y1;
				while (y2==y1) {
					y2 = (int)boost::math::iround(randUniformFloat(0,1,generator)*(float)(current_rowSize-1));
				}
				int shortest_distance = minCircularDistance(y1,y2,current_rowSize,generator);
				moveCells(pushingCellIdx,shortest_distance,0);
			}
		}
		return;
	}
	// should never reach this point
}

/*!
 * Change germline geometry, pushing cells appropriately
 */
void GermLine::restructureGeometry() {
	// if no change in geometry, then don't need to do anything
	if (geometry_->geometryConstant(cellDivisionsElapse_,cellDivisionsElapse_+1)) {
		return;
	}

	// get current and future geometry
	Array1D<int> cellPerRow_t1,cellPerRow_t2;
	geometry_->getCurrentGeometry(cellDivisionsElapse_,&cellPerRow_t1);
	geometry_->getCurrentGeometry(cellDivisionsElapse_+1,&cellPerRow_t2);

	// grow rows
	for (int r=0; r<cellPerRow_t1.size(); r++) {
		for (int i=cellPerRow_t1(r);i<cellPerRow_t2(r);i++) {
			(*mapPosition_)[(unsigned int)r]->push_back(-1);
		}
	}

	// shrink rows, push cells that are squeezed out
	for (unsigned int r=0; r<(unsigned int)cellPerRow_t1.size(); r++) {  // shrink germLine and push cells that are squeezed out
		for (int i=cellPerRow_t1(r);i>cellPerRow_t2(r);i--) {
			// get index of cell that will be squeezed out
			int pushingCellIdx = (*mapPosition_)[r]->back();
			// if there is no cell there, don't need to do any cell movement
			if (pushingCellIdx==-1) {
				(*mapPosition_)[r]->pop_back(); // alter germline geometry, don't need to worry about moving cells since there is no cell there
			}
			// if cell row is about to disappear, then cell will be squeezed out of mitotic region
			else if ((*mapPosition_)[r]->size()==1) {
				pushCellFromMR((unsigned int)pushingCellIdx);
				(*mapPosition_)[r]->pop_back();
			}
			// squeeze cell upward
			else {
				int newY=(*germCellVector_)[(unsigned int)pushingCellIdx]->getY()-1;  // new y-position after cell is pushed upwards
				int cellBeingPushed = (*germCellVector_)[(unsigned int)pushingCellIdx]->setPosition((int)r,newY,cellCycle_->calculateCellCycleLength((int)r,cellDivisionsElapse_));
				(*mapPosition_)[r]->pop_back(); // alter geometry
				if (cellBeingPushed!=-1) {
					moveCells((unsigned int)cellBeingPushed,0,0); // move cells
				}
			}
		}
	}
}
