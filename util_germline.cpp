/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#include "util_germline.h"
using namespace std;

/*!
 * print cell cycle phase indices
 */
void printCellPhaseIndices(Array2D<float> *G1idx, Array2D<float> *Sidx, Array2D<float> *G2idx, Array2D<float> *Midx) {
	int dimx,dimy;
	G1idx->getDimensions(dimx,dimy);
	printf("G1= [");
	for (int i=0;i<dimx;i++) {
		printf("%.3f",G1idx->mean_x(i));
		if (i==dimx-1) {
			printf("]\n");
		} else {
			printf(", ");
		}
	}
	printf("S = [");
	for (int i=0;i<dimx;i++) {
		printf("%.3f",Sidx->mean_x(i));
		if (i==dimx-1) {
			printf("]\n");
		} else {
			printf(", ");
		}
	}
	printf("G2= [");
	for (int i=0;i<dimx;i++) {
		printf("%.3f",G2idx->mean_x(i));
		if (i==dimx-1) {
			printf("]\n");
		} else {
			printf(", ");
		}
	}
	printf("M = [");
	for (int i=0;i<dimx;i++) {
		printf("%.3f",Midx->mean_x(i));
		if (i==dimx-1) {
			printf("]\n");
		} else {
			printf(", ");
		}
	}
}

/*!
 * Returns true if the maximum size of the mitotic region (over time) exceeds threshold
 */
bool mitoticRegionTooBig(ProgramOptionParser *opt) {
	GeometryProfile gP;
	opt->getGeometryProfile(&gP);
	int maxSizeOfMR = gP.maxSizeOfMR();
	if (maxSizeOfMR>opt->maxMitoticRegionSize) {
		return true;
	} else {
		return false;
	}
}

/*!
 * Sort two vectors simultaneously based on v2
 */
void sort2(Array1D<float> *v1, Array1D<float> *v2) {
	#if defined INCLUDESANITYCHECKS
	// Sanity check
	int sz1 = v1->size();
	int sz2 = v2->size();
	if (sz1!=sz2) {
		printf("Cannot sort vectors that don't have same size (sz1=%d, sz2=%d)",sz1,sz2);
		exit(1);
	}
	#endif

	// sort v2 and get sorted indices
	Array1D<int> sortedIdx;
	v2->sort(&sortedIdx);

	// sort v1 based on v2
	Array1D<float> v1Copy(*v1);
	for (int i=0;i<v2->size();i++) {
		(*v1)(i) = v1Copy(sortedIdx(i));
	}
}


/*!
 * Calculate the doubling time for a given intrinsic growth rate r
 */
float calculateDoublingTime(float r) {
	return logf(2)/logf(r);
}

/*!
 * Calculate generation fitness using rigorous selection coefficient approach
 */
float calculateGenerationFitness_selectionCoefficient(Array2D<float> *simulatedFecunditySchedule,Array2D<float> *fecundityConstraints, float pd, float deleteriousMutationRate, float maximalGenerationRate, float maximalPedigreeDepth) {
        float rm = calculateGenerationRate(simulatedFecunditySchedule,fecundityConstraints);  // generation rate in mutant worm excluding effects of pedigree depth reduction
        float mq = logf(rm);  // calculate malthusian parameter corresponding to generation rate

        float r0 = maximalGenerationRate;                    // generation rate in a wild-type worm
        float mr = logf(r0);

        float sgr = (mq - mr)/mr;

        int numProgeny = 0;
        int dimx,dimy;
        simulatedFecunditySchedule->getDimensions(dimx,dimy);
        for (int i=1;i<dimy;i++) {
             numProgeny+=(*simulatedFecunditySchedule)(0,i);
        }

        float s = expf( sgr * logf((float)numProgeny)) - 1;  // selection coefficient WITHOUT TAKING PEDIGREE DEPTH REDUCTION INTO ACCOUNT

        float U = deleteriousMutationRate; // number of deleterious mutations in a worm with flat cell cycle profile
        float p0 = maximalPedigreeDepth;     // average pedigree depth in a flat profile worm
        float p = pd;                    // average pedigree depth in an mutant worm
        float dU = U*(p-p0)/p0;          // change in the average number of deleterious mutations.  Positive if mutant worm has more mutations

        float s2 = -dU/2;

        float generationRate = s + s2;
        return generationRate;
}

/*!
 * Calculate the fitness metric based on generation rate.
 * Generation rate is affected by two processes:
 *   1) Fecundity schedule
 *   2) Number of mutations
 */
float calculateGenerationFitness(Array2D<float> *simulatedFecunditySchedule,Array2D<float> *fecundityConstraints, float pd, float deleteriousMutationRate, float maximalGenerationRate, float maximalPedigreeDepth) {
        float rm = calculateGenerationRate(simulatedFecunditySchedule,fecundityConstraints);  // generation rate in mutant worm excluding effects of pedigree depth reduction
        float r0 = maximalGenerationRate;                    // generation rate in a worm with flat cell cycle profile
        float U = deleteriousMutationRate; // number of deleterious mutations in a worm with flat cell cycle profile
        float p0 = maximalPedigreeDepth;     // average pedigree depth in a worm with flat cell cycle profile
        float p = pd;                                    // average pedigree depth in an mutant worm
        float dU = U*(p-p0)/p0;          // change in the average number of deleterious mutations.  Positive if mutant worm has more mutations

        r0 = rm;
        float generationRate = rm - dU/2*r0;
        return generationRate;
}

/*!
 * Calculate the hourly generation rate at steady based on simulated fecundity schedule
 */
float calculateGenerationRate(Array2D<float> *simulatedFecunditySchedule,Array2D<float> *fecundityConstraints) {

	// calculate transition matrix m
	Array2D<float>::array_index dummy,dimy_fecundity;
	simulatedFecunditySchedule->getDimensions(dummy,dimy_fecundity);
	int dim_L2start = (int)ceilf((*fecundityConstraints)(0,dimy_fecundity-1))+1;
        int dim = dim_L2start + TIME_EGGLAY_L2;
	Eigen::MatrixXd m(dim,dim);
	// initialize m to zeros
	for (int i=0;i<dim;i++) {
		for (int j=0;j<dim;j++) {
			m(i,j)=0;
		}
	}

	// fill first row of matrix
	for (int i=0;i<dimy_fecundity-1;i++) {
		float numOocytes = simulatedFecunditySchedule->mean(i+1);
		int idxLow = calculateSmallestIntegerGreaterThan((*fecundityConstraints)(0,i)) + TIME_EGGLAY_L2;
		int idxHigh = calculateLargestIntegerLessEqualThan((*fecundityConstraints)(0,i+1)) + TIME_EGGLAY_L2;
		float divisor = (float)(idxHigh-idxLow+1);
		for (int j=idxLow;j<=idxHigh;j++) {
			m(0,j) = numOocytes/divisor;
		}
	}

	// fill rest of matrix
	for (int i=0;i<dim-1;i++) {
		m(i+1,i) = 1;
	}

	// calculate eigenvalues of transition matrix
	Eigen::VectorXcd eivals = m.eigenvalues();

	// calculate dominant eigenvalue
	float dominantEigenvalue = -INF;
	for (int i=0;i<dim;i++) {
		float a = (float)eivals(i).real();
		float b = (float)eivals(i).imag();
		if (b<EPS) {
			if (a>dominantEigenvalue) {
				dominantEigenvalue=a;
			}
		}
	}

	return dominantEigenvalue;
}

/*!
 * Print simulation results to standard output
 */
void printSimulationResults(PedigreeDepthStorage *pdStorage, int numSperm, bool cellProductionConstraintsMet, bool fecundityScheduleConstraintsMet, bool allConstraintsMet, ProgramOptionParser *opt) {
	// print fitness metric values
	printf("\n******************************\n");
	printf("Fitness metrics:\n");
	printf("pedigree depth sperm (mean) = %f\n",pdStorage->pd_sperm->mean());
	printf("pedigree depth oocyte (mean) = %f\n",pdStorage->pd_oocyte->mean());
	printf("pedigree depth sperm and oocyte (mean) = %f\n",pdStorage->pd_spermOocyte->mean());
	printf("pedigree depth mitotic region (mean) = %f\n",pdStorage->pd_mitoticRegion->mean());
	printf("pedigree depth meiotic region (mean) = %f\n",pdStorage->pd_meioticRegion->mean());
	printf("pedigree depth whole germline (mean) = %f\n",pdStorage->pd_wholeGermline->mean());
	printf("pedigree depth sperm (median) = %f\n",pdStorage->pd_sperm->perctile(0.5));
	printf("pedigree depth oocyte (median) = %f\n",pdStorage->pd_oocyte->perctile(0.5));
	printf("pedigree depth sperm and oocyte (median) = %f\n",pdStorage->pd_spermOocyte->perctile(0.5));
	printf("pedigree depth mitotic region (median) = %f\n",pdStorage->pd_mitoticRegion->perctile(0.5));
	printf("pedigree depth meiotic region (median) = %f\n",pdStorage->pd_meioticRegion->perctile(0.5));
	printf("pedigree depth whole germline (median) = %f\n",pdStorage->pd_wholeGermline->perctile(0.5));
	printf("pedigree depth sperm and oocyte (sum/numGametes) = %f\n",pdStorage->pd_spermOocyteSum->mean());
	printf("pedigree depth of germline at egg laying = %f\n",pdStorage->pd_GermlineAtEggLaying->mean());
	printf("pedigree depth of first N cells that entered meiotic region = %f\n",pdStorage->pd_AllCellsEnterMeioticRegion->mean());
	// print generation rates
	float r = calculateGenerationRate(pdStorage->simulatedFecunditySchedule,pdStorage->fecundityConstraints);
	float rFitness = calculateGenerationFitness(pdStorage->simulatedFecunditySchedule,pdStorage->fecundityConstraints,pdStorage->pd_spermOocyte->mean(),opt->deleteriousMutationRate,opt->maximalGenerationRate,opt->maximalPedigreeDepth);
	printf("generation rate without pedigree depth contribution (hours^-1) = %f\n",r);
	printf("negative generation rate with pedigree depth contribution (hours^-1, negative so minimization will optimize fitness) = %f\n",-rFitness); // negative so that minimization will optimize fitness
	printf("\n");

	// print timers
	printf("Timers:\n");
	printf("Number of cell divisions = %f\n",pdStorage->num_cellDivisions->mean());
	printf("Number of cells that enter meiotic region = %f\n",pdStorage->num_meioticCellsProduced->mean());
	printf("Time of first cell exiting meiotic region = %f\n",pdStorage->time_firstMeioticExit->mean());
	printf("\n");

	// print simulated fecundity schedule
	Array2D<float>::array_index dummy,dimy_fecundity;
	pdStorage->fecundityConstraints->getDimensions(dummy,dimy_fecundity);
	printf("Simulated fecundity schedule:\n");
	for (int i=0;i<dimy_fecundity;i++) {
		float meanGametesProduced = pdStorage->simulatedFecunditySchedule->mean(i);
		if (i==0) {
			printf("%f sperm produced at %f hours (%d expected)\n",meanGametesProduced,(*pdStorage->fecundityConstraints)(0,i),numSperm);
		} else {
			printf("%f oocytes produced at %f hours (%d expected)\n",meanGametesProduced,(*pdStorage->fecundityConstraints)(0,i),(int)(*pdStorage->fecundityConstraints)(i));
		}
	}
	printf("\n");

	// print simulated cell exiting meiotic region schedule
	printf("Simulated cells that hit end of meiotic region (and undergo apoptosis/oogenesis):\n");
	for (int i=0;i<dimy_fecundity;i++) {
		float meanCellDivisions = pdStorage->simulatedCellProductionSchedule_atFecundityBins->mean(i);
		float meanNumCellsHitMeioticRegion = pdStorage->simulatedMeioticExitSchedule->mean(i);
		printf("%f cell divisions, %f cells hit end of meiotic region at at %f hours\n",meanCellDivisions,meanNumCellsHitMeioticRegion,(*pdStorage->fecundityConstraints)(0,i));
	}
	printf("\n");

	printf("theoretical oocyte probability = step*0:");
	float total_numCellsHitMeioticRegion=0;
	for (int i=1;i<dimy_fecundity;i++) {
		float numCellsHitMeioticRegion       =  pdStorage->simulatedMeioticExitSchedule->mean(i);
		float expectedNumberOocytes = (*pdStorage->fecundityConstraints)(i);
		printf("%.3f",expectedNumberOocytes/numCellsHitMeioticRegion);
		if (i!=dimy_fecundity-1) {
			total_numCellsHitMeioticRegion += numCellsHitMeioticRegion;
			printf("/%.2f:",total_numCellsHitMeioticRegion);
		}
	}
	printf("\n\n");

	// print simulated cell production schedule
	printf("cell production schedule:\n");
	Array2D<float>::array_index dimx_cellProduction,dimy_cellProduction;
	pdStorage->simulatedCellProductionSchedule->getDimensions(dimx_cellProduction,dimy_cellProduction);
	for (int i=0;i<dimy_cellProduction;i++) {
		float time = pdStorage->simulatedCellProductionSchedule->mean(i);
		float constraint = (*pdStorage->cellProductionConstraints)(0,i);
		int cellDivisions = (int)(*pdStorage->cellProductionConstraints)(i);
		printf("%d cell divisions at %f hours (constraint=%f)\n",cellDivisions,time,constraint);
	}
	printf("\n");

	// print whether constraints were met
	printf("cell production constraints met: %d\n", cellProductionConstraintsMet);
	printf("fecundity constraints met: %d\n", fecundityScheduleConstraintsMet);
	printf("all constraints met: %d\n", allConstraintsMet);

	// pedigree depth for individual runs
	printf("\n****** Individual Output ******\n");
	printf("sperm/oocyte mean pedigree depth = ");
	pdStorage->pd_spermOocyte->display();

	// print whether cell production constraint is met for individual runs
	printf("cell production constraint = ");
	for (int x=0;x<dimx_cellProduction;x++) {
		bool constraintsSatisfied=true;
		for (int y=0;y<dimy_cellProduction;y++) {
			if ((*pdStorage->simulatedCellProductionSchedule)(x,y)>(*pdStorage->cellProductionConstraints)(0,y)) {
				constraintsSatisfied=false;
				break;
			}
		}
		if (opt->constraintMethod.compare("none")==0){
			constraintsSatisfied=true;
		}
		printf("%d",constraintsSatisfied);
		if (x!=dimx_cellProduction-1) {
			printf(",");
		}
	}
	printf("\n");

	// print whether fecundity constraint is met for individual runs
	Array2D<float>::array_index dimx_simulatedFecundity, dimy_simulatedFecundity;
	printf("fecundity constraint = ");
	pdStorage->simulatedFecunditySchedule->getDimensions(dimx_simulatedFecundity,dimy_simulatedFecundity);
	for (Array2D<float>::array_index x=0;x<dimx_simulatedFecundity;x++) {
		bool constraintsSatisfied=true;
		int sum_sim = 0;
		int sum_sched = 0;
		for (int y=0;y<dimy_simulatedFecundity;y++) {
			if (y==0) {  			// check sperm constraint
				int numSperm_sim = boost::math::iround((*pdStorage->simulatedFecunditySchedule)(x,y));
				int numSperm_sched = opt->numSperm; // boost::math::iround((*pdStorage->fecundityConstraints)(0,y));
				if (numSperm_sim!=numSperm_sched) {
					constraintsSatisfied=false;
					break;
				}
			} else {				// check oocyte constraint
				sum_sim+=ceil((*pdStorage->simulatedFecunditySchedule)(x,y));
				sum_sched+=boost::math::iround((*pdStorage->fecundityConstraints)(y));
				if (sum_sim<sum_sched) {
					constraintsSatisfied=false;
					break;
				}
			}
		}
		if (opt->constraintMethod.compare("none")==0){
			constraintsSatisfied=true;
		}
		printf("%d",constraintsSatisfied);
		if (x!=dimx_simulatedFecundity-1) {
			printf(",");
		}
	}
	printf("\n");

	for (int i=0;i<dimy_cellProduction;i++) {
		printf("cell production time by simulation (constraint=%f) = ",(*pdStorage->cellProductionConstraints)(0,i));
		for (int j=0;j<dimx_cellProduction;j++) {
			printf("%f",(*pdStorage->simulatedCellProductionSchedule)(j,i));
			if (j!=dimx_cellProduction-1) {
				printf(",");
			}
		}
		printf("\n");
	}

	int sumConstraint=0;
	Array1D<float> sumGamete;
	sumGamete.resize(dimx_cellProduction);
	for (int i=0;i<dimy_simulatedFecundity;i++) {
		int constraint;
		if (i==0) {
			constraint = opt->numSperm;
		} else {
			constraint = boost::math::iround((*pdStorage->fecundityConstraints)(i));
		}
		sumConstraint+=constraint;
		printf("fecundity time by simulation (constraint=%d) = ",sumConstraint);
		for (int j=0;j<dimx_cellProduction;j++) {
			sumGamete(j)+=(*pdStorage->simulatedFecunditySchedule)(j,i);
			printf("%d",boost::math::iround(sumGamete(j)));
			if (j!=dimx_cellProduction-1) {
				printf(",");
			}
		}
		printf("\n");
	}
}

/*!
 *  badly hacked together constraint information for Ted's cell production optimizer
 */
void tedConstraints(int &t11, float &t12, int &t13, Array1D<int> *t2, int &v11, float &v12, int &v13, Array1D<int> *v2, Array1D<int> *geometryDiff, bool &constraint1, bool &constraint2, ProgramOptionParser *opt) {
	GeometryProfile gP;
	opt->getGeometryProfile(&gP);
	t11 = 400; // maximum number of cells in MR
	v11 = gP.maxSizeOfMR();		 // actual number of cells in MR
	t12 = -0.1f;			 // threshold for monotonicity
	v12 = -gP.isMonotonicIncreasing();  // -1 if geometry is monotonically increasing, 0 otherwise
	gP.getCurrentGeometry(0,geometryDiff);
	t13 = 100;				// maximum difference in number of cells between row 1, row 8
	v13 = abs((*geometryDiff)(0)-(*geometryDiff)(geometryDiff->size()-1));
	for (int i=0;i<geometryDiff->size()-1;i++) {
		(*geometryDiff)(i) = abs((*geometryDiff)(i)-(*geometryDiff)(i+1));
	}
	geometryDiff->resize(geometryDiff->size()-1);
	t2->resize(geometryDiff->size());
	t2->fill(-10);
	v2->resize(geometryDiff->size());
	for (int i=0;i<v2->size();i++) {
		(*v2)(i)=-(*geometryDiff)(i);
	}

	constraint1 = (v11<t11) & (v12<t12) & (v13<t13);
	constraint2=true;
	for (int i=0;i<v2->size();i++) {
		if ((*v2)(i)<(*t2)(i)) {
			constraint2=false;
		}
	}
}



/*!
 * Print simulation results to standard output in a format that is readable by Ted's MCMC setup
 */
void printSimulationResultsForTed_cellProduction(PedigreeDepthStorage *pdStorage, ProgramOptionParser *opt) {
	// print fitness metric values
	printf("\n\n\n\n");
	printf("time for N cells to exit the mitotic region = %f\n",pdStorage->time_nthMitoticExit->mean());
	printf("\n\n\n\n");
	printf("\n\n\n\n");
	printf("\n\n\n");

	// print simulated fecundity schedule
	printf("\n");
	printf("\n\n\n\n\n\n\n");
	printf("\n");

	// print simulated cell exiting meiotic region schedule
	printf("\n");

	printf("\n\n\n\n\n\n\n");
	printf("\n");

	// print simulated cell production schedule
	printf("\n");
	printf("\n\n\n");
	printf("\n");

	// get constraints
	int t11, t13,v11,v13;
	float t12,v12;
	Array1D<int> t2, v2, geometryDiff;
	bool constraint1, constraint2;
	tedConstraints(t11, t12, t13, &t2, v11, v12, v13, &v2, &geometryDiff, constraint1, constraint2, opt);


	// print aggregate max num cells constraint
	printf("constraint set #1: %d\n", constraint1);

	// print aggregate monotonic geometry constraint
	printf("constraint set #2: %d\n", constraint2);

	// print aggregated constraint
	printf("all constraints met: %d\n", constraint1 && constraint2);

	// pedigree depth for individual runs
	printf("\n****** Individual Output ******\n");
	printf("time for N cells to exit MR = ");
	pdStorage->time_nthMitoticExit->display();

	// print whether cell production constraint is met for individual runs
	printf("aggregated constraint set #1 = ");
	for (int x=0;x<opt->numSim;x++) {
		printf("%d",constraint1);
		if (x!=opt->numSim-1) {
			printf(",");
		}
	}
	printf("\n");

	// print whether monotonic constraint is met for individual runs
	printf("aggregated constraint set #2 = ");
	for (int x=0;x<opt->numSim;x++) {
		bool m=true;
		for (int i=0;i<v2.size();i++) {
			if (v2(i)<t2(i)) {
				m=false;
			}
		}
		printf("%d",m);
		if (x!=opt->numSim-1) {
			printf(",");
		}
	}
	printf("\n");

	// print number of cells in simulation
	printf("num cells in MR by simulation (constraint set #1, accept if < constraint=%d) = ", t11);
	for (int j=0;j<opt->numSim;j++) {
		printf("%d",v11);
		if (j!=opt->numSim-1) {
			printf(",");
		}
	}
	printf("\n");

	// print whether monotonic
	printf("geometry is monotonic simulation (constraint set #1, accept if < constraint=%f) = ", t12);
	for (int x=0;x<opt->numSim;x++) {
		printf("%f",v12);
		if (x!=opt->numSim-1) {
			printf(",");
		}
	}
	printf("\n");

	// print difference in number of cells in first/last cell rows

	printf("row 1-8 difference (constraint set #1, accept if < constraint=%d) = ", t13);
	for (int x=0;x<opt->numSim;x++) {
		printf("%d",v13);
		if (x!=opt->numSim-1) {
			printf(",");
		}
	}
	printf("\n");

	// print difference in number of cells in adjacent cell rows
	for (int i=0;i<geometryDiff.size();i++) {
		printf("row %d-%d difference (constraint set #2, accept if > constraint=%d) = ",i,i+1,t2(1));
		for (int x=0;x<opt->numSim;x++) {
			printf("%d",v2(i));
			if (x!=opt->numSim-1) {
				printf(",");
			}
		}
		printf("\n");
	}

}


/*!
 * Print simulation results to standard output when MR constraints are not satisfied.
 * These results fixed to very high pedigree depths / non-satisfied constraints
 */
void printSimulationResultsForTed_pedigreeDepth_badMR(PedigreeDepthStorage *pdStorage, ProgramOptionParser *opt) {

	// print pedigree depth (i.e., value that Ted's optimizer is minimizing)
	printf("999\n");

	// print whether cell production constraints (cpc) met
	printf("0\n");

	// print whether fecundity constraints (fc) met
	printf("0\n");

	// print whether all constraints met
	printf("0\n");

	// print pedigree depth of individual simulations runs
	Array1D<float> temp;
	Array1D<float>::array_index numRuns = pdStorage->pd_AllCellsEnterMeioticRegion->size();
	temp.resize(numRuns);
	temp.fill(999);
	temp.display();

	// print whether cell production constraint is met for individual simulation runs
	Array2D<float>::array_index dimx_cellProduction,dimy_cellProduction;
	pdStorage->simulatedCellProductionSchedule->getDimensions(dimx_cellProduction,dimy_cellProduction);
	for (int x=0;x<dimx_cellProduction;x++) {
		printf("0");
		if (x!=dimx_cellProduction-1) {
			printf(",");
		}
	}
	printf("\n");

	// print whether fecundity constraint is met for individual simulation runs
	Array2D<float>::array_index dimx_simulatedFecundity, dimy_simulatedFecundity;
	pdStorage->simulatedFecunditySchedule->getDimensions(dimx_simulatedFecundity,dimy_simulatedFecundity);
	for (Array2D<float>::array_index x=0;x<dimx_simulatedFecundity;x++) {
		printf("%d",false);
		if (x!=dimx_simulatedFecundity-1) {
			printf(",");
		}
	}
	printf("\n");

	// print cell production time constraints for individual runs
	for (int i=0;i<dimy_cellProduction;i++) {
		printf("%f",(*pdStorage->cellProductionConstraints)(0,i));
		if (i!=dimy_cellProduction-1) {
			printf(",");
		}
	}
	printf("\n");

	// print cell production times for individual runs
	for (int i=0;i<dimy_cellProduction;i++) {
		for (int j=0;j<dimx_cellProduction;j++) {
			printf("999");
			if (j!=dimx_cellProduction-1) {
				printf(",");
			}
		}
		printf("\n");
	}

	// print fecundity constraints for individual runs
	int sumConstraint=0;
	for (int i=0;i<dimy_simulatedFecundity;i++) {
		int constraint;
		if (i==0) {
			constraint = opt->numSperm;
		} else {
			constraint = boost::math::iround((*pdStorage->fecundityConstraints)(i));
		}
		sumConstraint+=constraint;
		printf("%d",sumConstraint);
		if (i!=dimy_simulatedFecundity-1) {
			printf(",");
		}
	}
	printf("\n");

	// print fecundity for individual runs
	Array1D<float> sumGamete;
	sumGamete.resize(dimx_cellProduction);
	for (int i=0;i<dimy_simulatedFecundity;i++) {
		for (int j=0;j<dimx_cellProduction;j++) {
			sumGamete(j)+=(*pdStorage->simulatedFecunditySchedule)(j,i);
			printf("%d",boost::math::iround(sumGamete(j)));
			if (j!=dimx_cellProduction-1) {
				printf(",");
			}
		}
		printf("\n");
	}
}


/*!
 * Print simulation results to standard output
 */
void printSimulationResultsForTed_pedigreeDepth(PedigreeDepthStorage *pdStorage, bool cellProductionConstraintsMet, bool fecundityScheduleConstraintsMet, bool allConstraintsMet, ProgramOptionParser *opt) {

	// print pedigree depth (i.e., value that Ted's optimizer is minimizing)
	float valueBeingOptimized = 0;
	string option_str(opt->fitnessMetric.c_str());
	if      (option_str.find("meanPedigreeDepthMeioticCells:")!=string::npos) {valueBeingOptimized = pdStorage->pd_AllCellsEnterMeioticRegion->mean();}
	else if (strcmp(option_str.c_str(),"spermOocyteMean")==0) {valueBeingOptimized = pdStorage->pd_spermOocyte->mean();}
	else if (strcmp(option_str.c_str(),"spermOocyteSum")==0) {valueBeingOptimized = pdStorage->pd_spermOocyteSum->mean();}
	else if (strcmp(option_str.c_str(),"generationRate")==0) {float r = calculateGenerationFitness(pdStorage->simulatedFecunditySchedule,pdStorage->fecundityConstraints,pdStorage->pd_spermOocyte->mean(),opt->deleteriousMutationRate,opt->maximalGenerationRate,opt->maximalPedigreeDepth);valueBeingOptimized = -r;}
        else if (strcmp(option_str.c_str(),"generationRate2")==0) {float r = calculateGenerationFitness_selectionCoefficient(pdStorage->simulatedFecunditySchedule,pdStorage->fecundityConstraints,pdStorage->pd_spermOocyte->mean(),opt->deleteriousMutationRate,opt->maximalGenerationRate,opt->maximalPedigreeDepth);valueBeingOptimized = -r;}

	else {std::cerr << "Invalid choice of fitnessScoringMethod" << endl;exit(1);}
	printf("%f\n",valueBeingOptimized);

	// print whether cell production constraints (cpc) met
	printf("%d\n", cellProductionConstraintsMet);

	// print whether fecundity constraints (fc) met
	printf("%d\n", fecundityScheduleConstraintsMet);

	// print whether all constraints met
	printf("%d\n", allConstraintsMet);

	// print pedigree depth of individual simulations runs
	if      (option_str.find("meanPedigreeDepthMeioticCells:")!=string::npos) {pdStorage->pd_AllCellsEnterMeioticRegion->display();}
	else if (strcmp(option_str.c_str(),"spermOocyteMean")==0) {pdStorage->pd_spermOocyte->display();}
	else if (strcmp(option_str.c_str(),"spermOocyteSum")==0) {pdStorage->pd_spermOocyteSum->display();}
	else if (strcmp(option_str.c_str(),"generationRate")==0) {
		Array1D<float> output;
		output.resize(pdStorage->pd_spermOocyte->size());
        float r = calculateGenerationFitness(pdStorage->simulatedFecunditySchedule,pdStorage->fecundityConstraints,pdStorage->pd_spermOocyte->mean(),opt->deleteriousMutationRate,opt->maximalGenerationRate,opt->maximalPedigreeDepth);
		for (int i=0;i<pdStorage->pd_spermOocyte->size();i++) {
			output(i) = -r;
		}
		output.display();
	}
        else if (strcmp(option_str.c_str(),"generationRate2")==0) {
                Array1D<float> output;
                output.resize(pdStorage->pd_spermOocyte->size());
        float r = calculateGenerationFitness_selectionCoefficient(pdStorage->simulatedFecunditySchedule,pdStorage->fecundityConstraints,pdStorage->pd_spermOocyte->mean(),opt->deleteriousMutationRate,opt->maximalGenerationRate,opt->maximalPedigreeDepth);
                for (int i=0;i<pdStorage->pd_spermOocyte->size();i++) {
                        output(i) = -r;
                }
                output.display();
        }

	else {std::cerr << "Invalid choice of fitnessScoringMethod" << endl;exit(1);}

	// print whether cell production constraint is met for individual simulation runs
	Array2D<float>::array_index dimx_cellProduction,dimy_cellProduction;
	pdStorage->simulatedCellProductionSchedule->getDimensions(dimx_cellProduction,dimy_cellProduction);
	for (int x=0;x<dimx_cellProduction;x++) {
		bool constraintsSatisfied=true;
		for (int y=0;y<dimy_cellProduction;y++) {
			if ((*pdStorage->simulatedCellProductionSchedule)(x,y)>(*pdStorage->cellProductionConstraints)(0,y)) {
				constraintsSatisfied=false;
				break;
			}
		}
		if (opt->constraintMethod.compare("none")==0){
			constraintsSatisfied=true;
		}
		printf("%d",constraintsSatisfied);
		if (x!=dimx_cellProduction-1) {
			printf(",");
		}
	}
	printf("\n");

	// print whether fecundity constraint is met for individual simulation runs
	Array2D<float>::array_index dimx_simulatedFecundity, dimy_simulatedFecundity;
	pdStorage->simulatedFecunditySchedule->getDimensions(dimx_simulatedFecundity,dimy_simulatedFecundity);
	for (Array2D<float>::array_index x=0;x<dimx_simulatedFecundity;x++) {
		bool constraintsSatisfied=true;
		int sum_sim = 0;
		int sum_sched = 0;
		for (int y=0;y<dimy_simulatedFecundity;y++) {
			if (y==0) {  			// check sperm constraint
				int numSperm_sim = boost::math::iround((*pdStorage->simulatedFecunditySchedule)(x,y));
				int numSperm_sched = opt->numSperm; // boost::math::iround((*pdStorage->fecundityConstraints)(0,y));
				if (numSperm_sim!=numSperm_sched) {
					constraintsSatisfied=false;
					break;
				}
			} else {				// check oocyte constraint
				sum_sim+=ceil((*pdStorage->simulatedFecunditySchedule)(x,y));
				sum_sched+=boost::math::iround((*pdStorage->fecundityConstraints)(y));
				if (sum_sim<sum_sched) {
					constraintsSatisfied=false;
					break;
				}
			}
		}
		if (opt->constraintMethod.compare("none")==0){
			constraintsSatisfied=true;
		}
		printf("%d",constraintsSatisfied);
		if (x!=dimx_simulatedFecundity-1) {
			printf(",");
		}
	}
	printf("\n");

	// print cell production time constraints for individual runs
	for (int i=0;i<dimy_cellProduction;i++) {
		printf("%f",(*pdStorage->cellProductionConstraints)(0,i));
		if (i!=dimy_cellProduction-1) {
			printf(",");
		}
	}
	printf("\n");

	// print cell production times for individual runs
	for (int i=0;i<dimy_cellProduction;i++) {
		for (int j=0;j<dimx_cellProduction;j++) {
			printf("%f",(*pdStorage->simulatedCellProductionSchedule)(j,i));
			if (j!=dimx_cellProduction-1) {
				printf(",");
			}
		}
		printf("\n");
	}

	// print fecundity constraints for individual runs
	int sumConstraint=0;
	for (int i=0;i<dimy_simulatedFecundity;i++) {
		int constraint;
		if (i==0) {
			constraint = opt->numSperm;
		} else {
			constraint = boost::math::iround((*pdStorage->fecundityConstraints)(i));
		}
		sumConstraint+=constraint;
		printf("%d",sumConstraint);
		if (i!=dimy_simulatedFecundity-1) {
			printf(",");
		}
	}
	printf("\n");

	// print fecundity for individual runs
	Array1D<float> sumGamete;
	sumGamete.resize(dimx_cellProduction);
	for (int i=0;i<dimy_simulatedFecundity;i++) {
		for (int j=0;j<dimx_cellProduction;j++) {
			sumGamete(j)+=(*pdStorage->simulatedFecunditySchedule)(j,i);
			printf("%d",boost::math::iround(sumGamete(j)));
			if (j!=dimx_cellProduction-1) {
				printf(",");
			}
		}
		printf("\n");
	}
}



/*!
 * Check whether cell division constraints are met
 */
bool checkCellProductionConstraints(PedigreeDepthStorage *pdStorage, ProgramOptionParser *opt) {
	// check whether cell production constraints are met
	bool cellProductionConstraintsMet=true;
	if (opt->constraintMethod.compare("mean")==0){
		Array2D<float>::array_index dimx,dimy;
		pdStorage->simulatedCellProductionSchedule->getDimensions(dimx,dimy);
		for (int i=0;i<dimy;i++) {
			float time = pdStorage->simulatedCellProductionSchedule->mean(i);
			float constraint = (*pdStorage->cellProductionConstraints)(0,i);
			if (time>constraint) {
				cellProductionConstraintsMet=false;
			}
		}
	} else if (opt->constraintMethod.compare("median50")==0 || opt->constraintMethod.compare("median90")==0) {
		// count number of simulations that satisfy constraints
		Array2D<float>::array_index dimx,dimy;
		pdStorage->simulatedCellProductionSchedule->getDimensions(dimx,dimy);
		float numSatisfied=0;
		for (int x=0;x<dimx;x++) {
			bool constraintsSatisfied=true;
			for (int y=0;y<dimy;y++) {
				if ((*pdStorage->simulatedCellProductionSchedule)(x,y)>(*pdStorage->cellProductionConstraints)(0,y)) {
					constraintsSatisfied=false;
					break;
				}
			}
			if (constraintsSatisfied) {
				numSatisfied++;
			}
		}
		// check whether enough simulations satisfy constraints
		float threshold = -1;
		if (opt->constraintMethod.compare("median50")==0) {threshold = 0.5f;}
		if (opt->constraintMethod.compare("median90")==0) {threshold = 0.9f;}
		if (numSatisfied/(float)dimx < threshold) {
			cellProductionConstraintsMet = false;
		}
	} else if (opt->constraintMethod.compare("none")==0){
		cellProductionConstraintsMet = true;
	} else {
		printf("error: invalid choice of constraintMethod (%s)\n",opt->constraintMethod.c_str());
		exit(1);
	}
	return cellProductionConstraintsMet;
}

/*!
 * Check whether fecundity constraints are met
 */
bool checkFecundityConstraints(PedigreeDepthStorage *pdStorage, ProgramOptionParser *opt) {
	bool fecundityScheduleConstraintsMet = true;
	float totalGametesProduced=0;
	float totalGametesProduced_simulated=0;
	if (opt->constraintMethod.compare("mean")==0){
		Array2D<float>::array_index dimx,dimy;
		pdStorage->simulatedFecunditySchedule->getDimensions(dimx,dimy);
		for (int y=0;y<dimy;y++) {
			float meanGametesProduced_simulated = pdStorage->simulatedFecunditySchedule->mean(y);
			float meanGametesProduced = (*pdStorage->fecundityConstraints)(y);
			if (y==0) {
				meanGametesProduced = (float)opt->numSperm;
			}
			totalGametesProduced+=meanGametesProduced;
			totalGametesProduced_simulated+=meanGametesProduced_simulated;

			// ceiling to prevent underflow
			int round_totalGametesProduced = boost::math::iround(totalGametesProduced);
			int round_totalGametesProduced_simulated = boost::math::iround(ceil(totalGametesProduced_simulated));

			if (round_totalGametesProduced_simulated<round_totalGametesProduced) {
				fecundityScheduleConstraintsMet=false;
				break;
			}
		}
	} else if (opt->constraintMethod.compare("median50")==0 || opt->constraintMethod.compare("median90")==0) {
		// count number of simulations that satisfies constraints
		Array2D<float>::array_index dimx,dimy;
		pdStorage->simulatedFecunditySchedule->getDimensions(dimx,dimy);
		float numSatisfied = 0;
		for (int x=0;x<dimx;x++) {
			float cumulativeProgeny = 0;
			float cumulativeConstraint = 0;
			bool constraintSatisfied = true;
			for (int y=0;y<dimy;y++) {
				float constraint;
				if (y==0) {
					constraint = (*pdStorage->fecundityConstraints)(0,0);
				} else {
					constraint =  (*pdStorage->fecundityConstraints)(y);
				}
				cumulativeConstraint += constraint;
				cumulativeProgeny += (*pdStorage->simulatedFecunditySchedule)(x,y);
				if (cumulativeProgeny<cumulativeConstraint) {
					constraintSatisfied = false;
				}
			}
			if (constraintSatisfied) {
				numSatisfied++;
			}
		}

		float threshold = -1;
		if (opt->constraintMethod.compare("median50")==0) {threshold = 0.5f;}
		if (opt->constraintMethod.compare("median90")==0) {threshold = 0.9f;}
		if (numSatisfied/(float)dimx < threshold) {
			fecundityScheduleConstraintsMet = false;
		}
	} else if (opt->constraintMethod.compare("none")==0) {
		fecundityScheduleConstraintsMet = true;
	} else {
		printf("error: invalid choice of constraintMethod (%s)\n",opt->constraintMethod.c_str());
		exit(1);
	}
	return fecundityScheduleConstraintsMet;
}

/*!
 * Aggregate fitness metric across simulations
 */
float calculateFitness(const char *option, PedigreeDepthStorage *pdStorage, ProgramOptionParser *opt) {
	float rn = 0;
	string option_str(option);
	if (strcmp(option,"oocyteMedian")==0) {
		rn = pdStorage->pd_oocyte->perctile(0.5);
	} else if (strcmp(option,"oocyteMean")==0) {
		rn = pdStorage->pd_oocyte->mean();
	} else if (strcmp(option,"spermOocyteMedian")==0) {
		rn = pdStorage->pd_spermOocyte->perctile(0.5);
	} else if (strcmp(option,"spermOocyteMean")==0) {
		rn = pdStorage->pd_spermOocyte->mean();
	} else if (strcmp(option,"spermOocyteSum")==0) {
		rn = pdStorage->pd_spermOocyteSum->mean();
	} else if (strcmp(option,"mitoticRegionMean")==0) {
		rn = pdStorage->pd_mitoticRegion->mean();
	} else if (strcmp(option,"wholeGermlineMean")==0) {
		rn = pdStorage->pd_wholeGermline->mean();
	} else if (strcmp(option,"germlineAtEggLaying")==0) {
		rn = pdStorage->pd_GermlineAtEggLaying->mean();
	} else if (strcmp(option,"generationRate")==0) {
		float r = calculateGenerationFitness(pdStorage->simulatedFecunditySchedule,pdStorage->fecundityConstraints,pdStorage->pd_spermOocyte->mean(),opt->deleteriousMutationRate,opt->maximalGenerationRate,opt->maximalPedigreeDepth);
		rn = -r; // negative since optimizer finds minimum
	} else if (strcmp(option,"generationRate2")==0) {
                float r = calculateGenerationFitness_selectionCoefficient(pdStorage->simulatedFecunditySchedule,pdStorage->fecundityConstraints,pdStorage->pd_spermOocyte->mean(),opt->deleteriousMutationRate,opt->maximalGenerationRate,opt->maximalPedigreeDepth);
                rn = -r; // negative since optimizer finds minimum 
	} else if (option_str.find("meanPedigreeDepthMeioticCells:")!=string::npos) {
		rn = pdStorage->pd_AllCellsEnterMeioticRegion->mean();
	} else if (option_str.find("tedNIPS:")!=string::npos) {
		rn = pdStorage->time_nthMitoticExit->mean();
	} else {
		std::cerr << "Invalid choice of fitnessScoringMethod" << endl;
		exit(1);
	}
	return rn;
}


/*!
 * Used for histogram binning
 */
int getClosestIdx(Array1D<float> *vec, float val) {
	Array1D<float>::array_index sz=vec->size();
	float closestDistance=INF;
	int closestIdx=-1;
	for (Array1D<float>::array_index i=0;i<sz;i++) {
		if (abs((*vec)(i)-val)<closestDistance) {
			closestDistance=abs((*vec)(i)-val);
			closestIdx=(int)i;
		}
	}
	return closestIdx;
}
