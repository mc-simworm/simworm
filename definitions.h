/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#define TIME_EGGLAY_L2 26				// time from egg layihng to L2, used to calculate transition matrix for generation rate
#define SZRESULTSTATS 10				// size of vector that is returned to refining grid search optimizer
#define INF 99999						// value of a very big number
#define CCRESCALEINTERVAL 5 			// rescale cell cycle lengths every n divisions
#define EPS 0.000001f					// value of a very small number
#define PEDIGREE_DEPTH_0 10.124672f      // pedigree depth of experimentally measured  profile {4.23},{2.83},{6.71},{4.43},{6.71},{4.43}
#define M_SEARCH_RANGE 5
//#define MAX_MITOTIC_REGION_SIZE 4000	// maximum number of cells in the mitotic region, return -1 if exceeds
//#define FITNESS_0 1.089677f		  	    // ceiling on hourly rate of increase, based on experimentally measured profile // 1.093593
//#define FITNESS_COST_PER_MUTATION 0.0066f // percent decrease in r per mutation per generation.  Drawn from Vassilieva 1999 (The Rate of Spontaneous Mutation for Life-History Traits).  (1-(1.35-0.0023)/1.35)/0.26
//#define NUMBER_OF_MUTATIONS_0 0.26f 	 	 // Number of mutations per genome per generation
//#define DELETERIOUS_MUTATION_RATE 0.03f 	 	 // Number of deleterious mutations per genome per generation.  Drawn from Drawn from Vassilieva 1999 (The Rate of Spontaneous Mutation for Life-History Traits).  (1-(1.35-0.0023)/1.35)/0.26
