/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/


#include "run.h"

void getVersion(char c[]){
	pid_t p=getpid();
	sprintf(c,"pid=%i, "
        		"commit = %s; "
        		"Added functionality for cells to enter pre-meiotic state in mitotic region, see program option --meioticFractionProfile"
        		,p, gitVersion);
}

/*!
 * Entry point for simworm shared library
 */
float costFunction(int argc, char *argv[], float resultStats[]) {
    return costFunction_pedigreeDepth(argc, argv, resultStats);
}
