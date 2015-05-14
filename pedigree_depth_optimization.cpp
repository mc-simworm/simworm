/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#include "definitions.h"
#include "costFunction.h"

/*!
 * Entry point for simworm binary
 */
int main(int argc, char *argv[]) {
	float float_array[SZRESULTSTATS]={0};
	costFunction_pedigreeDepth(argc,argv,float_array);
	return 0;
}

