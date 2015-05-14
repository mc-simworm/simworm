/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#include "costFunction.h"
#include "gitVersion.h"


extern "C"
void getVersion(char c[]);

extern "C"
float costFunction(int argc, char *argv[], float resultStats[]);
