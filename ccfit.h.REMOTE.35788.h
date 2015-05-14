/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#include "ProgramOptionParser.h"
#include "GermLine.h"
#include "FitData.h"
#if defined USE_LIBDISPATCH 
#include "dispatch/dispatch.h"
#endif




void ccfit(int argc, char *argv[], FitData *fD, int chaseTimeIdx);
