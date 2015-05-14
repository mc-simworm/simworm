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
<<<<<<< HEAD
#if defined USELIBDISPATCH
#include <dispatch/dispatch.h>
=======
#if defined USELIBDISPATCH 
#include "dispatch/dispatch.h"
>>>>>>> 43d79ab0b3d557b2228c8868f6a72fa692e0283a
#endif


void ccfit(int argc, char *argv[], FitData *fD, int chaseTimeIdx);
