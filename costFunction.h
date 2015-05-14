/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/

#include <stdio.h>
#include "definitions.h"
#include "ProgramOptionParser.h"
#include "util.h"
#include "GermLine.h"
#include "util_germline.h"
#include <climits>
#if defined USE_LIBDISPATCH
#include <dispatch/dispatch.h>
#endif

extern "C"
float costFunction_pedigreeDepth(int argc, char *argv[], float resultStats[]);

void * run_thread(void *arg) __attribute__ ((noreturn));
