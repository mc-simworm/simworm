/*******************************************************************************
  * simworm v0.1
  * Copyright (c) 2009-2013 Michael Chiang.
  * All rights reserved. This program and the accompanying materials
  * are made available under a dual license: the two-clause BSD license or
  * the GNU Public License v2.
  ******************************************************************************/
#include "ProgramOptionParser.h"
using namespace std;

/*
 * Constructor
 */
ProgramOptionParser::ProgramOptionParser(string mode) {
	pedigreeDepthOrCellCycleFit=mode;
}

/*
 * Destructor
 */
ProgramOptionParser::~ProgramOptionParser() {
}

/*
 * Parse program options
 */
void ProgramOptionParser::parseArgs (int argc, char *argv[]) {
		// program option list
		PO::options_description desc("Options");
		desc.add_options()
		    // these are optional parameters implemented for user convenience
			("help,h", "Display help")
			("debug", "Include this option if you want to write debugging output to standard output.  Example usage: --debug")
			("verbose", "Include this option if you want to write simulation results to the standard output.  Example usage: --verbose")
			("tedOptimizerCellProduction", "Include this option if you want to write simulation results to the standard output in a form parsable by Ted's optimizer.  Example usage: --tedOptimizerCellProduction")
			("tedOptimizerPedigreeDepth", "Include this option if you want to write simulation results to the standard output in a form parsable by Ted's optimizer.  Example usage: --tedOptimizerPedigreeDepth")

			// these are mandatory parameters that must be specified whether you are running pedigree depth simulations or cell cycle fits
			("numSim", PO::value<int>(), "Number of simulations to run.  Example usage: --numSim 200")
			("noise", PO::value<float>(), "Noise in cell cycle lengths.  Example usage: --noise 0.3")
			("randomSeed", PO::value<int>(), "Unsigned int (>=0) used to seed random number generator.  Use -1 to use random seed.  Example usage: --randomSeed -1")
			("cellCycleProfile", PO::value<string>(), "Formatted cell cycle parameters.  Example usage: --cellCycleLength {3},{3},{4.3},{3},{6.6},{4.5},{6.6},{4.5}*0:0,14/327:1,15/1200:1,11/3000:1,11*linear*23")
			("mitoticRegionGeometry", PO::value<string>(), "Formatted mitotic region geometry.  Example usage: --geometry 327:5,5,6,7,7,8,8,9,9,9,10,10,10,10,10,10,11,11,11,11,11,10,10/1200:6,7,8,10,10,11,13,13,14,15,16,16,16,17,17,17,16,16,15,0,0,0,0/3000:9,10,13,15,16,17,19,19,20,21,21,21,20,0,0,0,0,0,0,0,0,0,0")
			("truncatedMitoticIndex", PO::value<string>(), "Formatted truncated mitotic index used to determine threshold for pre-meiotic G2 arrest.  A value of x means if the local mitotic index exceeds x, then cells will arrest in pre-meiotic G2.  Example usage: --meioticFractionProfile 327:1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.039700,0.036000,0.034500,0.036500,0.016500,0.009600,0.013300,0")

			// These are parameters specific to pedigree depth simulations.  You do not have to specify these parameters if running cell cycle fits
			("cellProductionConstraints", PO::value<string>(), "Cell production constraints.  Example usage: --constraints 327:30/1200:47.75/3000:95.75")
			("oocyteProbability", PO::value<string>(), "Probability a germ cell becomes an oocyte instead of undergoing apoptosis when it exits the meiotic region.  Example usage: --oocyteProbability step*0:0.5739/30:0.0948/638.458:0.0686")
			("numSperm", PO::value<int>(), "Number of sperm to produce.  Example usage: --numSperm 33")
			("numOocyte", PO::value<int>(), "Number of oocytes to produce.  Example usage: --numOocyte 131")
			("meioticRegionGeometry", PO::value<string>(), "Serialized meiotic region geometry.  Example usage: --meioticRegionGeometry 327:119/1200:900/3000:1250")
			("fitnessMetric", PO::value<string>(), "Fitness metric returned by main().  Commonly used arguments: spermOocyteMean, meanPedigreeDepthMeioticCells:N, tedNIPS:N.  Example usage: --fitnessMetric spermOocyteMean")
			("constraintMethod", PO::value<string>(), "Method used to determine whether simulations fulfill constraint. Choices are mean, percentile50, percentile90, none.  Example usage: --constraintMethod mean'")
			("forwardMovementProbability", PO::value<float>(), "Probability of a cell moving forward instead of side to side.  Example usage: --forwardMovementProbability 0.5")
			("fecunditySchedule", PO::value<string>(), "Serialized oocyte production schedule.  Example usage: --fecunditySchedule 0:30/17:47.75/58:71.75/42:95.75/12:119.75/1:143.75/1:167.75")
			("deleteriousMutationRate", PO::value<float>(), "Deleterious mutation rate per genome per generation.  Used to calculate fitness (combining generation rate + pedigree depth).  Example usage: --deleteriousMutationRate 0.03")
			("maximalGenerationRate", PO::value<float>(), "Generation rate (units = hours^-1) of flat profile, used to calculate aggregate fitness.  Example usage: --maximalGenerationRate 1.089677")
			("pedigreeDepthIncrementMethod", PO::value<int>(), "Method of incrementing pedigree depth.  0 for +1, 1 for +1/ccLengthExample, 2 for +1/normalizedCcLength.  Example usage: --pedigreeDepthIncrementMethod 0")
			("maxMitoticRegionSize", PO::value<int>(), "Maximum mitotic region size, program returns -1 if at any point MR simulation exceeds N cells.  Example usage: --maxMitoticRegionSize 4000")
			("filledGermline", "This is a hacked in parameter to do pedigree depth simulations starting from a filled germline.  This was used to calculate how many cells are produced in eight hours by the mitotic region to compare against Fox 2011.  If not specified, pedigree depth simulations default to start with one cell.  Has no effect on cell cycle fits.  Example usage: --filledGermline")
			("maximalPedigreeDepth", PO::value<float>(), "Pedigree depth of flat profile, used to calculate aggregate fitness.  Example usage: --maximalPedigreeDepth 0")

			// These are parameters specific to cell cycle fit simulations
			("chaseTime", PO::value<string>(), "Chase time after EdU pulse")
			("mitoticPhaseIndices", PO::value<string>(), "Mitotic phase indices")
			("dnaBin", PO::value<string>(), "DNA bins")
			("fitData", PO::value<string>(), "fitData, order is chase time, dna bin, emm, flm")
			("dnaNoise", PO::value<float>(), "Noise in DNA quantification")
			;
		PO::variables_map vm;
		PO::store(PO::parse_command_line(argc, argv, desc), vm);
		PO::notify(vm);

		/*
		 * Read optional program options
		 */
		if (vm.count("help")) {
			cout << desc << endl;
			cout << "*************" << endl;
			cout << "Example usage:" << endl;
			cout << "./simworm_darwin --numSim 10 --noise 0.01 --mitoticRegionGeometry 327:5,5,6,7,7,8,8,9,9,9,10,10,10,10,10,10,11,11,11,11,11,10,10/1200:6,7,8,10,10,11,13,13,14,15,16,16,16,17,17,17,16,16,15,0,0,0,0/3000:9,10,13,15,16,17,19,19,20,21,21,21,20,0,0,0,0,0,0,0,0,0,0 --randomSeed -1 --cellCycleProfile {3},{3},{4.3},{3},{6.6},{4.5},{6.6},{4.5}*0:1,15/327:1,15/1200:1,11/3000:1,11*linear --cellProductionConstraints 327:30/1200:47.75/3000:95.75 --oocyteProbability step*0:0.5739/30:0.0948/638.458:0.0686 --numSperm 33 --numOocyte 131 --fecunditySchedule 0:30/17:47.75/58:71.75/42:95.75/12:119.75/1:143.75/1:167.75 --meioticRegionGeometry 327:119/1200:900/3000:1250 --fitnessMetric spermOocyteMean --constraintMethod mean --forwardMovementProbability 0.5 --deleteriousMutationRate 0.03 --maximalGenerationRate 1.089677 --pedigreeDepthIncrementMethod 0 --maxMitoticRegionSize 4000 --meioticFractionProfile 327:1,15,23/1200:1,11,19 --verbose" << endl;
			exit(1);
		}
		if (vm.count("verbose")) {
			verbose = true;
		} else {
			verbose = false;
		}
		if (vm.count("filledGermline")) {
			filledGermline = true;
		} else {
			filledGermline = false;
		}
		if (vm.count("tedOptimizerCellProduction")) {
			verbose_tedCellProduction = true;
		} else {
			verbose_tedCellProduction = false;
		}
		if (vm.count("tedOptimizerPedigreeDepth")) {
			verbose_tedPedigreeDepth = true;
		} else {
			verbose_tedPedigreeDepth = false;
		}
		if (vm.count("debug")) {
			debug = true;
		} else {
			debug = false;
		}
		if (vm.count("maximalPedigreeDepth")) {
			maximalPedigreeDepth=vm["maximalPedigreeDepth"].as<float>();
		} else {
			maximalPedigreeDepth = PEDIGREE_DEPTH_0;
		}

		/*
		 * Read required program options
		 */
		cellCycleProfile = vm["cellCycleProfile"].as<string>();
		mitoticRegionGeometry = vm["mitoticRegionGeometry"].as<string>();
		numSim = vm["numSim"].as<int>();
		if (debug) {numSim=1;}
		noise = vm["noise"].as<float>();
		randomSeed = vm["randomSeed"].as<int>();
		truncatedMitoticIndex = vm["truncatedMitoticIndex"].as<string>();
		mitoticPhaseIndices=vm["mitoticPhaseIndices"].as<string>();

		/*
		 * Read program options specific to pedigree depth simulations
		 */
		if (pedigreeDepthOrCellCycleFit.compare("pd")==0) {
			cellProductionConstraints = vm["cellProductionConstraints"].as<string>();
			oocyteProbability = vm["oocyteProbability"].as<string>();
			fecunditySchedule = vm["fecunditySchedule"].as<string>();
			numSperm = vm["numSperm"].as<int>();
			numOocyte = vm["numOocyte"].as<int>();
			fitnessMetric = vm["fitnessMetric"].as<string>();
			constraintMethod = vm["constraintMethod"].as<string>();
			forwardMovementProbability = vm["forwardMovementProbability"].as<float>();
			meioticRegionGeometry=vm["meioticRegionGeometry"].as<string>();
			deleteriousMutationRate = vm["deleteriousMutationRate"].as<float>();
			maximalGenerationRate = vm["maximalGenerationRate"].as<float>();
			pedigreeDepthIncrementMethod = vm["pedigreeDepthIncrementMethod"].as<int>();
			maxMitoticRegionSize = vm["maxMitoticRegionSize"].as<int>();

			// Minimum number of divisions to run should be last cell production constraint
			Array2D<float> c;
			getCellProductionConstraints(&c);
			Array2D<float>::array_index dimx,dimy;
			c.getDimensions(dimx,dimy);
			numDivisions=c(dimy-1)+1;

			// Minimum time to run should be last bin in fecundity schedule
			Array2D<float> f;
			getFecundityConstraints(&f);
			f.getDimensions(dimx,dimy);
			minTimeStr=boost::lexical_cast<string>(f(0,dimy-1));
			minTime=0;

			// Cell cycle fit program options.  These should not be specified in pedigree depth simulations
			//phaseIndices="none";
			dnaBin="none";
			fitData="none";
			dnaNoise=0;
		}

		/*
		 * Read program options specific to cell cycle fits
		 */
		else if (pedigreeDepthOrCellCycleFit.compare("ccfit")==0) {
			// Pedigree depth program options.  These should not be specified in cell cycle fit simulations
			cellProductionConstraints = "0:0";
			oocyteProbability = "linear*0:0";
			fecunditySchedule = "0:0,0";
			numSperm = 0;
			numOocyte = 0;
			numDivisions = 0;
			meioticRegionGeometry ="0:0";
			constraintMethod = "mean";
			forwardMovementProbability = 0.5;
			deleteriousMutationRate = 0;
			pedigreeDepthIncrementMethod = 0;
			maxMitoticRegionSize = INF;
			fitnessMetric = "none";
			maximalGenerationRate = 0;

			// ccfit parameters
			minTimeStr=vm["chaseTime"].as<string>();
			setMinTime(0);
			dnaBin=vm["dnaBin"].as<string>();
			if (vm.count("fitData")) {
				fitData=vm["fitData"].as<string>();
			} else {
				fitData="none";
			}
			dnaNoise=vm["dnaNoise"].as<float>();
		}
		else {
			cerr << "invalid mode choice" << endl;
		}
}

/*
 * Parse formatted DNA bins into a 1D array
 */
void ProgramOptionParser::getDnaBin(Array1D<float> *dnaBins) {
	dnaBins->parseString(dnaBin);
}

/*
 * Parse formatted chase times into a 1D array
 */
void ProgramOptionParser::getChaseTimes(Array1D<float> *chaseTimes) {
	chaseTimes->parseString(minTimeStr);
}

/*
 * Set the minimum running time of a cell cycle fit
 */
void ProgramOptionParser::setMinTime(int i) {
	Array1D<float> arg;
	arg.parseString(minTimeStr);
	minTime=arg(i);
}

/*
 * Display parsed program options
 */
void ProgramOptionParser::display() {
	cout << "Verbose: " << debug << endl;
	cout << "Number of simulations to run: " << numSim << endl;
	cout << "Cell cycle noise: " << noise << endl;

	// display cell cycle profile
	CellCycleProfile cc;
	getCellCycleProfile(&cc);
	cout << "Cell cycle profile:" << endl;
	cc.display();

	// display geometry
	GeometryProfile gg;
	getGeometryProfile(&gg);
	cout << "Geometry profile:" << endl;
	gg.display();
}

/*!
 * Parse serialized meiotic fraction profile into MeioticFractionProfile object
 */
void ProgramOptionParser::getMeioticFractionProfile(MeioticFractionProfile *mfp) {
	mfp->parseFormattedParameters(truncatedMitoticIndex,mitoticRegionGeometry);
}

/*
 * Parse ser9alized cell cycle profile into CellCycleProfile object
 */
void ProgramOptionParser::getCellCycleProfile(CellCycleProfile *cc) {
	int numCellRows = getNumCellRows();
	cc->parseFormattedParameters(cellCycleProfile, numCellRows);
}

/*
 * Parse seralized mitotic region geometry into GeometryProfile object
 */
void ProgramOptionParser::getGeometryProfile(GeometryProfile *gg) {
	gg->parseFormattedParameters(mitoticRegionGeometry);
}


/*
 * Get number of cell rows in the mitotic region including rows with 0 cells
 */
int ProgramOptionParser::getNumCellRows() {
	GeometryProfile gg;
	getGeometryProfile(&gg);
	return gg.numCellRows();
}

/*
 * Get number of chase times for cell cycle fits
 */
int ProgramOptionParser::getNumChaseTimes() {
	Array1D<float> chaseTimes;
	getChaseTimes(&chaseTimes);
	return (int)chaseTimes.size();
}

/*
 * Parse formatted phase indices into PhaseIndexProfile object
 */
void ProgramOptionParser::getPhaseIndicesProfile(PhaseIndexProfile *pp) {
	pp->parseFormattedParameters(mitoticPhaseIndices,mitoticRegionGeometry);
}

/*
 * Parse formatted cell production constraints into 2D matrix
 */
void ProgramOptionParser::getCellProductionConstraints(Array2D<float> *cellProductionConstraintsMatrix) {
	cellProductionConstraintsMatrix->parseString(cellProductionConstraints);
}

/*
 * Parse formatted meiotic region sizes into 2D matrix
 */
void ProgramOptionParser::getMeioticRegionSize(Array2D<int> *meioticRegionSz) {
	meioticRegionSz->parseString(meioticRegionGeometry);
}

/*
 * Parse formatted oocyte probabilities into 2D matrix
 */
void ProgramOptionParser::getOocyteProbability(Array2D<float> *probability) {
	vector<string> splitStr;
	string oocyteProbability_copy(oocyteProbability); // boost split appears to modify oocyteProbability when splitting which behaves badly when multithreading, so make a local copy to split
	boost::split(splitStr,oocyteProbability_copy,boost::is_any_of("*"));
	probability->parseString(splitStr[1]);
}

/*
 * Return interpolation method used for oocyte probabilities
 */
string ProgramOptionParser::oocyteProbabilityInterpolationMethod() {
	vector<string> splitStr;
	string oocyteProbability_copy(oocyteProbability); // boost split appears to modify oocyteProbability when splitting which behaves badly when multithreading, so make a local copy to split
	boost::split(splitStr,oocyteProbability_copy,boost::is_any_of("*"));
	#if defined __INCLUDESANITYCHECKS
	// sanity check
	if (strcmp(splitStr[0].c_str(),"linear")!=0 && strcmp(splitStr[0].c_str(),"step")!=0) {
		cerr << "Invalid choice of method in oocyteProbability" << endl;
		exit(1);
	}
	#endif
	return splitStr[0];
}

/*
 * Parse formatted fecundity schedule into 2D matrix
 */
void ProgramOptionParser::getFecundityConstraints(Array2D<float> *fecundityScheduleMatrix){
	fecundityScheduleMatrix->parseString(fecunditySchedule);
}

