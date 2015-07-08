#include "../hdr/QECorr.h"
#include "../hdr/ThermodynamicFunction.h"
#include "../hdr/InputReader.h"
#include "../hdr/PeriodicTable.h"
#include <string>
#include <iostream>
#include <stdlib.h>

int main(int argc, char **argv) {
	std::string inputfile = DEFAULT_INPUT_FILE;
	std::string outfile = DEFAULT_OUTPUT_FILE;
	int loglvl = 0; // no log files by default
	bool showprogress = true;
	int i = 1;
	while (i < argc) {
		std::string parameter;
		parameter = argv[i];
		if (!parameter.compare("--input")) {
			++i; inputfile = argv[i];
		}
		if (!parameter.compare("--output")) {
			++i; outfile = argv[i];
		}
		if (!parameter.compare("--loglvl")) {
			++i; loglvl = atoi(argv[i]);
		}
		if (!parameter.compare("--showprogress")) {
			++i; parameter = argv[i];
			if (!parameter.compare("off"))
				showprogress = false;
		}
		++i;
	}
	std::string element;
	std::string vol;
	std::string temp;
	std::string output;
	double tolerance;
	std::cout << "reading input file: in/" + inputfile << std::endl;
	readInput(inputfile, element, vol, temp, output, tolerance);
	PeriodicTable PT;
	QECorr qe(PT[element].Z, PT[element].mass);
	if (loglvl > 0) {
		qe.setPrintMainLogOn();
		if (loglvl > 1) {
			qe.setPrintPointLogOn();
		}
	}
	if (showprogress) qe.setShowProgressOn();
	qe.setParameters(vol, temp);
	qe.setTolerance(tolerance);
	qe.calculateData(output);
	std::cout << "writing output file: out/" + inputfile << std::endl;
	qe.printOutput("out/" + outfile);
	return 0;
}