#include "../hdr/QECorr.h"
#include "../hdr/ThermodynamicFunction.h"
#include <iostream>
#include <string>

int main() {
	QECorr qe(13.0, 27.0);
	std::string vol = "D[g/cm^3,lin](1.0, 3.0, 3)";
	std::string temp = "T[atm,lin](1.0, 3.0, 3)";
	qe.setPrintMainLogOn();
	//qe.setShowProgressOn();
	qe.setPrintPointLogOn();
	qe.setParameters(vol, temp);
	qe.setTolerance(1e-7);
	std::string input = "D[g/cm^3] T[eV] DP[GPa,lin] DET";
	qe.calculateData(input);
	qe.printOutput("out/QECorr.dat");
	return 0;
}