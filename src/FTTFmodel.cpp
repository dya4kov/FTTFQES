#include "../hdr/FTTF.h"
#include "../hdr/ThermodynamicFunction.h"
#include <iostream>
#include <string>

int main() {
	FTTF m(13.0, 27.0);
	std::string vol = "D[g/cm^3,lin](1.0, 3.0, 3)";
	std::string temp = "T[atm,lin](1.0, 3.0, 3)";
	m.setPrintMainLogOn();
	m.setShowProgressOn();
	//m.setPrintPointLogOn();
	m.setParameters(vol, temp);
	m.setTolerance(1e-7);
	std::string input = "D[g/cm^3] T[eV] P[GPa,lin] E";
	m.calculateData(input);
	m.printOutput("out/FTTFmodel.dat");
	return 0;
}