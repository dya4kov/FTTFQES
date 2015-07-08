#include "../hdr/FTTFQEpotential.h"
#include "../hdr/ThermodynamicFunction.h"
#include <iostream>
#include <stdlib.h>
#include "../hdr/stringUtils.h"

int main(int argc, char** argv) {
	FTTFQEpotential psi;
	double vol = 1.0;
	double temp = 1.0;
	double Z = 1.0;
	double eps = 1e-6;
	std::string outfile = "out/psi(";
	std::string params = "";
	int i = 1;
	while (i < argc) {
		std::string parameter;
		parameter = argv[i];
		if (!parameter.compare("--output")) {
			++i; outfile = argv[i];
			outfile = split(outfile, '.')[0];
			outfile = "out/" + outfile + "(";
		}
		if (!parameter.compare("--log")) {
			++i; parameter = argv[i];
			if (!parameter.compare("on"))
				psi.setPrintLogOn();
		}
		if (!parameter.compare("--V")) {
			++i; vol = atof(argv[i]);
			params = params + "_V=" + argv[i] + "_";
		}
		if (!parameter.compare("--T")) {
			++i; temp = atof(argv[i]);
			params = params + "_T=" + argv[i] + "_";
		}
		if (!parameter.compare("--Z")) {
			++i; Z = atof(argv[i]);
			params = params + "_Z=" + argv[i] + "_";
		}
		if (!parameter.compare("--eps")) {
			++i; eps = atof(argv[i]);
			params = params + "_eps=" + argv[i] + "_";
		}
		++i;
	}
	outfile = outfile + params + ").dat";
	PhysQ V(vol);
	PhysQ T(temp);
	psi.setParameters(V, T, Z);
	psi.setTolerance(eps);
	std::cout << "psi(1) = " << psi(1) << std::endl;
	std::cout << "psi'(1) = " << psi.derivative(1) << std::endl;
	psi.printData(outfile.c_str());
	return 0;
}