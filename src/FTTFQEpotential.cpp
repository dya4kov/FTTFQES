#include "../hdr/FTTFQEpotential.h"
#include "../hdr/ThermodynamicFunction.h"
#include <iostream>

int main() {
	FTTFQEpotential psi;
	Volume V(1.0);
	Temperature T(1.0);
	psi.setParameters(V, T);
	psi.setTolerance(1e-6);
	std::cout << psi.valueAt_1() << std::endl;
	// std::cout << psi.derivativeAt_0() << std::endl;
	psi.printData("out/psi.dat");
	return 0;
}