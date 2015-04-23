#include "../hdr/FTTFpotential.h"
#include "../hdr/ThermodynamicFunction.h"
#include <iostream>

int main() {
	FTTFpotential phi;
	Volume V(1.0);
	Temperature T(1.0);
	phi.setParameters(V, T, 1.0);
	phi.setTolerance(1e-6);
	std::cout << phi.valueAt_1() << std::endl;
	std::cout << phi.derivativeAt_0() << std::endl;
	phi.printData("out/phi.dat");
	return 0;
}