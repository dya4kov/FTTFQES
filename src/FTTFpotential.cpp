#include "../hdr/FTTFpotential.h"
#include "../hdr/ThermodynamicFunction.h"
#include <iostream>

int main() {
	FTTFpotential phi;
	phi.setPrintLogOn();
	PhysQ V(1.0);
	PhysQ T(1.0);
	phi.setParameters(V, T, 1.0);
	phi.setTolerance(1e-6);
	std::cout << phi(1) << std::endl;
	std::cout << phi.derivative(1) << std::endl;
	phi.printData("out/phi.dat");
	phi.setParameters(V, T, 2.0);
	std::cout << phi(1) << std::endl;
	phi.printData("out/phi_new.dat");
	phi.setPrintLogOff();
	return 0;
}