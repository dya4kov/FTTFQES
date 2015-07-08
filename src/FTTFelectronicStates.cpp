#include "../hdr/FTTFelectronicStates.h"
#include "../hdr/ThermodynamicFunction.h"
#include <iostream>

int main() {
	FTTFelectronicStates ES(7);
	ES.setPrintMainLogOn();
	PhysQ V(1.2);
	PhysQ T(100.2);
	ES.setParameters(V, T, 13.0);
	std::cout << ES.Nlow() << std::endl;
	std::cout << ES.NTFlow() << std::endl;
	std::cout << ES.DNlow() << std::endl;
	return 0;
}