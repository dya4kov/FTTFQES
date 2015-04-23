#include "../hdr/FTTFelectronicStates.h"
#include "../hdr/ThermodynamicFunction.h"
#include <iostream>

int main() {
	FTTFelectronicStates ES(7);
	Volume V(1.0);
	Temperature T(1.0);
	ES.setParameters(V, T, 1.0);
	std::cout << ES.getDN() << std::endl;
	return 0;
}