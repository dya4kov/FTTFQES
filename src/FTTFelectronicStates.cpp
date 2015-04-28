#include "../hdr/FTTFelectronicStates.h"
#include "../hdr/ThermodynamicFunction.h"
#include <iostream>

int main() {
	FTTFelectronicStates ES(7);
	Volume V(1.2);
	Temperature T(100.2);
	ES.setParameters(V, T, 13.0);
	std::cout << ES.N() << std::endl;
	std::cout << ES.NTF() << std::endl;
	std::cout << ES.DN() << std::endl;
	return 0;
}