#include "../hdr/Yfunction.h"
#include "../lib/numeric/hdr/numericTypes.h"
#include <iostream>

int main() {
	Yfunction Y(1e-10);
	std::cout.setf(std::ios::scientific | std::ios::left | std::ios::showpos);
    std::cout.width(20);
    std::cout.precision(20);
    std::cout.fill(' ');
	std::cout << "Y(0) = " << Y(100.0) << std::endl;
	Y.printData("./out/Y.dat");
	return 0;
}