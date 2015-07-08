#include "../hdr/PeriodicTable.h"

int main() {
	PeriodicTable pt;
	std::cout << pt["Fe"].Z << std::endl;
	std::cout << pt["Fe"].mass << std::endl;
	pt.print();
	return 0;
}