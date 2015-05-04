#include "../hdr/Printer.h"
#include <string>

int main() {
	Printer p;
	std::stringstream ss;
	p.printString(ss, "abcde", 15, right, '=');
	std::string s;
	s = ss.str();
	std::cout << s;
	std::cout << std::endl;
}