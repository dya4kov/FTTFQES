#include "../hdr/FTTFmodel.h"
#include "../hdr/ThermodynamicFunction.h"
#include <iostream>
#include <string>

int main() {
	FTTFmodel m(13.0, 27.0);
	Int Dsize = 100;
	Int Tsize = 1;

	DensVec D(Dsize);
	TempVec T(Tsize);

	for (Int d = 0; d < Dsize; ++d) {
		D[d].setValue(0.1*(d + 1), gOverCmc);
	}
	for (Int t = 0; t < Tsize; ++t) {
		T[t].setValue(1.0*(t + 1));
	}
	m.setParameters(D, T);
	m.setTolerance(1e-7);
	std::string input = "T[K] P[GPa,lin] E TEoS";
	m.prepareOutput(input);
	m.printOutput("out/FTTFmodel.dat");
	return 0;
}