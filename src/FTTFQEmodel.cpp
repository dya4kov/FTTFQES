#include "../hdr/FTTFQEmodel.h"
#include "../hdr/ThermodynamicFunction.h"
#include <iostream>
#include <string>

int main() {
	FTTFQEmodel m(1.0, 1.0);
	Int Dsize = 1;
	Int Tsize = 20;

	DensVec D(Dsize);
	TempVec T(Tsize);


	for (Int d = 0; d < Dsize; ++d) {
		D[d].setValue(3.0*(d + 1), gOverCmc);
	}
	for (Int t = 0; t < Tsize; ++t) {
		T[t].setValue(1.0*(t + 1), eV);
	}
	m.setParameters(D, T);
	m.setTolerance(1e-7);
	std::string input = "T[eV] P[GPa] DP[GPa] P0[GPa] M[eV] DM[eV] M0[eV] ";
	m.prepareOutput(input);
	m.printOutput("out/H_5gcc.dat");
	return 0;
}