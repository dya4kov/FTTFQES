#include "../lib/numeric/hdr/ODE/ODEsolver.h"
#include "../lib/numeric/hdr/ODE/steppers/ODEstepperPD853.h"
#include "../lib/numeric/hdr/numericTypes.h"
#include <iostream>

struct rhs_van {
	Double eps;
	static const Int dim = 2;
	rhs_van(Double epss) : eps(epss) {}
	void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
		dydx[0] = y[1];
		dydx[1] = ((1.0-y[0]*y[0])*y[1]-y[0])/eps;
	}
};

int main() {
	const Double atol = 1.0e-3, rtol = atol, h1 = 0.01, hmin = 0.0, x1 = 0.0, x2 = 2.0;
	DoubleVec ystart(rhs_van::dim);
	ystart[0] = 2.0;
	ystart[1] = 0.0;
	ODEdata out(40);
	rhs_van d(1.0e-3);
	ODEsolver<ODEstepperPD853<rhs_van> > ode(h1,hmin);
	ode.SetTolerance(atol, rtol);
	ode.SetOutput(out);
	ode.Integrate(d, ystart, x1, x2);
	for (Int i = 0; i < out.Count(); ++i) {
		std::cout << out.xSave[i] << " " << out.ySave[0][i] << " " <<
			out.ySave[1][i] << std::endl;
	}
	return 0;
}