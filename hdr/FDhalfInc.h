#pragma once
#include "../lib/numeric/hdr/ODE/ODEsolver.h"
#include "../lib/numeric/hdr/ODE/ODEdata.h"
#include "../lib/numeric/hdr/ODE/steppers/ODEstepperPD853.h"

class FDhalfInc {
public:
	struct RHS {
		static const Int dim = 1;
		Double p;
		RHS() {}
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
			dydx[0] = sqrt(x)/(1 + exp(x - p));
		}
	};
	Double operator() (const Double x, const Double from, const Double to); 
	FDhalfInc(const Double tolerance) : data(-1), // save all steps
				  solver(1e-7),
				  rhs() 
	{
		solver.SetOutput(data);
        solver.SetTolerance(tolerance, tolerance);
	}
private:
    RHS rhs;
    ODEsolver<ODEstepperPD853<RHS> > solver;
    ODEdata data;
};

Double FDhalfInc::operator() (const Double x, const Double from, const Double to) {
	DoubleVec y0(0.0, 1);
	Double xFrom;
	Double xTo;
	rhs.p = x;
	y0[0] = 0;
	xFrom = from > 0 ? from : 0;
	xTo = to;
	solver.Integrate(rhs, y0, xFrom, xTo);
	return y0[0];
}