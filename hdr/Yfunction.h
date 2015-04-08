#pragma once
#include <fstream>
#include "../lib/numeric/hdr/ODE/ODEsolver.h"
#include "../lib/numeric/hdr/ODE/ODEdata.h"
#include "../lib/numeric/hdr/ODE/steppers/ODEstepperPD853.h"
#include "../lib/numeric/hdr/numericTypes.h"
#include "../lib/numeric/hdr/mathUtils.h"
#include "../lib/numeric/hdr/SF/FermiDirac.h"

class Yfunction {
public:
	Yfunction(Double tolerance = 1e-7) : 
				  calculatedData(-1), // save all steps
                  solver(1e-6, 0.0) {
      	solver.SetOutput(calculatedData);
		solver.SetTolerance(0.0, tolerance);
    }
	Double operator() (const Double x) {
		return 0.5*FDhalf(x)*FDmhalf(x) + 1.5*Integral(x);
	}
	Double derivative(const Double x) {
		static Double f;
		f = FDmhalf(x);
		return 1.75*f*f + 0.5*FDhalf(x)*FDdmhalf(x);
	}
	void printData(const char* filename) {
		std::ofstream out;
	    Int digitsToPrint = 6;
	    out.open(filename, std::ios::out);
	    Int size = calculatedData.Count();
	    out.setf(std::ios::left);
	    out.width(20);
	    out.fill(' ');
	    out << "x";
	    out.setf(std::ios::left);
	    out.width(20);
	    out.fill(' ');
	    out << "int(I^2_{-1/2}(x))";
	    out << std::endl;
	    for (Int i = 0; i < size; ++i) {
		    out.setf(std::ios::scientific | std::ios::left | std::ios::showpos);
		    out.width(20);
		    out.precision(digitsToPrint);
		    out.fill(' ');
	        out << calculatedData.xSave[i];
	        out.setf(std::ios::scientific | std::ios::left | std::ios::showpos);
		    out.width(20);
		    out.precision(digitsToPrint);
		    out.fill(' ');
	        out << calculatedData.ySave[0][i];
	        out << std::endl;
    	}
	}
	struct rhsYfunc {
		static const Int dim = 1;
		FermiDirac<Mhalf> FDmhalf;
		rhsYfunc() {}
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
			static Double f;
			f = FDmhalf(x);
			dydx[0] = f*f;
		}
	};
private:
	Double Integral(const Double x) {
		DoubleVec y(0.0, rhsYfunc::dim);
    	solver.Integrate(rhs, y, -100.0, x);
    	return y[0];
	}

	FermiDirac<DMhalf> FDdmhalf;
	FermiDirac<Half> FDhalf;
	FermiDirac<Mhalf> FDmhalf;

	rhsYfunc rhs;
    ODEsolver<ODEstepperPD853<rhsYfunc> > solver;
    ODEdata calculatedData;
};