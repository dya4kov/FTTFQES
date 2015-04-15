#pragma once
#include <fstream>
#include "../lib/numeric/hdr/ODE/ODEsolver.h"
#include "../lib/numeric/hdr/ODE/ODEdata.h"
#include "../lib/numeric/hdr/ODE/steppers/ODEstepperPD853.h"
#include "../lib/numeric/hdr/numericTypes.h"
#include "../lib/numeric/hdr/mathUtils.h"
#include "../lib/numeric/hdr/SF/FermiDirac.h"
#include "../hdr/ThermodynamicFunction.h"

/**
 * @brief This class implements Interface for calculating Thomas-Fermi potential.
 */
class FTTFpotential {
public:
   /**
	* @brief Constructor of Thomas-Fermi potential.
	* @details Allocates memory for potential solver
	* and loads table of precalculated values of potential at
	* @f$ x = 1 @f$.
	*/
    FTTFpotential(void);
   /**
	* @brief Value of potential @f$ \phi @f$ at @f$ x = 1 @f$.
	* @details Solving boundary problem by shooting method with 
	* good start values of @f$ \phi(1) @f$ from table for fast 
	* convergency.
	* @return @f$ \phi(1) @f$
	*/
    Double valueAt_1();
	/**
	* @brief Value of potential @f$ \phi @f$ at some poInt @f$ x @f$.
	* @return @f$ \phi(x) @f$
	*/
    Double valueAt_x(Double x);
	/**
	* @brief Value of derivative of potential @f$ \phi'(0) @f$.
	* @return @f$ \phi'(0) @f$
	*/
    Double derivativeAt_0();
	/**
	* @brief Set volume and temperature values for calculating potential.
	*/
    void setParameters(const Volume &V, const Temperature &T, const Double Z);
	/**
	* @brief Set tolerance for calculating potential.
	*/
    void setTolerance(Double eps);
    /**
	* @brief Print output data.
	*/
    void printData(const char* filename);
    /**
    * @brief Boundary problem for the Thomas-Fermi potential calculation
    * @details Thomas-Fermi potential is calculated from the following boundary problem:
    * @f{eqnarray*}{
    *  \frac{d^2\phi}{dx^2} &=& a x I_{1/2}
    *	  \left(
    *		\frac{\phi(x)}{x}
    *	  \right), \\
    *  \phi(0) &=& \frac{1}{T}\sqrt{\frac{4\pi}{3V}}, \\
    *  \phi(1) &=& \phi'(1).
    * @f}
    * Here @f$ a = 2^{7/6}3^{2/3}\pi^{-5/3}T^{1/2}V^{2/3} @f$.
    */
    struct rhsFTTFpotential {
		Double a;
		static const Int dim = 2;
		FermiDirac<Half> FDhalf;
		rhsFTTFpotential(Double _a) : a(_a) {}
		void updateParameter(const Double _a) { a = _a; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
			dydx[0] = y[1];
    		if (x > 0) dydx[1] = a*x*FDhalf(y[0]/x);
    		else dydx[1] = 1e+10;
		}
	};

private:
    void setPhiTable();
    void calculate();
    void setInitialShotParameters(Double& p1, Double& p2);

    const Int vTableSize;
    const Int tTableSize;
    const Double lgV0;
    const Double lgT0;
    const Double lgVStep;
    const Double lgTStep;
    DoubleMat phiTable;
    
    Double phi_0;
    Double phi_1;
    Double Vol;
    Double Temp;
    Double eps;

    rhsFTTFpotential rhs;
    ODEsolver<ODEstepperPD853<rhsFTTFpotential> > solver;

    ODEdata calculatedData;

    bool calculated;
};

FTTFpotential::FTTFpotential(void) : 
                vTableSize(200), 
                tTableSize(175), 
                lgV0(-9.9), 
                lgT0(-7.4),
                lgVStep(0.1),
                lgTStep(0.1),
                phiTable(vTableSize, tTableSize),
                calculatedData(-1), // save all steps
                solver(1e-6, 0.0),
                rhs(0.0)
{
    setPhiTable();
    Vol = 0;
    Temp = 0;
}

void FTTFpotential::setPhiTable() {
    std::ifstream data("res/TFPotential.dat", std::ios::in);
    Double currentValue;
    for (Int t = 0; t < tTableSize; ++t) {
        for (Int v = 0; v < vTableSize; ++v) {
            data >> currentValue;
            phiTable[v][t] = currentValue;
        }
    }
    data.close();
}

void FTTFpotential::setInitialShotParameters(Double& p1, Double& p2) {
    Double phi1_max;
    Double phi1_min;
    Int v;
    Int t;
    bool v_is_calculated = false;
    bool t_is_calculated = false;

    v = 0;
    t = 0;

    v = (Int) ceil((log10(Vol) - lgV0)*10);
    if (abs(log10(Vol) - lgV0 - v*lgVStep) < 1e-10)
    {
        v_is_calculated = true;
    }
    else if (abs(log10(Vol) - lgV0 - (v - 1)*lgVStep) < 1e-10)
    {
        v_is_calculated = true;
        v--;
    }

    t = (Int) ceil((log10(Temp) - lgT0)*10);
    if (abs(log10(Temp) - lgT0 - t*lgTStep) < 1e-10)
    {
        t_is_calculated = true;
    }
    else if (abs(log10(Temp) - lgT0 - (t - 1)*lgTStep) < 1e-10)
    {
        t_is_calculated = true;
        t--;
    }

    if (v_is_calculated && t_is_calculated)
    {
        p1 = p2 = phiTable[v][t];
    }
    else 
    {
        if (v == 0 || v_is_calculated) { 
            phi1_min = phiTable[v][t];
            phi1_max = phiTable[v][t - 1];
        }
        else if (t == 0 || t_is_calculated) {
            phi1_min = phiTable[v][t];
            phi1_max = phiTable[v - 1][t];
        }
        else {
            phi1_min = phiTable[v][t];
            phi1_max = phiTable[v - 1][t - 1];
        }
        p1 = phi1_min;
        p2 = phi1_max;
    }
}

void FTTFpotential::calculate() {
	DoubleVec phiCurrent(0.0, rhsFTTFpotential::dim);
    Double shotParameter_1;
    Double shotParameter_2;
    Double xFrom = 1.0;
    Double xTo = 0.0;
    setInitialShotParameters(shotParameter_1, shotParameter_2);

    if (abs(shotParameter_1 - shotParameter_2)/abs(shotParameter_1) < eps/10) {
        shotParameter_1 = shotParameter_1 - 1e-2*abs(shotParameter_1);
        shotParameter_2 = shotParameter_1 + 1e-2*abs(shotParameter_2);
    }

    solver.SetOutput(calculatedData);
    solver.SetTolerance(0.0, eps/10);
    
    while (abs(shotParameter_2 - shotParameter_1)/abs(shotParameter_1) > eps) {
        phiCurrent[0] = (shotParameter_1 + shotParameter_2)/2;
        phiCurrent[1] = phiCurrent[0];
        solver.Integrate(rhs, phiCurrent, xFrom, xTo);
        if (phiCurrent[0] - phi_0 > 0)
            shotParameter_2 -= (shotParameter_2 - shotParameter_1)/2;
        else shotParameter_1 += (shotParameter_2 - shotParameter_1)/2;
    }
    phi_1 = (shotParameter_1 + shotParameter_2)/2;
    calculated = true;
}

Double FTTFpotential::valueAt_1() {
    if (!calculated) calculate();
    return phi_1;
}

Double FTTFpotential::valueAt_x(Double x) {
    // if (!calculated) calculate();
    // Double initials[2] = { phi_1, phi_1 };
    Double result = 0;
    // solver->setInitials(initials);
    // solver->setLimits(1.0, x);
    // result = solver->driverApply().y[0];
    // solver->reset();
    return result;
}

Double FTTFpotential::derivativeAt_0() {
    if (!calculated) calculate();
    Int lastPoint = calculatedData.Count() - 1;
    return calculatedData.ySave[1][lastPoint];
}

void FTTFpotential::setParameters(const Volume &V, const Temperature &T, const Double Z) {
    Volume V1;
    Temperature T1;  
    V1.setValue(V()*Z);
    T1.setValue(T()*pow(Z, -4.0/3.0));
    if (abs(log10(V1()) - log10(Vol)) > 1e-10
        || abs(log10(T1()) - log10(Temp)) > 1e-10) 
    { 
        calculated = false;
        Vol = V1();
        Temp = T1();
        Double a =   pow(2.0, 7.0/6.0)
                   * pow(3.0, 2.0/3.0)
                   * pow(M_PI, -5.0/3.0)
                   * sqrt(Temp)*pow(Vol, 2.0/3.0);
        rhs.updateParameter(a);
        phi_0 = pow(4.0*M_PI/3.0/Vol, 1.0/3.0)/Temp;
    }
}

void FTTFpotential::printData(const char* filename) {
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
    out << "phi";
    out.setf(std::ios::left);
    out.width(20);
    out.fill(' ');
    out << "dphi";
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
        out.setf(std::ios::scientific | std::ios::left | std::ios::showpos);
	    out.width(20);
	    out.precision(digitsToPrint);
	    out.fill(' ');
        out << calculatedData.ySave[1][i];
    	out << std::endl;
    }
}

void FTTFpotential::setTolerance(Double _eps) {
	eps = _eps;
}