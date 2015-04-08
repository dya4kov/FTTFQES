#pragma once
#include <fstream>
#include <iostream>
#include "../lib/numeric/hdr/ODE/ODEsolver.h"
#include "../lib/numeric/hdr/ODE/ODEdata.h"
#include "../lib/numeric/hdr/ODE/steppers/ODEstepperPD853.h"
#include "../lib/numeric/hdr/numericTypes.h"
#include "../lib/numeric/hdr/mathUtils.h"
#include "../lib/numeric/hdr/SF/FermiDirac.h"
#include "../hdr/ThermodynamicFunction.h"
#include "../hdr/FTTFpotential.h"
#include "../hdr/Yfunction.h"
/**
* @brief The class implements Interface for calculating correction to the Thomas-Fermi potential.
*/
class FTTFQEpotential
{
public:
	/**
	* @brief Constructor of correction to the Thomas-Fermi potential.
	* @details Allocates memory for potential correction solver
	* and loads table of precalculated values of potential correction at
	* @f$ x = 1 @f$.
	*/
    FTTFQEpotential(void);
	/**
	* @brief Value of correction to potential @f$ \psi @f$ at @f$ x = 1 @f$.
	* @details Solving boundary problem by shooting method with
	* good start values of @f$ \psi(1) @f$ from table for fast
	* convergency.
	* @return @f$ \psi(1) @f$
	*/
    Double valueAt_1();
	/**
	* @brief Value of correction to potential @f$ \psi @f$ at some poInt @f$ x @f$.
	* @return @f$ \psi(x) @f$
	*/
    Double valueAt_x(Double x);
	/**
	* @brief Value of derivative of correction to potential @f$ \psi'(0) @f$.
	* @return @f$ \psi'(0) @f$
	*/
    Double derivativeAt_0();
	/**
	* @brief Set volume and temperature values for calculating correction to potential.
	*/
    void setParameters(const Volume& V, const Temperature &T);
	/**
	* @brief Set precision for calculating correction to potential.
	*/
    void setTolerance(Double eps);
    /**
	* @brief Print output data.
	*/
    void printData(const char* filename);
    /**
    * @brief Boundary problem for calculating the quantum and exchange correction to the Thomas-Fermi potential
    * @details correction to the Thomas-Fermi potential is calculated from the following boundary problem:
	* @f{eqnarray*}{
	*  \frac{d^2\phi}{dx^2} &=& a x I_{1/2}
	*	  \left(
	*		\frac{\phi(x)}{x}
	*	  \right), \\
	*  \frac{d^2\psi}{dx^2} &=& a \left[
	*		       \frac12 \psi(x) I_{-1/2}\left(\frac{\phi(x)}{x}\right) + xY'(x)
	*           \right], \\
	*  \phi(1) &=& \phi'(1) = \phi_1, \\
	*  \psi(1) &=& \psi'(1), \\
	*  \psi(0) &=& 0.
	* @f}
	* Here @f$ a = 2^{7/6}3^{2/3}\pi^{-5/3}T^{1/2}V^{2/3} @f$, @f$ \phi @f$ - Thomas-Fermi potential, 
	* @f$ \psi @f$ - correction to Thomas-Fermi potential.
    */
	struct rhsFTTFQEpotential {
		Double a;
		static const Int dim = 4;
		FermiDirac<Half> FDhalf;
		FermiDirac<Mhalf> FDmhalf;
		Yfunction Y;
		rhsFTTFQEpotential(Double _a) : a(_a) {}
		void updateParameter(const Double _a) { a = _a; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
			dydx[0] = y[1];
		    dydx[2] = y[3];
		    if (x > 0) {
		        dydx[1] = a*x*FDhalf(y[0]/x);
		        dydx[3] = a*(0.5*FDmhalf(y[0]/x)*y[2] + x*Y.derivative(y[0]/x));
		    }
		    else {
		        dydx[1] = dydx[3] = 1e+10;
		    }
		}
	};

private:

    void setPsiTable();
    void calculate();
    void setInitialShotParameters(Double& p1, Double& p2);

    const Int vTableSize;
    const Int tTableSize;
    const Double lgV0;
    const Double lgT0;
    const Double lgVStep;
    const Double lgTStep;
    DoubleMat psiTable;

    Double psi_1;
    Double Vol;
    Double Temp;
    Double eps;

    rhsFTTFQEpotential rhs;
    ODEsolver<ODEstepperPD853<rhsFTTFQEpotential> > solver;

    ODEdata calculatedData;

    FTTFpotential TFpot;

    bool calculated;
};

FTTFQEpotential::FTTFQEpotential(void) :
                        vTableSize(191), 
                        tTableSize(171), 
                        lgV0(-9.0), 
                        lgT0(-7.0),
                        lgVStep(0.1),
                        lgTStep(0.1),
                        psiTable(vTableSize, tTableSize),
                		calculatedData(-1), // save all steps
                		solver(1e-6, 0.0),
                		rhs(0.0)
{
    setPsiTable();
    Vol = 0;
    Temp = 0;
}

void FTTFQEpotential::setPsiTable() {
    std::ifstream data("res/TFPotentialCorrection.dat", std::ios::in);
    Double currentValue;
    for (Int t = 0; t < tTableSize; ++t) {
        for (Int v = 0; v < vTableSize; ++v) {
            data >> currentValue;
            psiTable[v][t] = currentValue;
        }
    }
    data.close();
}

void FTTFQEpotential::setInitialShotParameters(Double& p1, Double& p2) {
    Double psi1_min = 0;
    Double psi1_max = 0;
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
        p1 = p2 = psiTable[v][t];
    }
    else 
    {
        if (v == 0 || v_is_calculated) { 
            psi1_max = psiTable[v][t];
            psi1_min = psiTable[v][t - 1];
        }
        else if (t == 0 || t_is_calculated) {
            psi1_max = psiTable[v][t];
            psi1_min = psiTable[v - 1][t];
        }
        else {
            psi1_max = psiTable[v][t];
            psi1_min = psiTable[v - 1][t - 1];
        }
        p1 = psi1_min;
        p2 = psi1_max;
    }
}

void FTTFQEpotential::calculate() {
    DoubleVec psiCurrent(0.0, rhsFTTFQEpotential::dim);
    Double shotParameter_1;
    Double shotParameter_2;
    Double xFrom = 1.0;
    Double xTo = 0.0;
    setInitialShotParameters(shotParameter_1, shotParameter_2);

    if (abs(shotParameter_1 - shotParameter_2)/abs(shotParameter_1) < eps/10) {
        shotParameter_1 = shotParameter_1 - 1e-2*abs(shotParameter_1);
        shotParameter_2 = shotParameter_2 + 1e-2*abs(shotParameter_2);
    }

    solver.SetOutput(calculatedData);
    solver.SetTolerance(0.0, eps/10);

    Int nStep = 0;

    while (abs(shotParameter_2 - shotParameter_1)/abs(shotParameter_1) > eps) {
        psiCurrent[0] = TFpot.valueAt_1();
	    psiCurrent[1] = psiCurrent[0];
        psiCurrent[2] = (shotParameter_1 + shotParameter_2)/2;
        psiCurrent[3] = psiCurrent[2];
        solver.Integrate(rhs, psiCurrent, xFrom, xTo);
        if (psiCurrent[2] > 0)
            shotParameter_2 -= (shotParameter_2 - shotParameter_1)/2;
        else shotParameter_1 += (shotParameter_2 - shotParameter_1)/2;
    }
    psi_1 = (shotParameter_1 + shotParameter_2)/2;
    calculated = true;
}

Double FTTFQEpotential::valueAt_1() {
    if (!calculated) calculate();
    return psi_1;
}

Double FTTFQEpotential::valueAt_x(Double x) {
    // if (!calculated) calculate();
    // Double initials[2] = { psi_1, psi_1 };
    Double result = 0.0;
    // solver->setInitials(initials);
    // solver->setLimits(1.0, x);
    // result = solver->driverApply().y[0];
    // solver->reset();
    return result;
}

Double FTTFQEpotential::derivativeAt_0() {
    if (!calculated) calculate();
    Int lastPoint = calculatedData.Count() - 1;
    return calculatedData.ySave[3][lastPoint];
}

void FTTFQEpotential::setParameters(const Volume& V, const Temperature& T) {
	if (abs(log10(V()) - log10(Vol)) > 1e-10
        || abs(log10(T()) - log10(Temp)) > 1e-10) 
    { 
        calculated = false;
        Vol = V();
        Temp = T();
        Double a =   pow(2.0, 7.0/6.0)
                   * pow(3.0, 2.0/3.0)
                   * pow(M_PI, -5.0/3.0)
                   * sqrt(Temp)*pow(Vol, 2.0/3.0);
        rhs.updateParameter(a);
        TFpot.setParameters(V, T);
    }
}

void FTTFQEpotential::setTolerance(Double _eps) {
	TFpot.setTolerance(_eps);
    eps = _eps;
}

void FTTFQEpotential::printData(const char* filename) {
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
    out.setf(std::ios::left);
    out.width(20);
    out.fill(' ');
    out << "psi";
    out.setf(std::ios::left);
    out.width(20);
    out.fill(' ');
    out << "dpsi";
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
        out.setf(std::ios::scientific | std::ios::left | std::ios::showpos);
	    out.width(20);
	    out.precision(digitsToPrint);
	    out.fill(' ');
        out << calculatedData.ySave[2][i];
        out.setf(std::ios::scientific | std::ios::left | std::ios::showpos);
	    out.width(20);
	    out.precision(digitsToPrint);
	    out.fill(' ');
        out << calculatedData.ySave[3][i];
    	out << std::endl;
    }
}