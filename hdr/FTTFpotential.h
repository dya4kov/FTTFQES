#pragma once
#include <fstream>
#include "../lib/numeric/hdr/ODE/ODEsolver.h"
#include "../lib/numeric/hdr/ODE/ODEdata.h"
#include "../lib/numeric/hdr/ODE/steppers/ODEstepperPD853.h"
#include "../lib/numeric/hdr/numericTypes.h"
#include "../lib/numeric/hdr/mathUtils.h"
#include "../lib/numeric/hdr/SF/FermiDirac.h"
#include "../hdr/ThermodynamicFunction.h"
#include "../hdr/Printer.h"
#include "../hdr/Timer.h"

/**
 * @brief This class implements interface for calculating Thomas-Fermi potential.
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
    * @brief Value of potential @f$ \phi(x) @f$.
    * @return @f$ \phi(x) @f$
    */
    Double operator() (const Double x);
	/**
    * @brief Value of derivative @f$ \phi'(x) @f$.
    * @return @f$ \phi'(x) @f$
    */
    Double derivative(const Double x);
	/**
	* @brief Set volume, temperature, and atomic number values for calculating potential.
	*/
    void setParameters(const Volume &V, const Temperature &T, const Double Z);
	/**
	* @brief Set tolerance for calculating potential.
	*/
    void setTolerance(Double eps);
    /**
	* @brief Print all steps for potential and its derivative into file 
    * from x = 1 to x = 0.
	*/
    void printData(const char* filename);
    /**
    * @brief Print N points for potential and its derivative into file
    * from x = 1 to x = 0.
    */
    void printData(const char* filename, const Int Npoints);
    /**
    * @brief Print all steps for potential and its derivative into file 
    * from x = 1 to x = 0.
    */
    DoubleMat getData();
    /**
    * @brief Print N points for potential and its derivative into file
    * from x = 1 to x = 0.
    */
    DoubleMat getData(const Int Npoints);
    /**
    * @brief Turn on writing log file.
    */
    void setLogOn();
    /**
    * @brief Turn off writing log file.
    */
    void setLogOff();
    /**
    * @brief Get output log as string.
    */
    std::string getLogOutput();
    /**
    * @brief Print log output into file.
    */
    void printLogOutput();
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
    void prepareOut(Int Npoints);

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
    Int precision;

    rhsFTTFpotential rhs;
    ODEsolver<ODEstepperPD853<rhsFTTFpotential> > solver;

    ODEdata calculatedData;

    bool calculated;

    Printer printer;
    Timer timer;
    bool logOn;
    std::stringstream LOG;
    std::stringstream OUT;
    DoubleMat data;
};

FTTFpotential::FTTFpotential(void) : 
                vTableSize(200), 
                tTableSize(175), 
                lgV0(-9.9), 
                lgT0(-7.4),
                lgVStep(0.1),
                lgTStep(0.1),
                phiTable(vTableSize, tTableSize),
                calculatedData(),
                solver(1e-6, 0.0),
                rhs(0.0)
{
    setPhiTable();
    solver.SetOutput(calculatedData);
    Vol = 0;
    Temp = 0;
    precision = 6;
    eps = 1e-6;
    logOn = false;
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
    Double error;
    Double shotParameter_1;
    Double shotParameter_2;
    Double xFrom = 1.0;
    Double xTo = 0.0;

    DoubleVec timeSteps(0.0, 100);
    Double timePerStep;
    Double overallTime;
    Double averageTime;
    Int nStep = 0;

    if (logOn) {
        LOG << "Begin calculation using shooting method." << std::endl;
    }

    setInitialShotParameters(shotParameter_1, shotParameter_2);
    if (abs(shotParameter_1 - shotParameter_2)/abs(shotParameter_1) < eps/10) {
        shotParameter_1 = shotParameter_1 - 1e-2*abs(shotParameter_1);
        shotParameter_2 = shotParameter_2 + 1e-2*abs(shotParameter_2);
    }
    
    if (logOn) {
        LOG << "Selected initial test values: " << std::endl;
        printer.printString(LOG, "\\phi_1^{left} = ", 24, right);
        printer.printSciDouble(LOG, shotParameter_1, precision, 15, left);
        LOG << std::endl;
        printer.printString(LOG, "\\phi_1^{right} = ", 24, right);
        printer.printSciDouble(LOG, shotParameter_2, precision, 15, left);
        LOG << std::endl;
        LOG << "Begin shooting:" << std::endl;
        printer.printString(LOG, "nStep",           10, left);
        printer.printString(LOG, "\\phi_1^{left}",  20, left);
        printer.printString(LOG, "\\phi_1^{right}", 20, left);
        printer.printString(LOG, "\\phi_1^{test}",  20, left);
        printer.printString(LOG, "\\phi_0^{test}",  20, left);
        printer.printString(LOG, "error",           20, left);
        printer.printString(LOG, "time[ms]",        20, left);
        LOG << std::endl;
    }

    error = abs(shotParameter_2 - shotParameter_1)/abs(shotParameter_1 + shotParameter_2);
    while (error > eps) {
        if (logOn) { // write LOG
            printer.printInt(LOG, nStep, 10, left);
            printer.printSciDouble(LOG, shotParameter_1, precision, 20, left);
            printer.printSciDouble(LOG, shotParameter_2, precision, 20, left);
            printer.printSciDouble(LOG, 0.5*(shotParameter_1 + shotParameter_2), precision, 20, left);
            timer.start();
        }

        phiCurrent[0] = 0.5*(shotParameter_1 + shotParameter_2);
        phiCurrent[1] = phiCurrent[0];
        solver.Integrate(rhs, phiCurrent, xFrom, xTo);
        if (phiCurrent[0] - phi_0 > 0)
            shotParameter_2 -= (shotParameter_2 - shotParameter_1)/2;
        else shotParameter_1 += (shotParameter_2 - shotParameter_1)/2;

        if (logOn) { // write LOG
            timer.stop();
            timePerStep = timer.getElapsedTimeInMilliSec();
            timeSteps[nStep] = timePerStep;
            printer.printSciDouble(LOG, phiCurrent[0], precision, 20, left);
            printer.printSciDouble(LOG, error, precision, 20, left);
            printer.printSciDouble(LOG, timePerStep, precision, 20, left);
            LOG << std::endl;
            timer.reset();
            ++nStep;
        }
        // update error
        error = abs(shotParameter_2 - shotParameter_1)/abs(shotParameter_1 + shotParameter_2);
    }
    phi_1 = (shotParameter_1 + shotParameter_2)/2;
    if (logOn) { // write last line in table, overall timing, and potential output
        timer.start();
        phiCurrent[0] = phi_1;
        phiCurrent[1] = phi_1;
        ODEdata outData(-1);
        solver.SetOutput(outData);
        solver.Integrate(rhs, phiCurrent, xFrom, xTo);
        timer.stop();
        timePerStep = timer.getElapsedTimeInMilliSec();
        timeSteps[nStep] = timePerStep;
        printer.printInt(LOG, nStep, 10, left);
        printer.printSciDouble(LOG, shotParameter_1, precision, 20, left);
        printer.printSciDouble(LOG, shotParameter_2, precision, 20, left);
        printer.printSciDouble(LOG, 0.5*(shotParameter_1 + shotParameter_2), precision, 20, left);
        printer.printSciDouble(LOG, phiCurrent[0], precision, 20, left);
        printer.printSciDouble(LOG, error, precision, 20, left);
        printer.printSciDouble(LOG, timePerStep, precision, 20, left);
        LOG << std::endl;
        LOG << "Finally selected \\phi_{1} = ";
        printer.printSciDouble(LOG, phi_1, precision, 20, left);
        LOG << std::endl;
        overallTime = 0.0;
        for (Int i = 0; i <= nStep; ++i) {
            overallTime += timeSteps[i];
        }
        averageTime = overallTime/(nStep + 1);
        LOG << "Overall time = ";
        printer.printSciDouble(LOG, overallTime, precision, 15, left);
        LOG << " ms" << std::endl;
        LOG << "Average time per step = ";
        printer.printSciDouble(LOG, averageTime, precision, 15, left);
        LOG << " ms" << std::endl;
        LOG << "output steps for potential:" << std::endl;
        Int size = outData.Count();
        printer.printString(LOG, "x");
        printer.printString(LOG, "phi");
        printer.printString(LOG, "dphi");
        LOG << std::endl;
        for (Int i = 0; i < size; ++i) {
            printer.printSciDouble(LOG, outData.xSave[i], 6);
            printer.printSciDouble(LOG, outData.ySave[0][i], 6);
            printer.printSciDouble(LOG, outData.ySave[1][i], 6);
            LOG << std::endl;
        }
        solver.SetOutput(calculatedData);
    }
    calculated = true;
}

Double FTTFpotential::operator() (const Double x) {
    if (!calculated) calculate();
    if (x < 0.0) { 
        std::cerr << "cannot calculate potential at x < 0" << std::endl; 
        exit(0);
    }
    if (x < 1.0) {
        DoubleVec phiCurrent(0.0, rhsFTTFpotential::dim);
        Double xFrom = 1.0;
        Double xTo = x;
        phiCurrent[0] = phi_1;
        phiCurrent[1] = phi_1;
        solver.Integrate(rhs, phiCurrent, xFrom, xTo);
        return phiCurrent[0];
    }
    else return phi_1;
}

Double FTTFpotential::derivative(const Double x) {
    if (!calculated) calculate();
    if (x < 0.0) { 
        std::cerr << "cannot calculate potential at x < 0" << std::endl; 
        exit(0);
    }
    if (x < 1.0) {
        DoubleVec phiCurrent(0.0, rhsFTTFpotential::dim);
        Double xFrom = 1.0;
        Double xTo = x;
        phiCurrent[0] = phi_1;
        phiCurrent[1] = phi_1;
        solver.Integrate(rhs, phiCurrent, xFrom, xTo);
        return phiCurrent[1];
    }
    else return phi_1;
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
        if (logOn) { // write LOG
            LOG << "FTTFpotential accepted new parameters:" << std::endl;
            printer.printString(LOG, "V_{Z}[Atomic] = ", 22, right);
            printer.printSciDouble(LOG, V(), precision, 15, left);
            LOG << std::endl;
            printer.printString(LOG, "T_{Z}[Hartree] = ", 23, right);
            printer.printSciDouble(LOG, T(), precision, 15, left);
            LOG << std::endl;
            printer.printString(LOG, "Z = ", 10, right);
            printer.printDouble(LOG, Z, 3, 15, left);
            LOG << std::endl;
            LOG << "Recalculate parameters:" << std::endl;
            printer.printString(LOG, "V_{1}[Atomic] = V_{Z}*Z = ", 32, right);
            printer.printSciDouble(LOG, Vol, precision, 15, left);
            LOG << std::endl;
            printer.printString(LOG, "T_{1}[Hartree] = T_{Z}*Z^{-4/3} = ", 40, right);
            printer.printSciDouble(LOG, Temp, precision, 15, left);
            LOG << std::endl;
            printer.printString(LOG, "a = 2^{7/6}*3^{2/3}*\\pi^{-5/3}*T^{1/2}*V^{2/3} = ", 55, right);
            printer.printSciDouble(LOG, a, precision, 15, left);
            LOG << std::endl;
            printer.printString(LOG, "\\phi_{0} = T^{-1}*(4*\\pi/(3V))^{1/3} = ", 45, right);
            printer.printSciDouble(LOG, phi_0, precision, 15, left);
            LOG << std::endl;
        }
    }
}

void FTTFpotential::printData(const char* filename) {
    prepareOut(-1);
    std::ofstream out;
    out.open(filename, std::ios::out);
    out << OUT.str();
}

void FTTFpotential::printData(const char* filename, const Int Npoints) {
    prepareOut(Npoints);
    std::ofstream out;
    out.open(filename, std::ios::out);
    out << OUT.str();
}

DoubleMat FTTFpotential::getData() {
    prepareOut(-1);
    return data;
}

DoubleMat FTTFpotential::getData(const Int Npoints) {
    prepareOut(Npoints);
    return data;
}   

void FTTFpotential::prepareOut(Int Npoints) {
    if (!calculated) calculate();
    // additional integration with saving steps
    Double xFrom = 1.0;
    Double xTo = 0.0;
    DoubleVec phiCurrent(0.0, rhsFTTFpotential::dim);
    phiCurrent[0] = phi_1;
    phiCurrent[1] = phi_1;
    ODEdata outData(Npoints);
    solver.SetOutput(outData);
    solver.Integrate(rhs, phiCurrent, xFrom, xTo);
    // write output data
    Int size = outData.Count();
    data = DoubleMat(0.0, rhsFTTFpotential::dim + 1, size);
    printer.printString(OUT, "x");
    printer.printString(OUT, "phi");
    printer.printString(OUT, "dphi");
    OUT << std::endl;
    for (Int i = 0; i < size; ++i) {
        printer.printSciDouble(OUT, data[0][i] = outData.xSave[i],    precision);
        printer.printSciDouble(OUT, data[1][i] = outData.ySave[0][i], precision);
        printer.printSciDouble(OUT, data[2][i] = outData.ySave[1][i], precision);
        OUT << std::endl;
    }
    solver.SetOutput(calculatedData);
}

void FTTFpotential::setTolerance(Double _eps) {
	eps = _eps;
    solver.SetTolerance(0.0, eps/10);
    if (logOn) { // write LOG
        LOG << "FTTFpotential accepted new tolerance:" << std::endl;
        printer.printString(LOG, "eps = ", 12, right);
        printer.printSciDouble(LOG, eps, precision, 15, left);
        LOG << std::endl;
    }
    precision = static_cast<int>(-log10(eps));
    calculated = false;
}

void FTTFpotential::setLogOn() {
    logOn = true;
}

void FTTFpotential::setLogOff() {
    logOn = false;
}

std::string FTTFpotential::getLogOutput() {
    return LOG.str();
}

void FTTFpotential::printLogOutput() {
    std::ofstream out;
    std::string filename = "out/phi_log_";
    filename += timer.getCurrentDatetime();
    filename += ".txt";
    out.open(filename.c_str(), std::ios::out);
    out << LOG.str();
}
