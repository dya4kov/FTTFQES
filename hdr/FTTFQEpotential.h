#pragma once
#include <fstream>
#include <sstream>
#include "../lib/numeric/hdr/ODE/ODEsolver.h"
#include "../lib/numeric/hdr/ODE/ODEdata.h"
#include "../lib/numeric/hdr/ODE/steppers/ODEstepperPD853.h"
#include "../lib/numeric/hdr/numericTypes.h"
#include "../lib/numeric/hdr/mathUtils.h"
#include "../lib/numeric/hdr/SF/FermiDirac.h"
#include "../hdr/ThermodynamicFunction.h"
#include "../hdr/FTTFpotential.h"
#include "../hdr/Printer.h"
#include "../hdr/Timer.h"
#include "../hdr/Yfunction.h"

#define POTENTIAL_CORR_DATA_FILE "res/TFPotentialCorrection.dat"

/**
* @brief The class implements interface for calculating correction to the Thomas-Fermi potential.
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
    FTTFQEpotential();
    ~FTTFQEpotential();
	/**
	* @brief Value of correction to potential @f$ \psi(x) @f$.
	* @return @f$ \psi(x) @f$
	*/
    Double operator() (const Double x);
    /**
    * @brief Value of derivative of correction to potential @f$ \psi'(x) @f$.
    * @return @f$ \psi'(x) @f$
    */
    Double derivative(const Double x);
	/**
	* @brief Set volume and temperature values for calculating correction to potential.
	*/
    void setParameters(const PhysQ& V, const PhysQ &T, const Double Z);
	/**
	* @brief Set precision for calculating correction to potential.
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
    * @brief Write log to specified stream.
    */
    void setLogStream(std::ofstream* _LOG);
    /**
    * @brief Disable writing log to specified stream.
    */
    void clearLogStream();
    /**
    * @brief Enable self-printing log output into file.
    */
    void setPrintLogOn();
    /**
    * @brief Disable self-printing log output into file.
    */
    void setPrintLogOff();
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
    void performOutput(Int Npoints);

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
    Int precision;

    rhsFTTFQEpotential rhs;
    ODEsolver<ODEstepperPD853<rhsFTTFQEpotential> > solver;

    ODEdata calculatedData;

    FTTFpotential phi;

    bool calculated;

    Printer printer;
    Timer timer;
    bool printLogOn;
    bool logStreamIsSet;
    std::ofstream* LOG;
    std::ofstream* OUT;
    DoubleMat data;
};

FTTFQEpotential::FTTFQEpotential(void) :
                        vTableSize(191), 
                        tTableSize(171), 
                        lgV0(-9.0), 
                        lgT0(-7.0),
                        lgVStep(0.1),
                        lgTStep(0.1),
                        psiTable(vTableSize, tTableSize),
                		calculatedData(), // save all steps
                		solver(1e-6, 0.0),
                		rhs(0.0)
{
    setPsiTable();
    solver.SetOutput(calculatedData);
    precision = 6;
    eps = 1e-6;
    logStreamIsSet = false;
    printLogOn = false;
    Vol = 0;
    Temp = 0;
    OUT = NULL;
    LOG = NULL;
}

FTTFQEpotential::~FTTFQEpotential() {
    if (printLogOn) setPrintLogOff();
}

void FTTFQEpotential::setPsiTable() {
    std::ifstream datafile("res/TFPotentialCorrection.dat", std::ios::in);
    Double currentValue;
    for (Int t = 0; t < tTableSize; ++t) {
        for (Int v = 0; v < vTableSize; ++v) {
            datafile >> currentValue;
            psiTable[v][t] = currentValue;
        }
    }
    datafile.close();
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
    Double error;
    Double shotParameter_1;
    Double shotParameter_2;
    Double xFrom = 1.0;
    Double xTo = 0.0;
    Double phi_1 = phi(1);

    DoubleVec timeSteps(0.0, 100);
    Double timePerStep;
    Double overallTime;
    Double averageTime;
    Int nStep = 0;

    if (logStreamIsSet || printLogOn) {
        (*LOG) << "Begin calculation of correction using shooting method." << std::endl;
    }

    setInitialShotParameters(shotParameter_1, shotParameter_2);
    if (abs(shotParameter_1 - shotParameter_2)/abs(shotParameter_1) < eps/10) {
        shotParameter_1 = shotParameter_1 - 1e-2*abs(shotParameter_1);
        shotParameter_2 = shotParameter_2 + 1e-2*abs(shotParameter_2);
    }

    if (logStreamIsSet || printLogOn) {
        (*LOG) << "Selected initial test values: " << std::endl;
        printer.printString((*LOG), "\\psi_1^{left} = ", 24, right);
        printer.printSciDouble((*LOG), shotParameter_1, precision, 15, left);
        (*LOG) << std::endl;
        printer.printString((*LOG), "\\psi_1^{right} = ", 24, right);
        printer.printSciDouble((*LOG), shotParameter_2, precision, 15, left);
        (*LOG) << std::endl;
        (*LOG) << "Begin shooting:" << std::endl;
        printer.printString((*LOG), "nStep",           10, left);
        printer.printString((*LOG), "\\psi_1^{left}",  20, left);
        printer.printString((*LOG), "\\psi_1^{right}", 20, left);
        printer.printString((*LOG), "\\psi_1^{test}",  20, left);
        printer.printString((*LOG), "\\psi_0^{test}",  20, left);
        printer.printString((*LOG), "error",           20, left);
        printer.printString((*LOG), "time[ms]",        20, left);
        (*LOG) << std::endl;
    }

    error = abs(shotParameter_2 - shotParameter_1)/abs(shotParameter_1 + shotParameter_2);
    while (error > eps) {
        if (logStreamIsSet || printLogOn) { // write (*LOG)
            printer.printInt((*LOG), nStep, 10, left);
            printer.printSciDouble((*LOG), shotParameter_1, precision, 20, left);
            printer.printSciDouble((*LOG), shotParameter_2, precision, 20, left);
            printer.printSciDouble((*LOG), 0.5*(shotParameter_1 + shotParameter_2), precision, 20, left);
            timer.start();
        }

        psiCurrent[0] = phi_1;
	    psiCurrent[1] = phi_1;
        psiCurrent[2] = (shotParameter_1 + shotParameter_2)/2;
        psiCurrent[3] = psiCurrent[2];
        solver.Integrate(rhs, psiCurrent, xFrom, xTo);
        if (psiCurrent[2] > 0)
            shotParameter_2 -= (shotParameter_2 - shotParameter_1)/2;
        else shotParameter_1 += (shotParameter_2 - shotParameter_1)/2;

        if (logStreamIsSet || printLogOn) { // write (*LOG)
            timer.stop();
            timePerStep = timer.getElapsedTimeInMilliSec();
            timeSteps[nStep] = timePerStep;
            printer.printSciDouble((*LOG), psiCurrent[2], precision, 20, left);
            printer.printSciDouble((*LOG), error, precision, 20, left);
            printer.printSciDouble((*LOG), timePerStep, precision, 20, left);
            (*LOG) << std::endl;
            timer.reset();
            ++nStep;
        }
        // update error
        error = abs(shotParameter_2 - shotParameter_1)/abs(shotParameter_1 + shotParameter_2);
    }
    psi_1 = (shotParameter_1 + shotParameter_2)/2;
    if (logStreamIsSet || printLogOn) { // write last line in table, overall timing, and potential output
        timer.start();
        psiCurrent[0] = phi_1;
        psiCurrent[1] = phi_1;
        psiCurrent[2] = psi_1;
        psiCurrent[3] = psi_1;
        ODEdata outData(-1);
        solver.SetOutput(outData);
        solver.Integrate(rhs, psiCurrent, xFrom, xTo);
        timer.stop();
        timePerStep = timer.getElapsedTimeInMilliSec();
        timeSteps[nStep] = timePerStep;
        printer.printInt((*LOG), nStep, 10, left);
        printer.printSciDouble((*LOG), shotParameter_1, precision, 20, left);
        printer.printSciDouble((*LOG), shotParameter_2, precision, 20, left);
        printer.printSciDouble((*LOG), 0.5*(shotParameter_1 + shotParameter_2), precision, 20, left);
        printer.printSciDouble((*LOG), psiCurrent[2], precision, 20, left);
        printer.printSciDouble((*LOG), error, precision, 20, left);
        printer.printSciDouble((*LOG), timePerStep, precision, 20, left);
        (*LOG) << std::endl;
        (*LOG) << "Finally selected \\psi_{1} = ";
        printer.printSciDouble((*LOG), psi_1, precision, 20, left);
        (*LOG) << std::endl;
        overallTime = 0.0;
        for (Int i = 0; i <= nStep; ++i) {
            overallTime += timeSteps[i];
        }
        averageTime = overallTime/(nStep + 1);
        (*LOG) << "Overall time = ";
        printer.printSciDouble((*LOG), overallTime, precision, 15, left);
        (*LOG) << " ms" << std::endl;
        (*LOG) << "Average time per step = ";
        printer.printSciDouble((*LOG), averageTime, precision, 15, left);
        (*LOG) << " ms" << std::endl;
        (*LOG) << "output steps for potential and its correction:" << std::endl;
        Int size = outData.Count();
        printer.printString((*LOG), "x");
        printer.printString((*LOG), "phi");
        printer.printString((*LOG), "dphi");
        printer.printString((*LOG), "psi");
        printer.printString((*LOG), "dpsi");
        (*LOG) << std::endl;
        for (Int i = 0; i < size; ++i) {
            printer.printSciDouble((*LOG), outData.xSave[i], precision);
            printer.printSciDouble((*LOG), outData.ySave[0][i], precision);
            printer.printSciDouble((*LOG), outData.ySave[1][i], precision);
            printer.printSciDouble((*LOG), outData.ySave[2][i], precision);
            printer.printSciDouble((*LOG), outData.ySave[3][i], precision);
            (*LOG) << std::endl;
        }
        solver.SetOutput(calculatedData);
    }
    calculated = true;
}

Double FTTFQEpotential::operator() (Double x) {
    if (!calculated) calculate();
    if (x < 0.0) { 
        std::cerr << "cannot calculate potential at x < 0" << std::endl; 
        exit(0);
    }
    if (x < 1.0) {
        DoubleVec psiCurrent(0.0, rhsFTTFQEpotential::dim);
        Double xFrom = 1.0;
        Double xTo = x;
        psiCurrent[0] = phi(1);
        psiCurrent[1] = phi(1);
        psiCurrent[2] = psi_1;
        psiCurrent[3] = psi_1;
        solver.Integrate(rhs, psiCurrent, xFrom, xTo);
        return psiCurrent[2];
    }
    else return psi_1;
}

Double FTTFQEpotential::derivative(const Double x) {
    if (!calculated) calculate();
    if (x < 0.0) { 
        std::cerr << "cannot calculate potential at x < 0" << std::endl; 
        exit(0);
    }
    if (x < 1.0) {
        DoubleVec psiCurrent(0.0, rhsFTTFQEpotential::dim);
        Double xFrom = 1.0;
        Double xTo = x;
        psiCurrent[0] = phi(1);
        psiCurrent[1] = phi(1);
        psiCurrent[2] = psi_1;
        psiCurrent[3] = psi_1;
        solver.Integrate(rhs, psiCurrent, xFrom, xTo);
        return psiCurrent[3];
    }
    else return psi_1;
}

void FTTFQEpotential::setParameters(const PhysQ& V, const PhysQ& T, const Double Z) {
    PhysQ V1;
    PhysQ T1;  
    V1.setValue(V()*Z);
    T1.setValue(T()*pow(Z, -4.0/3.0));
	if (abs(log10(V1()) - log10(Vol)) > 1e-10
        || abs(log10(T1()) - log10(Temp)) > 1e-10) 
    { 
        calculated = false;
        Vol = V1();
        Temp = T1();

        if (printLogOn) {
            std::stringstream filename;
            filename << "log/log_psi(V=" << Vol 
                     << ", T=" << Temp << ")_";
            filename << timer.getCurrentDatetime();
            filename << ".txt";
            if (LOG->is_open()) LOG->close();
            LOG->open(filename.str().c_str(), std::ios::out);
        }

        Double a =   pow(2.0, 7.0/6.0)
                   * pow(3.0, 2.0/3.0)
                   * pow(M_PI, -5.0/3.0)
                   * sqrt(Temp)*pow(Vol, 2.0/3.0);
        rhs.updateParameter(a);
        phi.setParameters(V, T, Z);
    }
}

void FTTFQEpotential::setTolerance(Double _eps) {
    eps = _eps;
	phi.setTolerance(eps);
    solver.SetTolerance(0.0, eps/10);
    calculated = false;
}

void FTTFQEpotential::printData(const char* filename) {
    OUT = new std::ofstream;
    OUT->open(filename, std::ios::out);
    performOutput(-1);
    OUT->close();
}

void FTTFQEpotential::printData(const char* filename, const Int Npoints) {
    OUT = new std::ofstream;
    OUT->open(filename, std::ios::out);
    performOutput(Npoints);
    OUT->close();
}

DoubleMat FTTFQEpotential::getData() {
    performOutput(-1);
    return data;
}

DoubleMat FTTFQEpotential::getData(const Int Npoints) {
    performOutput(Npoints);
    return data;
}   

void FTTFQEpotential::performOutput(Int Npoints) {
    if (!calculated) calculate();
    bool outIsSet = (OUT != NULL);
    // additional integration with saving steps
    Double xFrom = 1.0;
    Double xTo = 0.0;
    DoubleVec psiCurrent(0.0, rhsFTTFQEpotential::dim);
    psiCurrent[0] = phi(1);
    psiCurrent[1] = phi(1);
    psiCurrent[2] = psi_1;
    psiCurrent[3] = psi_1;
    ODEdata outData(-1);
    solver.SetOutput(outData);
    solver.Integrate(rhs, psiCurrent, xFrom, xTo);
    // print output data
    Int size = outData.Count();
    data = DoubleMat(0.0, rhsFTTFQEpotential::dim + 1, size);
    if (outIsSet) {
        printer.printString((*OUT), "x");
        printer.printString((*OUT), "phi");
        printer.printString((*OUT), "dphi");
        printer.printString((*OUT), "psi");
        printer.printString((*OUT), "dpsi");
        (*OUT) << std::endl;
        for (Int i = 0; i < size; ++i) {
            printer.printSciDouble((*OUT), data[0][i] = outData.xSave[i],    precision);
            printer.printSciDouble((*OUT), data[1][i] = outData.ySave[0][i], precision);
            printer.printSciDouble((*OUT), data[2][i] = outData.ySave[1][i], precision);
            printer.printSciDouble((*OUT), data[3][i] = outData.ySave[2][i], precision);
            printer.printSciDouble((*OUT), data[4][i] = outData.ySave[3][i], precision);
            (*OUT) << std::endl;
        }
    }
    else {
        for (Int i = 0; i < size; ++i) {
            data[0][i] = outData.xSave[i];
            data[1][i] = outData.ySave[0][i];
            data[2][i] = outData.ySave[1][i];
            data[3][i] = outData.ySave[2][i];
            data[4][i] = outData.ySave[3][i];
        }
    }
    solver.SetOutput(calculatedData);
}

void FTTFQEpotential::setLogStream(std::ofstream* _LOG) {
    if (printLogOn) setPrintLogOff();
    LOG = _LOG;
    phi.setLogStream(LOG);
    logStreamIsSet = true;
}

void FTTFQEpotential::clearLogStream() {
    if (logStreamIsSet) {
        LOG = NULL;
        logStreamIsSet = false;
        phi.clearLogStream();
    }
}

void FTTFQEpotential::setPrintLogOn() {
    if (!logStreamIsSet) {
        LOG = new std::ofstream;
        phi.setLogStream(LOG);
        printLogOn = true;
    }
}

void FTTFQEpotential::setPrintLogOff() {
    if (!logStreamIsSet) {
        if (printLogOn) {
            printLogOn = false;
            LOG->close();
            delete LOG;
            LOG = NULL;
        }
    }
}