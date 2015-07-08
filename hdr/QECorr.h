#pragma once
#include "../lib/numeric/hdr/ODE/ODEsolver.h"
#include "../lib/numeric/hdr/ODE/ODEdata.h"
#include "../lib/numeric/hdr/ODE/steppers/ODEstepperPD853.h"
#include "../lib/numeric/hdr/numericTypes.h"
#include "../lib/numeric/hdr/mathUtils.h"
#include "../lib/numeric/hdr/SF/FermiDirac.h"
#include "../hdr/ThermodynamicFunction.h"
#include "../hdr/FTTFpotential.h"
#include "../hdr/FTTFQEpotential.h"
#include "../hdr/Units.h"
#include "../hdr/Timer.h"
#include "../hdr/Printer.h"
#include "../hdr/stringUtils.h"
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>

#define DEFAULT_INPUT_FILE "in/FTTFQEinput.dat"
#define DEFAULT_OUTPUT_FILE "out/FTTFQEoutput.dat"
/**
* @brief This class implements interface for Thomas-Fermi model.
* @details The main formulas for calculation of thermodynamic quantities
* are presented below and expressed throgh the FTTF potential:
* - the formula for the pressure calculation: 
* @f[
*	P = \frac{(2T)^{5/2}}{6\pi^2}I_{3/2}(\phi(1)).
* @f]
* - the formula for the chemical potential calculation:
* @f[
*	\mu = T\phi(1).
* @f]
* - The methods for calculating energy and entropy are presented inside model class
* - Thermal parts are caculated in the following way:
*	@f[
*	   P_T = P - \left.P\right|_{T = 0}.
*   @f]
*   @f[
*	   E_T = E - \left.E\right|_{T = 0}.
*   @f]
*   @f[
*	   S_T = S - \left.S\right|_{T = 0}.
*   @f]
*   @f[
*	   \mu_T = \mu - \left.\mu\right|_{T = 0}.
*	@f]
*/
class QECorr {
public:
	/**
	* @brief A constructor of Thomas-Fermi model.
	*/
    QECorr(Double _Z = 1.0, Double _Mass = 1.0);
    ~QECorr();
    /**
	* @brief Set tolerance eps for the further calculations. Default is @f$ 10^{-6} @f$.
	*/
    void setTolerance(const Double eps);
    /**
	* @brief Set temperature and volume/density/concentration range for calculation.
	* @details Examples:
	* volRange: V[Unit,Scaling](Vstart, Vend, Npoints), here C, D can be used instead of V
	* tempRange: T[Unit,Scaling](Tstart, Tend, Npoints)
	*/
    void setParameters(std::string volRange, std::string tempRange);
    /**
	* @brief Calculate output data as defined in inputString.
	* @details Example of string: P[GPa,log] PT[GPa] E[eV,lin] S M[eV] ...
	* List of possible parameters (units and scaling are optional, by default
	* atomic units and linear scaling are used):
	* - P[Units, Scaling] - full FTTF pressure, units: Pa, MBar, GPa, Atomic;
	* - PT[Units, Scaling] - thermal FTTF pressure;
	* - PC[Units, Scaling] - cold FTTF pressure;
	* - E[Units, Scaling] - full FTTF energy, units: eV, Atomic;
	* - ET[Units, Scaling] - thermal FTTF energy;
	* - EC[Units, Scaling] - cold FTTF energy;
	* - S - full FTTF entropy;
 	* - ST - thermal FTTF entropy;
 	* - SC - cold FTTF entropy;
 	* - M[Units, Scaling] - full FTTF chemical potential, units: eV, Atomic;
	* - MT[Units, Scaling] - thermal FTTF chemcial potential;
	* - MC[Units, Scaling] - cold FTTF chemical potential;
	*/
    void calculateData(std::string inputString);
    /**
	* @brief Print output data in file.
	*/
    void printOutput(const char* filename);
    /**
	* @brief Print output data in file.
	*/
    void printOutput(std::string filename);
    /**
    * @brief Write main log to specified stream.
    */
    void setMainLogStream(std::ofstream* _mainLOG);
    /**
    * @brief Write log for particular point to specified stream.
    */
    void setPointLogStream(std::ofstream* _pointLOG);
    /**
    * @brief Enable self-printing main log into file.
    */
    void setPrintMainLogOn();
    /**
    * @brief Disable self-printing main log into file.
    */
    void setPrintMainLogOff();
    /**
    * @brief Enable self-printing log for particular point into file.
    */
    void setPrintPointLogOn();
    /**
    * @brief Disable self-printing log for particular point into file.
    */
    void setPrintPointLogOff();

    void setShowProgressOn();
    void setShowProgressOff();
	    
    PhysQvec& operator[] (std::string quantity);
    /**
	* @brief Problem for calculating energy.
	* @details The calculation goes through solving the ODE with parameters
	* @f$ p_0 = 2^{7/6}3^{2/3}\pi^{-5/3}T^{1/2}V^{2/3} @f$,
	* @f$ p_1 = \pi^{-3}VT^{2} @f$,
	* @f$ \Delta E_0 = -0.269900170 @f$. Correction to Thomas-Fermi energy 
	* is solution of the ODE:
	* @f{eqnarray*}{
	*  \frac{d^2\phi}{dx^2} &=& p_0 x I_{1/2}
	*	  \left(
	*		\frac{\phi(x)}{x}
	*	  \right), \\
	*  \frac{d^2\psi}{dx^2} &=& p_0 \left[
	*		       \frac12 \psi(x) I_{-1/2}\left(\frac{\phi(x)}{x}\right) + xY'(x)
	*           \right], \\
	*  \Delta E'(x) &=& -p_1 x \left[
	*						\frac12 \psi(x) I_{1/2}\left(\frac{\phi(x)}{x}\right) + xY(x)
	*					\right], \\
	* \phi(1) &=& \phi'(1), \\
	* \psi(1) &=& \psi'(1), \\
	* \Delta E(1) &=& \frac{1}{3\pi}\sqrt{\frac{T}{2}}\psi'(0) - \Delta E_0,
	* @f}
	* Finally, @f$ \Delta E = \Delta E(0) @f$.
	*/
    struct rhsQECorrEnergy {
		Double a;
		Double b;
		static const Int dim = 4;
		FermiDirac<Mhalf> FDmhalf;
		FermiDirac<Half> FDhalf;
		Yfunction Y;
		rhsQECorrEnergy(Double _a = 0, Double _b = 0) : a(_a), b(_b) {}
		void updateParameters(const Double _a, const Double _b) { a = _a; b = _b; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
    		static Int sgn;
    		dydx[0] = y[1];
    		dydx[2] = y[3];
    		if (x > 0) {
    		    dydx[1] = a*x*FDhalf(y[0]/x);
    		    dydx[3] = a*(FDmhalf(y[0]/x)/2.0*y[2] + x*Y.derivative(y[0]/x));
    		    dydx[4] = -b*x*(0.5*y[2]*FDhalf(y[0]/x) + x*Y(y[0]/x));
    		    dydx[4] > 0 ? sgn = 1 : sgn = -1;
    		}
    		else {
    		    dydx[1] = dydx[3] = 1e+10;
    		    dydx[4] = sgn*1e+10;
    		}
		}
	};
	/**
	* @brief Problem for calculating entropy.
	* @details The calculation goes through solving the ODE with parameters
	* @f$ p_0 = 2^{7/6}3^{2/3}\pi^{-5/3}T^{1/2}V^{2/3} @f$,
	* @f$ p_1 = \pi^{-3}VT @f$.
	* Correction to Thomas-Fermi entropy is solution of the ODE:
	* @f{eqnarray*}{
	*  \frac{d^2\phi}{dx^2} &=& p_0 x I_{1/2}
	*	  \left(
	*		\frac{\phi(x)}{x}
	*	  \right), \\
	*  \frac{d^2\psi}{dx^2} &=& p_0 \left[
	*		       \frac12 \psi(x) I_{-1/2}\left(\frac{\phi(x)}{x}\right) + xY'(x)
	*           \right], \\
	*  \Delta S'(x) &=& -p_1 x \left[
	*						\frac12 \psi(x) I_{1/2}\left(\frac{\phi(x)}{x}\right) + 2xY(x)
	*					\right], \\
	* \phi(1) &=& \phi'(1), \\
	* \psi(1) &=& \psi'(1), \\
	* \Delta S(1) &=& \frac{1}{3\pi\sqrt{2T}}\psi'(0),
	* @f}
	* Finally, @f$ \Delta S = \Delta S(0) @f$.
	*/
	struct rhsQECorrEntropy {
		Double a;
		Double b;
		static const Int dim = 4;
		FermiDirac<Half> FDhalf;
		FermiDirac<Mhalf> FDmhalf;
		Yfunction Y;
		rhsQECorrEntropy(Double _a = 0, Double _b = 0) : a(_a), b(_b) {}
		void updateParameters(const Double _a, const Double _b) { a = _a; b = _b; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
			static Int sgn;
    		dydx[0] = y[1];
    		dydx[2] = y[3];
    		if (x > 0) {
    		    dydx[1] = a*x*FDhalf(y[0]/x);
    		    dydx[3] = a*(FDmhalf(y[0]/x)/2.0*y[2] + x*Y.derivative(y[0]/x));
    		    dydx[4] = -b*x*(0.5*y[2]*FDhalf(y[0]/x) + 2*x*Y(y[0]/x));
    		    dydx[4] > 0 ? sgn = 1 : sgn = -1;
    		}
    		else {
    		    dydx[1] = dydx[3] = 1e+10;
    		    dydx[4] = sgn*1e+10;
    		}
		}
	};

private:
	Double Z;
	Double Mass;
	Double eps;
	std::map<std::string, Unit> unit;
	std::map<std::string, Scaling> scale;
	std::map<std::string, bool> need;
	std::map<std::string, bool> calculated;
	std::map<std::string, PhysQvec*> data;
	static const Int varNum = 16;
	std::string varName[varNum];
	Int inputVars;
	Int coldVarBound;
	Int inputVarBound;
	Int Vsize;
	Int Tsize;
	// calculate single point
	void calculate(std::string var, Int v, Int t);
	void  calculateDP(Int v, Int t); 
	void calculateDPT(Int v, Int t); 
	void calculateDPC(Int v); 
	void  calculateDE(Int v, Int t); 
	void calculateDET(Int v, Int t); 
	void calculateDEC(Int v); 
	void  calculateDS(Int v, Int t); 
	void calculateDST(Int v, Int t); 
	void calculateDSC(Int v); 
	void  calculateDM(Int v, Int t); 
	void calculateDMT(Int v, Int t); 
	void calculateDMC(Int v); 

	void calculateAll(); // uses fast calculation when prepare output
	void performOutput();
	void clearData();
	
    // ode solvers
    rhsQECorrEnergy rhsEnergy;
    ODEsolver<ODEstepperPD853<rhsQECorrEnergy> > energySolver;
    ODEdata energyData;

    rhsQECorrEntropy rhsEntropy;
    ODEsolver<ODEstepperPD853<rhsQECorrEntropy> > entropySolver;
    ODEdata entropyData;

    rhsQECorrEnergy rhsColdEnergy;
    ODEsolver<ODEstepperPD853<rhsQECorrEnergy> > coldEnergySolver;
    ODEdata coldEnergyData;

    rhsQECorrEntropy rhsColdEntropy;
    ODEsolver<ODEstepperPD853<rhsQECorrEntropy> > coldEntropySolver;
    ODEdata coldEntropyData;

    //transform quantities
    void transformDtoV();
    void transformCtoV();
    void transformVtoD();
    void transformVtoC();

    FermiDirac<ThreeHalf> FD3half;
    FermiDirac<Half> FDhalf;
	FermiDirac<Mhalf> FDmhalf;
	Yfunction Y;
    PhysQ coldT;

    FTTFpotential phi;
    FTTFpotential coldPhi;

    FTTFQEpotential psi;
    FTTFQEpotential coldPsi;

    // log and output
    Printer printer;
    Timer pointTimer;
    Timer mainTimer;
    bool printMainLogOn;
    bool printPointLogOn;
    bool mainLogStreamIsSet;
    bool pointLogStreamIsSet;
    std::ofstream* pointLOG;
    std::ofstream* mainLOG;
    std::ofstream* OUT;
    Int precision;

    bool isLogVol;
	bool isLogMdns;
	bool isLogVdns;

	bool showProgress;
};

QECorr::QECorr(Double _Z, Double _Mass) : 
	Z(_Z), Mass(_Mass), eps(1e-6),
	energyData(-1), // save all steps
	energySolver(1e-6, 0.0),
	rhsEnergy(0.0),
	coldEnergyData(-1), // save all steps
	coldEnergySolver(1e-6, 0.0),
	rhsColdEnergy(0.0),
	entropyData(-1), // save all steps
	entropySolver(1e-6, 0.0),
	rhsEntropy(0.0),
	coldEntropyData(-1), // save all steps
	coldEntropySolver(1e-6, 0.0),
	rhsColdEntropy(0.0),
	coldT(1e-6) {

	varName[0] = "T"; varName[1] = "V"; varName[2] = "D"; varName[3] = "C";  
	inputVarBound = 4; 
	varName[4] = "DPC"; varName[5] = "DEC"; varName[6] = "DSC"; varName[7] = "DMC";
	coldVarBound = 8;
	varName[8] = "DP"; varName[9] = "DE"; varName[10] = "DS"; varName[11] = "DM";
	varName[12] = "DPT"; varName[13] = "DET"; varName[14] = "DST"; varName[15] = "DMT";
	
	for (Int i = 0; i < varNum; ++i) {
		data[varName[i]] = NULL;
		unit[varName[i]] = Atomic;
		scale[varName[i]] = lin;
		need[varName[i]] = false;
		if (i >= inputVarBound) calculated[varName[i]] = false;
	}
	precision = static_cast<int>(-log10(eps));
	energySolver.SetOutput(energyData);
	coldEnergySolver.SetOutput(coldEnergyData);
	entropySolver.SetOutput(entropyData);
	coldEntropySolver.SetOutput(coldEntropyData);
	pointLOG = NULL;
	mainLOG = NULL;
	OUT = NULL;
    pointLogStreamIsSet = false;
    printPointLogOn = false;
    mainLogStreamIsSet = false;
    printMainLogOn = false;
    isLogVol = false;
    isLogMdns = false;
    isLogVdns = false;
    showProgress = false;
}

QECorr::~QECorr() {
	clearData();
	setPrintMainLogOff();
	setPrintPointLogOff();
}

void QECorr::calculate(std::string quantity, Int v = -1, Int t = -1) {
	if (v < 0 && t < 0) {
		calculateData(quantity);
	}
	if (v >= 0 && t < 0) {
			 if (!quantity.compare("DPC")) calculateDPC(v);
		else if (!quantity.compare("DEC")) calculateDEC(v);
		else if (!quantity.compare("DSC")) calculateDSC(v);
		else if (!quantity.compare("DMC")) calculateDMC(v);
	}
	if (v >= 0 && t >= 0) {
			 if (!quantity.compare("DP"))  calculateDP(v, t);
		else if (!quantity.compare("DE"))  calculateDE(v, t);
		else if (!quantity.compare("DS"))  calculateDS(v, t);
		else if (!quantity.compare("DM"))  calculateDM(v, t);
		else if (!quantity.compare("DPT")) calculateDPT(v, t);
		else if (!quantity.compare("DET")) calculateDET(v, t);
		else if (!quantity.compare("DST")) calculateDST(v, t);
		else if (!quantity.compare("DMT")) calculateDMT(v, t);
	}
}

void QECorr::setTolerance(const Double _eps) { 
	eps = _eps;
	precision = static_cast<int>(-log10(eps));
	if (printMainLogOn || mainLogStreamIsSet) {
		*mainLOG << "QECorr accepted new tolerance, eps = ";
		printer.printSciDouble(*mainLOG, eps, precision);
		*mainLOG << std::endl;
		psi.setLogStream(mainLOG);
		coldPsi.setLogStream(mainLOG);
	}
	phi.setTolerance(eps);
	psi.setTolerance(eps);
	if (printMainLogOn || mainLogStreamIsSet) 
		*mainLOG << "Cold: " << std::endl;
	coldPhi.setTolerance(eps);
	coldPsi.setTolerance(eps);
	energySolver.SetTolerance(0.0, eps/10);
	coldEnergySolver.SetTolerance(0.0, eps/10);
	entropySolver.SetTolerance(0.0, eps/10);
	coldEntropySolver.SetTolerance(0.0, eps/10);
	if (printMainLogOn || mainLogStreamIsSet) {
		psi.clearLogStream();
		coldPsi.clearLogStream();	
	}
}

void QECorr::setMainLogStream(std::ofstream* _mainLOG) {
    if (printMainLogOn) setPrintMainLogOff();
    mainLOG = _mainLOG;
    mainTimer.start();
    mainLogStreamIsSet = true;
}

void QECorr::setPrintMainLogOn() {
    if (!mainLogStreamIsSet) {
        mainLOG = new std::ofstream;
        mainTimer.start();
        std::stringstream filename;
        filename << "log/log_QECorr_Z(";
        filename << Z << ")_M("; 
        filename << Mass << ")_";
        filename << mainTimer.getCurrentDatetime();
        filename << ".txt";
        mainLOG->open(filename.str().c_str(), std::ios::out);
        *mainLOG << "QECorr log is started" << std::endl;
        *mainLOG << "Z = " << Z << std::endl;
        *mainLOG << "Atomic Mass = " << Mass << std::endl;
        printMainLogOn = true;
    }
}

void QECorr::setPrintMainLogOff() {
    if (!mainLogStreamIsSet) {
        if (printMainLogOn) {
            printMainLogOn = false;
            if (mainLOG->is_open()) mainLOG->close();
            delete mainLOG;
            mainLOG = NULL;
        }
    }
}

void QECorr::setPointLogStream(std::ofstream* _pointLOG) {
    if (printPointLogOn) setPrintPointLogOff();
    pointLOG = _pointLOG;
    pointLogStreamIsSet = true;
}

void QECorr::setPrintPointLogOn() {
    if (!pointLogStreamIsSet) {
        pointLOG = new std::ofstream;
        printPointLogOn = true;
    }
}

void QECorr::setPrintPointLogOff() {
    if (!pointLogStreamIsSet) {
        if (printPointLogOn) {
            printPointLogOn = false;
            if (pointLOG->is_open()) pointLOG->close();
            delete pointLOG;
            pointLOG = NULL;
        }
    }
}

void QECorr::setShowProgressOn() {
	showProgress = true;
}

void QECorr::setShowProgressOff() {
	showProgress = false;
}

void QECorr::clearData() {
	for (Int i = 0; i < varNum; ++i) {
		if (data[varName[i]] != NULL) delete data[varName[i]];
		if (i >= inputVarBound) { calculated[varName[i]] = false; }
	}
}

void QECorr::setParameters(std::string volRange, std::string tempRange) {
	Double localTime = 0;
	if (printMainLogOn || mainLogStreamIsSet) { 
		localTime = mainTimer.getElapsedTimeInMilliSec();
	}
	clearData();
	const char startUnit  = '[';
    const char endUnit    = ']';
	const char startRange = '(';
	const char endRange   = ')';
    const char spaceDelim = ' ';
    const char rangeDelim = ',';
    std::string volUnitScale; 
    std::string tempUnitScale;
    // main parameters
    bool isVolume = (volRange[0] == 'V');
    bool isMdns = (volRange[0] == 'D');
    bool isVdns = (volRange[0] == 'C');
    Unit volUnit = Atomic;
    Scaling volScale = lin;
    Unit tempUnit = Atomic;
    Scaling tempScale = lin;
	Double  volRangeStart = 1.0; Double  volRangeEnd = 1.0; Vsize = 1; 
	Double tempRangeStart = 1.0; Double tempRangeEnd = 1.0; Tsize = 1;
    // analyze volRange and tempRange whether there are [] and ()
    Int volUnitStartPos   = volRange.find(startUnit);
    Int volUnitEndPos     = volRange.find(endUnit);
    Int volRangeStartPos  = volRange.find(startRange);
    Int volRangeEndPos    = volRange.find(endRange);
    Int tempUnitStartPos  = tempRange.find(startUnit);
    Int tempUnitEndPos    = tempRange.find(endUnit);
	Int tempRangeStartPos = tempRange.find(startRange);
	Int tempRangeEndPos   = tempRange.find(endRange);
	Int length;

	// parse volUnit and tempUnit:
	std::vector<std::string> tokens;
	if (volUnitStartPos != std::string::npos && 
		volUnitEndPos   != std::string::npos) {
		length = volUnitEndPos - volUnitStartPos - 1;
		volUnitScale = std::string(volRange, volUnitStartPos + 1, length);
		tokens = split(volUnitScale, rangeDelim);
		if (tokens.size() > 0) {
			volUnit = stringToUnit(tokens[0]);
		}
		if (tokens.size() > 1) {
			volScale = stringToScale(tokens[1]);
		}
	}
	
	if (tempUnitStartPos != std::string::npos && 
		tempUnitEndPos   != std::string::npos) {
		length = tempUnitEndPos - tempUnitStartPos - 1;
		tempUnitScale = std::string(tempRange, tempUnitStartPos + 1, length);
		tokens = split(tempUnitScale, rangeDelim);
		if (tokens.size() > 0) {
			tempUnit = stringToUnit(tokens[0]);
		}
		if (tokens.size() > 1) {
			tempScale = stringToScale(tokens[1]);
		}	
	}
	
    if (volRangeStartPos != std::string::npos && 
		volRangeEndPos   != std::string::npos) {
		length = volRangeEndPos - volRangeStartPos - 1;
		volRange = std::string(volRange, volRangeStartPos + 1, length);
		tokens = split(volRange, rangeDelim);
		if (tokens.size() > 0) {
			volRangeStart = atof(tokens[0].c_str());
		}
		if (tokens.size() > 1) {
			volRangeEnd = atof(tokens[1].c_str());
		}
		if (tokens.size() > 2) {
			Vsize = atoi(tokens[2].c_str());	
		}
	}
	
	if (tempRangeStartPos != std::string::npos && 
		tempRangeEndPos   != std::string::npos) {
		length = tempRangeEndPos - tempRangeStartPos - 1;
		tempRange = std::string(tempRange, tempRangeStartPos + 1, length);
		tokens = split(tempRange, rangeDelim);
		if (tokens.size() > 0) {
			tempRangeStart = atof(tokens[0].c_str());
		}
		if (tokens.size() > 1) {
			tempRangeEnd = atof(tokens[1].c_str());
		}
		if (tokens.size() > 2) {
			Tsize = atoi(tokens[2].c_str());	
		}
	}
	/****************************LOG Section***************************/
	if (printMainLogOn || mainLogStreamIsSet) {                       //
		*mainLOG << "QECorr accepts new parameters:" << std::endl;//
		if (isVolume) *mainLOG << "Volume: V[";                       //
		if (isMdns) *mainLOG << "Mass density: D[";                   //
		if (isVdns) *mainLOG << "Number density: C[";                 //
		*mainLOG << unitToString(volUnit) << ","                      //
			<< scaleToString(volScale);                               //
		*mainLOG << "](";                                             //
		*mainLOG << volRangeStart << ",";                             //
		*mainLOG << volRangeEnd   << ",";                             //
		*mainLOG << Vsize         << ")";                             //
		*mainLOG << std::endl;                                        //
		*mainLOG << "Temperature: T[";                                //
		*mainLOG << unitToString(tempUnit) << ","                     //
			<< scaleToString(tempScale);                              //
		*mainLOG << "](";                                             //
		*mainLOG << tempRangeStart << ",";                            //
		*mainLOG << tempRangeEnd   << ",";                            //
		*mainLOG << Tsize          << ")";                            //
		*mainLOG << std::endl;                                        //
		*mainLOG << "elapsed time: " <<                               //
			mainTimer.getElapsedTimeInMilliSec() - localTime;         //
		*mainLOG << " ms" << std::endl;                               //
		*mainLOG << "Allocating memory for input parameters, ";       //
		localTime = mainTimer.getElapsedTimeInMilliSec();             //
	}                                                                 //
	/******************************************************************/
	// generate arrays
	PhysQvec* Volume = new PhysQvec(Vsize);
	PhysQvec* Temperature = new PhysQvec(Tsize);
	if (Vsize > 1) {
		Double V = volRangeStart;
		Double dV = (volRangeEnd - volRangeStart)/(Vsize - 1);
		for (Int v = 0; v < Vsize; ++v) {
			(*Volume)[v].setValue(V, volUnit, volScale);
			V += dV;
		}
	}
	else {
		(*Volume)[0].setValue(volRangeStart, volUnit, volScale);
	}

	if (Tsize > 1) {
		Double T = tempRangeStart;
		Double dT = (tempRangeEnd - tempRangeStart)/(Tsize - 1);
		for (Int t = 0; t < Tsize; ++t) {
			(*Temperature)[t].setValue(T, tempUnit, tempScale);
			T += dT;
		}
	}
	else {
		(*Temperature)[0].setValue(tempRangeStart, tempUnit, tempScale);
	}
	//assign arrays:
	data["T"] = Temperature;
	if (isVolume) {
		data["V"] = Volume; 
		transformVtoC();
		transformVtoD();
	}
	if (isMdns) {
		data["D"] = Volume;
		transformDtoV();
		transformVtoC();
	}
	if (isVdns) {
		data["C"] = Volume;
		transformCtoV();
		transformVtoD();
	}
	if (printMainLogOn || mainLogStreamIsSet) {
		*mainLOG << "elapsed time: " << mainTimer.getElapsedTimeInMilliSec() - localTime;
		*mainLOG << " ms" << std::endl;
		isLogVdns = isVdns;
		isLogMdns = isMdns;
		isLogVol = isVolume;
	}
}

void QECorr::transformDtoV() {
	PhysQvec* D = data["D"];
	PhysQvec* V = new PhysQvec(Vsize);
	Double currentV;
	for (Int i = 0; i < Vsize; ++i) {
		currentV = Mass/(Avogadro*(aVol*1e+6))/(*D)[i](gOverCmc);
		(*V)[i].setValue(currentV);
	}
	data["V"] = V;
}

void QECorr::transformCtoV() {
	PhysQvec* C = data["C"];
	PhysQvec* V = new PhysQvec(Vsize);
	Double currentV;
	for (Int i = 0; i < Vsize; ++i) {
		currentV = 1.0/(*C)[i]();
		(*V)[i].setValue(currentV);
	}
	data["V"] = V;
}

void QECorr::transformVtoD() {
	PhysQvec* V = data["V"];
	PhysQvec* D = new PhysQvec(Vsize);
	Double currentD;
	for (Int i = 0; i < Vsize; ++i) {
		currentD = Mass/(Avogadro*(aVol*1e+6))/(*V)[i]();
		(*D)[i].setValue(currentD, gOverCmc);
	}
	data["D"] = D;
}

void QECorr::transformVtoC() {
	PhysQvec* V = data["V"];
	PhysQvec* C = new PhysQvec(Vsize);
	Double currentC;
	for (Int i = 0; i < Vsize; ++i) {
		currentC = 1.0/(*V)[i]();
		(*C)[i].setValue(currentC);
	}
	data["C"] = C;
}

void QECorr::calculateData(std::string inputString) {
	Double localTime;
	/****************************LOG Section*******************/
	if (printMainLogOn || mainLogStreamIsSet) {               //
		localTime = mainTimer.getElapsedTimeInMilliSec();     //
		*mainLOG << "QECorr accepts physical";            //
		*mainLOG << " quantities to calculate:" << std::endl; //
	}                                                         //
	/**********************************************************/
	const char startUnit = '[';
    const char endUnit = ']';
    const char delim = ' ';
    const char unitScaleDelim = ',';
    // split inputString into tokens
    std::vector<std::string> tokens = split(inputString, delim);
    // analyze tokens
    std::vector<std::string>::iterator itoken;
    for (itoken = tokens.begin(); itoken != tokens.end(); ++itoken) {
    	std::string newQuantity, newUnit, newScale;
    	newQuantity = *itoken;
    	if (!isalpha(*(newQuantity.end() - 1))) {
			newUnit = split(newQuantity, startUnit)[1];
			newUnit = std::string(newUnit.begin(), newUnit.end() - 1);
			newQuantity = split(newQuantity, startUnit)[0];
			if (split(newUnit, unitScaleDelim).size() > 1){
				newScale = split(newUnit, unitScaleDelim)[1];
				newUnit  = split(newUnit, unitScaleDelim)[0];
			}
    	}
		need[newQuantity] = true;
		if (!newUnit.empty())  unit[newQuantity]  = stringToUnit(newUnit);
		if (!newScale.empty()) scale[newQuantity] = stringToScale(newScale);
		/******************LOG Section*****************/
		if (printMainLogOn || mainLogStreamIsSet) {   //
			*mainLOG << newQuantity;                  //
			if (!newUnit.empty())  *mainLOG << "["    //
				<< unitToString(unit[newQuantity]);   //
			if (!newScale.empty()) *mainLOG << ","    //
				<< scaleToString(scale[newQuantity]); //
			if (!newUnit.empty())  *mainLOG << "]";   //
			*mainLOG << std::endl;                    //
		}                                             //
		/**********************************************/
    }
    /****************************LOG Section****************************/
    if (printMainLogOn || mainLogStreamIsSet) {                        //
    	*mainLOG << "elapsed time: " <<                                //
    	mainTimer.getElapsedTimeInMilliSec() - localTime << std::endl; //
    	localTime = mainTimer.getElapsedTimeInMilliSec();              //
    	*mainLOG << "Allocate memory for output, ";                    //
    }                                                                  //
    /*******************************************************************/
    // generate arrays
    for (Int i = inputVarBound; i < coldVarBound; ++i) {
    	data[varName[i]] = new PhysQvec(Vsize);
    }
    for (Int i = coldVarBound; i < varNum; ++i) {
		data[varName[i]] = new PhysQvec(Vsize*Tsize);
    }
    /****************************LOG Section****************************/
    if (printMainLogOn || mainLogStreamIsSet) {                        //
    	*mainLOG << "elapsed time: " <<                                //
    	mainTimer.getElapsedTimeInMilliSec() - localTime << std::endl; //
    	*mainLOG << std::endl;                                         //
    }                                                                  //
    /*******************************************************************/
    calculateAll();
}

void QECorr::performOutput() {
	Double dataToPrint;
	for (Int i = 0; i < varNum; ++i) {
		std::string s;
		if (need[varName[i]]) {
			s = varName[i];
			s += "[";
			s += unitToString(unit[varName[i]]);
			s += ",";
			s += scaleToString(scale[varName[i]]);
			s += "]";
			printer.printString((*OUT), s);
		}
	}
	(*OUT) << std::endl;
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			for (Int i = 0; i < varNum; ++i) {
				if (need[varName[i]]) {
					if (data[varName[i]]->size() == Vsize)
						dataToPrint = (*data[varName[i]])[v](unit[varName[i]], scale[varName[i]]);
					if (data[varName[i]]->size() == Tsize)
						dataToPrint = (*data[varName[i]])[t](unit[varName[i]], scale[varName[i]]);
					if (data[varName[i]]->size() == Vsize*Tsize)
						dataToPrint = (*data[varName[i]])[v*Tsize + t](unit[varName[i]], scale[varName[i]]);
					printer.printSciDouble((*OUT), dataToPrint, precision);
				}
			}
			(*OUT) << std::endl;
		}
	}
}

void QECorr::printOutput(const char* filename) {
	OUT = new std::ofstream;
    OUT->open(filename, std::ios::out);
	performOutput();
    OUT->close();
    delete OUT;
}

void QECorr::printOutput(std::string filename) {
	printOutput(filename.c_str());
}

void QECorr::calculateAll() {
	Double localTime = 0;
	Double pointTime = 0;
	Int progress = 0;
	// if (showProgress) {
	// 	std::cout << "progress:  0\%";
	// }
	/****************************LOG Section**************************/
	if (printMainLogOn || mainLogStreamIsSet) {                      //
		*mainLOG << "Start calculation:" << std::endl;               //
		// print table head                                          //
		for (Int i = 0; i < varNum; ++i) {                           //
			if (need[varName[i]]) {                                  //
				printer.printString(                                 //
					*mainLOG,                                        //
					varName[i] + "[" +                               //
					unitToString(unit[varName[i]]) + "," +           //
					scaleToString(scale[varName[i]]) + "]",          //
					20, left);                                       //
			}                                                        //
		}                                                            //
		printer.printString(*mainLOG, "phi(1)", 20, left);           //
		printer.printString(*mainLOG, "phi'(0)", 20, left);          //
		printer.printString(*mainLOG, "phi_C(1)", 20, left);       //
		printer.printString(*mainLOG, "phi'_C(0)", 20, left);      //
		printer.printString(*mainLOG, "psi(1)", 20, left);           //
		printer.printString(*mainLOG, "psi'(0)", 20, left);          //
		printer.printString(*mainLOG, "psi_C(1)", 20, left);       //
		printer.printString(*mainLOG, "psi'_C(0)", 20, left);      //
		printer.printString(*mainLOG, "time[ms]", 20, left);         //
		*mainLOG << std::endl;                                       //
	}                                                                //
	/*****************************************************************/
	// calculate potential only one time per point
	for (Int v = 0; v < Vsize; ++v) {
		/****************************LOG Section****************************************/
		if (printMainLogOn || mainLogStreamIsSet) {                                    //
			localTime = mainTimer.getElapsedTimeInMilliSec();                          //
			pointTime = 0;                                                             //
		}                                                                              //
		if (printPointLogOn) {                                                         //
            std::stringstream filename;                                                //
            if (isLogVol)  filename << "log/log_QECorr(V=" << (*data["V"])[v]();         //
        	if (isLogMdns) filename << "log/log_QECorr(D=" << (*data["D"])[v](gOverCmc); //
    		if (isLogVdns) filename << "log/log_QECorr(C=" << (*data["C"])[v]();         //
			filename << ", T=" << (*data["T"])[0]() << ")_";                           //
            filename << pointTimer.getCurrentDatetime();                               //
            filename << ".txt";                                                        //
            if (pointLOG->is_open()) pointLOG->close();                                //
            pointLOG->open(filename.str().c_str(), std::ios::out);                     //
		}                                                                              //
		if (printPointLogOn || pointLogStreamIsSet) {                                  //
			pointTimer.start();                                                        //
			*pointLOG << "Calculation LOG for point:" << std::endl;                    //
			*pointLOG << "Volume[Atomic,lin] = ";                                      //
			printer.printSciDouble(*pointLOG, (*data["V"])[v](), precision, 15, left); //
			*pointLOG << std::endl;                                                    //
			*pointLOG << "Mass density[g*cm^{-3},lin] = ";                             //
			printer.printSciDouble(*pointLOG, (*data["D"])[v](), precision, 15, left); //
			*pointLOG << std::endl;                                                    //
			*pointLOG << "Number density[cm^{-3},lin] = ";                             //
			printer.printSciDouble(*pointLOG, (*data["C"])[v](), precision, 15, left); //
			*pointLOG << std::endl;                                                    //
			*pointLOG << "Temperature[Hartree,lin] = ";                                //
			printer.printSciDouble(*pointLOG, (*data["T"])[0](), precision, 15, left); //
			*pointLOG << std::endl;                                                    //
			coldPsi.setLogStream(pointLOG);                                            //
		}                                                                              //
		/*******************************************************************************/
		coldPhi.setParameters((*(data["V"]))[v], coldT, Z);
		coldPsi.setParameters((*(data["V"]))[v], coldT, Z);
		for (Int i = inputVarBound; i < coldVarBound; ++i) {
			calculated[varName[i]] = false;
			if (need[varName[i]]) calculate(varName[i], v);
		}
		for (Int t = 0; t < Tsize; ++t) {
			/****************************LOG Section****************************************/
			if (t > 0 && (printMainLogOn || mainLogStreamIsSet))  {                        //
				localTime = mainTimer.getElapsedTimeInMilliSec();                          //
				pointTime = 0;                                                             //
			}                                                                              //
			if (t > 0 && printPointLogOn) {                                                //
            	std::stringstream filename;                                                //
            	if (isLogVol)  filename << "log/log_QECorr(V=" << (*data["V"])[v]();         //
        		if (isLogMdns) filename << "log/log_QECorr(D=" << (*data["D"])[v](gOverCmc); //
    			if (isLogVdns) filename << "log/log_QECorr(C=" << (*data["C"])[v]();         //
				filename << ", T=" << (*data["T"])[t]() << ")_";                           //
            	filename << pointTimer.getCurrentDatetime();                               //
            	filename << ".txt";                                                        //
            	if (pointLOG->is_open()) pointLOG->close();                                //
            	pointLOG->open(filename.str().c_str(), std::ios::out);                     //
			}                                                                              //
			if (t > 0 && (printPointLogOn || pointLogStreamIsSet)) {                       //
				pointTimer.start();                                                        //
				*pointLOG << "Calculation LOG for point:";                                 //
				*pointLOG << "Volume[Atomic,lin] = ";                                      //
				printer.printSciDouble(*pointLOG, (*data["V"])[v](), precision, 15, left); //
				*pointLOG << std::endl;                                                    //
				*pointLOG << "Mass density[g*cm^{-3},lin] = ";                             //
				printer.printSciDouble(*pointLOG, (*data["D"])[v](), precision, 15, left); //
				*pointLOG << std::endl;                                                    //
				*pointLOG << "Number density[cm^{-3},lin] = ";                             //
				printer.printSciDouble(*pointLOG, (*data["C"])[v](), precision, 15, left); //
				*pointLOG << std::endl;                                                    //
				*pointLOG << "Temperature[Hartree,lin] = ";                                //
				printer.printSciDouble(*pointLOG, (*data["T"])[t](), precision, 15, left); //
				*pointLOG << std::endl;                                                    //
			}                                                                              //
			if (printPointLogOn || pointLogStreamIsSet) {                                  //
				psi.setLogStream(pointLOG);                                                //
			}                                                                              //
			/*******************************************************************************/
			phi.setParameters((*(data["V"]))[v], (*(data["T"]))[t], Z);
			psi.setParameters((*(data["V"]))[v], (*(data["T"]))[t], Z);
			if (showProgress) {
				std::stringstream ss;
				int progress = (100*(v*Tsize + t + 1))/(Vsize*Tsize);
				ss << progress;
				std::cout << "\r";
				std::cout << "progress: ";
				std::cout << ss.str();
				std::cout << "\%";
			}
			for (Int i = coldVarBound; i < varNum; ++i) {                    
				calculated[varName[i]] = false;
				if (need[varName[i]]) calculate(varName[i], v, t);
			}
			/****************************LOG Section********************************************/
			if (printPointLogOn || pointLogStreamIsSet) {                                      //
				*pointLOG << "QECorr calculation for point is finished. Elapsed time: ";         //
				printer.printSciDouble(*pointLOG,                                              //
					pointTimer.getElapsedTimeInMilliSec(),                                     //
					precision, 15, left);                                                      //
				*pointLOG << " ms" << std::endl;                                               //
				pointTimer.reset();                                                            //
				pointLOG->close();                                                             //
			}                                                                                  //
			if (printMainLogOn || mainLogStreamIsSet) {                                        //
				pointTime = mainTimer.getElapsedTimeInMilliSec() - localTime;                  //
				for (Int i = 0; i < varNum; ++i) {                                             //
					if (need[varName[i]]) {                                                    //
						if (i < coldVarBound && varName[i][0] != 'T')                          //
							printer.printSciDouble(                                            //
								*mainLOG,                                                      //
								(*(data[varName[i]]))[v](unit[varName[i]], scale[varName[i]]), //
								precision, 20, left);                                          //
						if (varName[i][0] == 'T')                                              //
							printer.printSciDouble(                                            //
								*mainLOG,                                                      //
								(*(data[varName[i]]))[t](unit[varName[i]], scale[varName[i]]), //
								precision, 20, left);                                          //
						if (i >= coldVarBound)                                                 //
							printer.printSciDouble(                                            //
								*mainLOG,                                                      //
								(*(data[varName[i]]))[v*Tsize + t](unit[varName[i]],           //
								scale[varName[i]]),                                            //
								precision, 20, left);                                          //
					}                                                                          //
				}                                                                              //
				printer.printSciDouble(*mainLOG, phi(1), precision, 20, left);                 //
				printer.printSciDouble(*mainLOG, phi.derivative(0), precision, 20, left);      //
				printer.printSciDouble(*mainLOG, coldPhi(1), precision, 20, left);             //
				printer.printSciDouble(*mainLOG, coldPhi.derivative(0), precision, 20, left);  //
				printer.printSciDouble(*mainLOG, psi(1), precision, 20, left);                 //
				printer.printSciDouble(*mainLOG, psi.derivative(0), precision, 20, left);      //
				printer.printSciDouble(*mainLOG, coldPsi(1), precision, 20, left);             //
				printer.printSciDouble(*mainLOG, coldPsi.derivative(0), precision, 20, left);  //
				printer.printSciDouble(*mainLOG, pointTime, 6, 20, left);                      //
				*mainLOG << std::endl;                                                         //
			}                                                                                  //
			/***********************************************************************************/
		}
	}
	if (showProgress) {
		std::cout << "\r";
		std::cout << "progress: finished" << std::endl;
	}
}

PhysQvec& QECorr::operator[] (std::string quantity) { 
	if (!calculated[quantity]) calculate(quantity);
	return *(data[quantity]);
}

/////////////////////////////////////////////////////////////////////////////////////////////////

void QECorr::calculateDP(Int v, Int t) {
	Double Ptime;
	Double currentDP = FDhalf(phi(1));
	Double psi_1 = psi(1);
	if (printPointLogOn || pointLogStreamIsSet) {
		Ptime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating DP:" << std::endl;
		*pointLOG << "I_{1/2}(\\phi(1)) = ";
	}
	Double currentY = Y(phi(1));
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentDP, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "Y(\\phi(1)) = ";
		printer.printSciDouble(*pointLOG, currentDP, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "DP_{TF} = T^2/(3\\pi^3)*(I_{1/2}(\\phi(1))*\\psi(1) + Y(\\phi(1))) = ";
	}
	currentDP = pow((*(data["T"]))[t](), 2.0)/3.0/M_PI/M_PI/M_PI*(currentDP*psi_1 + currentY);
	(*(data["DP"]))[v*Tsize + t].setValue(currentDP);
	calculated["DP"] = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentDP, precision, 20, left);
		*pointLOG << std::endl;
		Ptime = pointTimer.getElapsedTimeInMilliSec() - Ptime;
		*pointLOG << "Elapsed time = " << Ptime << " ms" << std::endl;
	}
}
void QECorr::calculateDPT(Int v, Int t) {
	if (!calculated["DP"])  calculateDP(v, t);
	if (!calculated["DPC"]) calculateDPC(v);
	Double PTtime;
	if (printPointLogOn || pointLogStreamIsSet) {
		PTtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating DPT:" << std::endl;
		*pointLOG << "DPT = DP - DPC = ";
	}
	Double currentDPT;
	currentDPT = (*(data["DP"]))[v*Tsize + t]() - (*(data["DPC"]))[v]();	
	(*(data["DPT"]))[v*Tsize + t].setValue(currentDPT);
	calculated["DPT"] = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentDPT, precision, 20, left);
		*pointLOG << std::endl;
		PTtime = pointTimer.getElapsedTimeInMilliSec() - PTtime;
		*pointLOG << "Elapsed time = " << PTtime << " ms" << std::endl;
	}
}
void QECorr::calculateDPC(Int v) {
	Double PCtime;
	Double currentDPC = FDhalf(coldPhi(1));
	Double coldPsi_1 = coldPsi(1);
	if (printPointLogOn || pointLogStreamIsSet) {
		PCtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating DPC:" << std::endl;
		*pointLOG << "T_{C}[Hartree] = ";
		printer.printSciDouble(*pointLOG, coldT(), precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "I_{1/2}(\\phi_{C}(1)) = ";
	}
	Double currentY = Y(coldPhi(1));
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentDPC, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "Y(\\phi_{C}(1)) = ";
		printer.printSciDouble(*pointLOG, currentDPC, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "DPC_{TF} = T^2/(3\\pi^3)*(I_{1/2}(\\phi_{C}(1))*\\psi_{C}(1) + Y(\\phi_{C}(1))) = ";
	}
	currentDPC = pow(2*coldT(), 5.0/2.0)/6.0/M_PI/M_PI*currentDPC;
	(*(data["DPC"]))[v].setValue(currentDPC);
	calculated["DPC"] = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentDPC, precision, 20, left);
		*pointLOG << std::endl;
		PCtime = pointTimer.getElapsedTimeInMilliSec() - PCtime;
		*pointLOG << "Elapsed time = " << PCtime << " ms" << std::endl;
	}
}
void QECorr::calculateDE(Int v, Int t) {
	static const Double DE0 = 0.269900170;
	Double Etime;
	Double currentDE;
    DoubleVec startDE(0.0, rhsQECorrEnergy::dim);
    Double a, b;
    Double V1, T1;
    Double xFrom = 1.0;
    Double xTo = 0.0;
    Double phi_1 = phi(1);
    Double psi_1 = psi(1);

    if (printPointLogOn || pointLogStreamIsSet) {
		Etime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating DE:" << std::endl;
		*pointLOG << "Transform parameters " << std::endl;
	}

	V1 = (*(data["V"]))[v]()*Z;
	T1 = (*(data["T"]))[t]()*pow(Z, -4.0/3.0);

    a =   pow(2.0, 7.0/6.0)
        * pow(3.0, 2.0/3.0)
        * pow(M_PI, -5.0/3.0)
        * sqrt(T1)*pow(V1, 2.0/3.0);
    b = V1*pow(T1, 2.0)/M_PI/M_PI/M_PI;

	if (printPointLogOn || pointLogStreamIsSet) {
		*pointLOG << "V_1 = V_Z*Z = ";
		printer.printSciDouble(*pointLOG, V1, precision, 20, left);
		*pointLOG << std::endl << "T_1 = T_Z*Z^{-4/3}";
		printer.printSciDouble(*pointLOG, T1, precision, 20, left);
		*pointLOG << std::endl << "DE_0 = 0.269900170" << std::endl;
		*pointLOG << "a = 2^{7/6}*3^{2/3}*\\pi^{-5/3}*\\T_1^{1/2}*V_1^{2/3} = ";
		printer.printSciDouble(*pointLOG, a, precision, 20, left);
		*pointLOG << std::endl << "b = V_1*T_1^2*\\pi^{3} = ";
		printer.printSciDouble(*pointLOG, b, precision, 20, left);
		*pointLOG << std::endl;
	}

    startDE[0] = phi_1;
	startDE[1] = phi_1;
	startDE[2] = psi_1;
    startDE[3] = psi_1;
    startDE[4] = sqrt(0.5*T1)/3.0/M_PI*psi.derivative(0) + DE0;

    if (printPointLogOn || pointLogStreamIsSet) {
    	*pointLOG << "Initial vector for ODE" << std::endl;
    	*pointLOG << "y_0(1) = \\phi(1) = ";
    	printer.printSciDouble(*pointLOG, startDE[0], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_1(1) = \\phi(1) = ";
    	printer.printSciDouble(*pointLOG, startDE[1], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_2(1) = \\psi(1) = ";
    	printer.printSciDouble(*pointLOG, startDE[2], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_3(1) = \\psi(1) = ";
    	printer.printSciDouble(*pointLOG, startDE[3], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_4(1) = 2^{1/2}*3^{-1}*\\pi^{-1}*T_{1}^{1/2}*\\psi(1) + DE_0 = ";
    	printer.printSciDouble(*pointLOG, startDE[4], precision, 20, left);
    	*pointLOG << std::endl;
    }

	rhsEnergy.updateParameters(a, b);
	energySolver.Integrate(rhsEnergy, startDE, xFrom, xTo);
	currentDE = startDE[4];

	(*(data["DE"]))[v*Tsize + t].setValue(currentDE*pow(Z, 5.0/3.0));

	calculated["DE"]  = true;

	if (printPointLogOn || pointLogStreamIsSet) {
		Int size = energyData.Count();
        printer.printString((*pointLOG), "x");
        printer.printString((*pointLOG), "phi(x)");
        printer.printString((*pointLOG), "dphi(x)");
        printer.printString((*pointLOG), "psi(x)");
        printer.printString((*pointLOG), "dpsi(x)");
        printer.printString((*pointLOG), "DE(x)");
        (*pointLOG) << std::endl;
        for (Int i = 0; i < size; ++i) {
            printer.printSciDouble((*pointLOG), energyData.xSave[i], precision);
            printer.printSciDouble((*pointLOG), energyData.ySave[0][i], precision);
            printer.printSciDouble((*pointLOG), energyData.ySave[1][i], precision);
            printer.printSciDouble((*pointLOG), energyData.ySave[2][i], precision);
            printer.printSciDouble((*pointLOG), energyData.ySave[3][i], precision);
            printer.printSciDouble((*pointLOG), energyData.ySave[4][i], precision);
            (*pointLOG) << std::endl;
        }
        *pointLOG << "Calculated energy DE_1 = DE(0) = ";
        printer.printSciDouble((*pointLOG), currentDE, precision);
        *pointLOG << std::endl;

		*pointLOG << "Transform to DE_Z = DE_1*Z^{5/3} = ";
		printer.printSciDouble((*pointLOG), currentDE*pow(Z, 5.0/3.0), precision);
        *pointLOG << std::endl;

        Etime = pointTimer.getElapsedTimeInMilliSec() - Etime;
        *pointLOG << "Elapsed time = " << Etime << " ms" << std::endl;
    }

}
void QECorr::calculateDET(Int v, Int t) {
	if (!calculated["DE"])  calculateDE(v, t);
	if (!calculated["DEC"]) calculateDEC(v);
	Double ETtime;
	if (printPointLogOn || pointLogStreamIsSet) {
		ETtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating DET:" << std::endl;
		*pointLOG << "DET = DE - DEC = ";
	}
	Double currentDET;
	currentDET = (*(data["DE"]))[v*Tsize + t]() - (*(data["DEC"]))[v]();
	(*(data["DET"]))[v*Tsize + t].setValue(currentDET);
	calculated["DET"] = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentDET, precision, 20, left);
		*pointLOG << std::endl;
		ETtime = pointTimer.getElapsedTimeInMilliSec() - ETtime;
		*pointLOG << "Elapsed time = " << ETtime << " ms" << std::endl;
	}
}
void QECorr::calculateDEC(Int v) {
	static const Double DE0 = 0.269900170;
	Double ECtime;
	Double currentDEC;
	DoubleVec startDEC(0.0, rhsQECorrEnergy::dim);
	Double a, b;
	Double V1, T1;
	Double xFrom = 1.0;
	Double xTo = 0.0;
	Double coldPhi_1 = coldPhi(1);
	Double coldPsi_1 = coldPsi(1);

	if (printPointLogOn || pointLogStreamIsSet) {
		ECtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating DEC:" << std::endl;
		*pointLOG << "Transform parameters " << std::endl;
	}
	
	T1 = coldT()*pow(Z, -4.0/3.0);
	V1 = (*(data["V"]))[v]()*Z;
	
	a =   pow(2.0, 7.0/6.0)
		* pow(3.0, 2.0/3.0)
		* pow(M_PI, -5.0/3.0)
		* sqrt(T1)*pow(V1, 2.0/3.0);
	b = V1*pow(T1, 2.0)/M_PI/M_PI/M_PI;

	if (printPointLogOn || pointLogStreamIsSet) {
		*pointLOG << "V_1 = V_Z*Z = ";
		printer.printSciDouble(*pointLOG, V1, precision, 20, left);
		*pointLOG << std::endl << "TC_1 = TC_Z*Z^{-4/3}";
		printer.printSciDouble(*pointLOG, T1, precision, 20, left);
		*pointLOG << std::endl << "DE_0 = 0.269900170" << std::endl;
		*pointLOG << "a = 2^{7/6}*3^{2/3}*\\pi^{-5/3}*\\TC_1^{1/2}*V_1^{2/3} = ";
		printer.printSciDouble(*pointLOG, a, precision, 20, left);
		*pointLOG << std::endl << "b = V_1*TC_1^2*\\pi^{3} = ";
		printer.printSciDouble(*pointLOG, b, precision, 20, left);
		*pointLOG << std::endl;
	}

	startDEC[0] = coldPhi_1;
	startDEC[1] = coldPhi_1;
	startDEC[2] = coldPsi_1;
    startDEC[3] = coldPsi_1;
    startDEC[4] = sqrt(0.5*T1)/3.0/M_PI*coldPsi.derivative(0) + DE0;

    if (printPointLogOn || pointLogStreamIsSet) {
    	*pointLOG << "Initial vector for ODE" << std::endl;
    	*pointLOG << "y_0(1) = \\phi_{C}(1) = ";
    	printer.printSciDouble(*pointLOG, startDEC[0], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_1(1) = \\phi_{C}(1) = ";
    	printer.printSciDouble(*pointLOG, startDEC[1], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_2(1) = \\psi_{C}(1) = ";
    	printer.printSciDouble(*pointLOG, startDEC[2], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_3(1) = \\psi_{C}(1) = ";
    	printer.printSciDouble(*pointLOG, startDEC[3], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_4(1) = 2^{1/2}*3^{-1}*\\pi^{-1}*TC_{1}^{1/2}*\\psi_{C}(1) + DE_0 = ";
    	printer.printSciDouble(*pointLOG, startDEC[4], precision, 20, left);
    	*pointLOG << std::endl;
    }

	rhsColdEnergy.updateParameters(a, b);

	coldEnergySolver.Integrate(rhsColdEnergy, startDEC, xFrom, xTo);
	currentDEC = startDEC[4];
	(*(data["DEC"]))[v].setValue(currentDEC*pow(Z, 5.0/3.0));

	calculated["DEC"] = true;

	if (printPointLogOn || pointLogStreamIsSet) {
		Int size = coldEnergyData.Count();
        printer.printString((*pointLOG), "x");
        printer.printString((*pointLOG), "phi_{C}(x)");
        printer.printString((*pointLOG), "dphi_{C}(x)");
        printer.printString((*pointLOG), "psi_{C}(x)");
        printer.printString((*pointLOG), "dpsi_{C}(x)");
        printer.printString((*pointLOG), "DEC(x)");
        (*pointLOG) << std::endl;
        for (Int i = 0; i < size; ++i) {
            printer.printSciDouble((*pointLOG), coldEnergyData.xSave[i], precision);
            printer.printSciDouble((*pointLOG), coldEnergyData.ySave[0][i], precision);
            printer.printSciDouble((*pointLOG), coldEnergyData.ySave[1][i], precision);
            printer.printSciDouble((*pointLOG), coldEnergyData.ySave[2][i], precision);
            printer.printSciDouble((*pointLOG), coldEnergyData.ySave[3][i], precision);
            printer.printSciDouble((*pointLOG), coldEnergyData.ySave[4][i], precision);
            (*pointLOG) << std::endl;
        }
        *pointLOG << "Calculated energy DEC_1 = DEC(0) = ";
        printer.printSciDouble((*pointLOG), currentDEC, precision);
        *pointLOG << std::endl;

		*pointLOG << "Transform to DEC_Z = DEC_1*Z^{5/3} = ";
		printer.printSciDouble((*pointLOG), currentDEC*pow(Z, 5.0/3.0), precision);
        *pointLOG << std::endl;

        ECtime = pointTimer.getElapsedTimeInMilliSec() - ECtime;
        *pointLOG << "Elapsed time = " << ECtime << " ms" << std::endl;
    }
}
void QECorr::calculateDS(Int v, Int t) {
	Double currentDS;
	Double Stime;
    DoubleVec startDS(0.0, rhsQECorrEntropy::dim);
    Double a, b;
    Double V1, T1;
    Double xFrom = 1.0;
    Double xTo = 0.0;
    Double phi_1 = phi(1);
    Double psi_1 = psi(1);

    if (printPointLogOn || pointLogStreamIsSet) {
		Stime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating DS:" << std::endl;
		*pointLOG << "Transform parameters " << std::endl;
	}

	V1 = (*(data["V"]))[v]()*Z;
	T1 = (*(data["T"]))[t]()*pow(Z, -4.0/3.0);

    a =   pow(2.0, 7.0/6.0)
        * pow(3.0, 2.0/3.0)
        * pow(M_PI, -5.0/3.0)
        * sqrt(T1)*pow(V1, 2.0/3.0);
    b = T1*V1/M_PI/M_PI/M_PI;

	if (printPointLogOn || pointLogStreamIsSet) {
		*pointLOG << "V_1 = V_Z*Z = ";
		printer.printSciDouble(*pointLOG, V1, precision, 20, left);
		*pointLOG << std::endl << "T_1 = T_Z*Z^{-4/3}";
		printer.printSciDouble(*pointLOG, T1, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "a = 2^{7/6}*3^{2/3}*\\pi^{-5/3}*\\T_1^{1/2}*V_1^{2/3} = ";
		printer.printSciDouble(*pointLOG, a, precision, 20, left);
		*pointLOG << std::endl << "b = V_1*T_1*\\pi^{-3} = ";
		printer.printSciDouble(*pointLOG, b, precision, 20, left);
		*pointLOG << std::endl;
	}
	
	startDS[0] = phi_1;
	startDS[1] = phi_1;
	startDS[2] = psi_1;
	startDS[3] = psi_1;
	startDS[4] = 1.0/(sqrt(2.0*T1)*3.0*M_PI)*psi.derivative(0);

	if (printPointLogOn || pointLogStreamIsSet) {
    	*pointLOG << "Initial vector for ODE" << std::endl;
    	*pointLOG << "y_0(1) = \\phi(1)";
    	printer.printSciDouble(*pointLOG, startDS[0], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_1(1) = \\phi(1)";
    	printer.printSciDouble(*pointLOG, startDS[1], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_2(1) = \\psi(1) = ";
    	printer.printSciDouble(*pointLOG, startDS[2], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_3(1) = \\psi(1) = ";
    	printer.printSciDouble(*pointLOG, startDS[3], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_4(1) = (3*\\pi*(2*T_1)^{1/2})^{-1}*psi'(1) = ";
    	printer.printSciDouble(*pointLOG, startDS[4], precision, 20, left);
    	*pointLOG << std::endl;
    }

	rhsEntropy.updateParameters(a, b);

	entropySolver.Integrate(rhsEntropy, startDS, xFrom, xTo);
	currentDS = startDS[4];
	(*(data["DS"]))[v*Tsize + t].setValue(currentDS*pow(Z, 1.0/3.0));
    	
	calculated["DS"]  = true;

	if (printPointLogOn || pointLogStreamIsSet) {
		Int size = energyData.Count();
        printer.printString((*pointLOG), "x");
        printer.printString((*pointLOG), "phi(x)");
        printer.printString((*pointLOG), "dphi(x)");
        printer.printString((*pointLOG), "psi(x)");
        printer.printString((*pointLOG), "dpsi(x)");
        printer.printString((*pointLOG), "DS(x)");
        (*pointLOG) << std::endl;
        for (Int i = 0; i < size; ++i) {
            printer.printSciDouble((*pointLOG), entropyData.xSave[i], precision);
            printer.printSciDouble((*pointLOG), entropyData.ySave[0][i], precision);
            printer.printSciDouble((*pointLOG), entropyData.ySave[1][i], precision);
            printer.printSciDouble((*pointLOG), entropyData.ySave[2][i], precision);
            printer.printSciDouble((*pointLOG), entropyData.ySave[3][i], precision);
            printer.printSciDouble((*pointLOG), entropyData.ySave[4][i], precision);
            (*pointLOG) << std::endl;
        }
        *pointLOG << "Calculated entropy DS_1 = DS(0) = ";
        printer.printSciDouble((*pointLOG), currentDS, precision);
        *pointLOG << std::endl;

		*pointLOG << "Transform to DS_Z = DS_1*Z^{1/3} = ";
		printer.printSciDouble((*pointLOG), currentDS*pow(Z, 1.0/3.0), precision);
        *pointLOG << std::endl;

        Stime = pointTimer.getElapsedTimeInMilliSec() - Stime;
        *pointLOG << "Elapsed time = " << Stime << " ms" << std::endl;
    }

}
void QECorr::calculateDST(Int v, Int t) {
	if (!calculated["DS"])  calculateDS(v, t);
	if (!calculated["DSC"]) calculateDSC(v);
	Double STtime;
	if (printPointLogOn || pointLogStreamIsSet) {
		STtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating DST:" << std::endl;
		*pointLOG << "DST = DS - DSC = ";
	}
	Double currentDST;
	currentDST = (*(data["DS"]))[v*Tsize + t]() - (*(data["DSC"]))[v]();
	(*(data["DST"]))[v*Tsize + t].setValue(currentDST);
	
	calculated["DST"] = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentDST, precision, 20, left);
		*pointLOG << std::endl;
		STtime = pointTimer.getElapsedTimeInMilliSec() - STtime;
		*pointLOG << "Elapsed time = " << STtime << " ms" << std::endl;
	}
}
void QECorr::calculateDSC(Int v) {
	Double SCtime;
	Double currentDSC;
	DoubleVec startDSC(0.0, rhsQECorrEntropy::dim);
	Double a, b;
	Double V1, T1;
	Double xFrom = 1.0;
	Double xTo = 0.0;
	Double coldPhi_1 = coldPhi(1);
	Double coldPsi_1 = coldPsi(1);

	if (printPointLogOn || pointLogStreamIsSet) {
		SCtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating DSC:" << std::endl;
		*pointLOG << "Transform parameters " << std::endl;
	}
	
	T1 = coldT()*pow(Z, -4.0/3.0);
	V1 = (*(data["V"]))[v]()*Z;
	a =   pow(2.0, 7.0/6.0)
	    * pow(3.0, 2.0/3.0)
	    * pow(M_PI, -5.0/3.0)
	    * sqrt(T1)*pow(V1, 2.0/3.0);
	b = T1*V1/M_PI/M_PI/M_PI;
	
	if (printPointLogOn || pointLogStreamIsSet) {
		*pointLOG << "V_1 = V_Z*Z = ";
		printer.printSciDouble(*pointLOG, V1, precision, 20, left);
		*pointLOG << std::endl << "TC_1 = TC_Z*Z^{-4/3}";
		printer.printSciDouble(*pointLOG, T1, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "a = 2^{7/6}*3^{2/3}*\\pi^{-5/3}*\\TC_1^{1/2}*V_1^{2/3} = ";
		printer.printSciDouble(*pointLOG, a, precision, 20, left);
		*pointLOG << std::endl << "b = V_1*TC_1*\\pi^{-3} = ";
		printer.printSciDouble(*pointLOG, b, precision, 20, left);
		*pointLOG << std::endl;
	}
	
	startDSC[0] = coldPhi_1;
	startDSC[1] = coldPhi_1;
	startDSC[2] = coldPsi_1;
	startDSC[3] = coldPsi_1;
	startDSC[4] = 1.0/(sqrt(2.0*T1)*3.0*M_PI)*coldPsi.derivative(0);

	if (printPointLogOn || pointLogStreamIsSet) {
    	*pointLOG << "Initial vector for ODE" << std::endl;
    	*pointLOG << "y_0(1) = \\phi_{C}(1)";
    	printer.printSciDouble(*pointLOG, startDSC[0], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_1(1) = \\phi_{C}(1)";
    	printer.printSciDouble(*pointLOG, startDSC[1], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_2(1) = \\psi_{C}(1) = ";
    	printer.printSciDouble(*pointLOG, startDSC[2], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_3(1) = \\psi_{C}(1) = ";
    	printer.printSciDouble(*pointLOG, startDSC[3], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_4(1) = (3*\\pi*(2*TC_1)^{1/2})^{-1}*psi'_{C}(1) = ";
    	printer.printSciDouble(*pointLOG, startDSC[4], precision, 20, left);
    	*pointLOG << std::endl;
    }

	rhsColdEntropy.updateParameters(a, b);

	coldEntropySolver.Integrate(rhsColdEntropy, startDSC, xFrom, xTo);
	currentDSC = startDSC[4];
	(*(data["DSC"]))[v].setValue(currentDSC*pow(Z, 1.0/3.0));

	calculated["DSC"] = true;

	if (printPointLogOn || pointLogStreamIsSet) {
		Int size = coldEnergyData.Count();
        printer.printString((*pointLOG), "x");
        printer.printString((*pointLOG), "phi_{C}(x)");
        printer.printString((*pointLOG), "dphi_{C}(x)");
        printer.printString((*pointLOG), "psi_{C}(x)");
        printer.printString((*pointLOG), "dpsi_{C}(x)");
        printer.printString((*pointLOG), "DSC(x)");
        (*pointLOG) << std::endl;
        for (Int i = 0; i < size; ++i) {
            printer.printSciDouble((*pointLOG), coldEntropyData.xSave[i], precision);
            printer.printSciDouble((*pointLOG), coldEntropyData.ySave[0][i], precision);
            printer.printSciDouble((*pointLOG), coldEntropyData.ySave[1][i], precision);
            printer.printSciDouble((*pointLOG), coldEntropyData.ySave[2][i], precision);
            printer.printSciDouble((*pointLOG), coldEntropyData.ySave[3][i], precision);
            printer.printSciDouble((*pointLOG), coldEntropyData.ySave[4][i], precision);
            (*pointLOG) << std::endl;
        }
        *pointLOG << "Calculated entropy DSC_1 = DSC(0) = ";
        printer.printSciDouble((*pointLOG), currentDSC, precision);
        *pointLOG << std::endl;

		*pointLOG << "Transform to DSC_Z = DSC_1*Z^{1/3} = ";
		printer.printSciDouble((*pointLOG), currentDSC*pow(Z, 1.0/3.0), precision);
        *pointLOG << std::endl;

        SCtime = pointTimer.getElapsedTimeInMilliSec() - SCtime;
        *pointLOG << "Elapsed time = " << SCtime << " ms" << std::endl;
    }
}
void QECorr::calculateDM(Int v, Int t) {
	Double Mtime;
	Double phi_1 = phi(1);
	Double psi_1 = psi(1);
	if (printPointLogOn || pointLogStreamIsSet) {
		Mtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating DM:" << std::endl;
	}

	Double currentDM = FDmhalf(phi_1);
	if (printPointLogOn || pointLogStreamIsSet) {
		*pointLOG << "I_{1/2}(\\phi(1)) = ";
		printer.printSciDouble(*pointLOG, currentDM, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "DM_{TF} = (0.5*T)^{1/2}/(3*\\pi)*(0.5*\\I_{1/2}(phi(1)) + \\psi(1)) = ";
	}
	currentDM = sqrt(0.5*(*(data["T"]))[t]())/3.0/M_PI
        	  * (0.5*currentDM + psi_1);
	(*(data["DM"]))[v*Tsize + t].setValue(currentDM);
	calculated["DM"]  = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentDM, precision, 20, left);
		*pointLOG << std::endl;
		Mtime = pointTimer.getElapsedTimeInMilliSec() - Mtime;
		*pointLOG << "Elapsed time = " << Mtime << " ms" << std::endl;
	}
}
void QECorr::calculateDMT(Int v, Int t) {
	if (!calculated["DM"]) calculateDM(v, t);
	if (!calculated["DMC"]) calculateDMC(v);
	Double MTtime;
	if (printPointLogOn || pointLogStreamIsSet) {
		MTtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating DMT:" << std::endl;
		*pointLOG << "DMT = DM - DMC = ";
	}
	
	Double currentDMT;
	currentDMT = (*(data["DM"]))[v*Tsize + t]() - (*(data["DMC"]))[v]();
	(*(data["DMT"]))[v*Tsize + t].setValue(currentDMT);
	calculated["DMT"] = true;

	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentDMT, precision, 20, left);
		*pointLOG << std::endl;
		MTtime = pointTimer.getElapsedTimeInMilliSec() - MTtime;
		*pointLOG << "Elapsed time = " << MTtime << " ms" << std::endl;
	}
}
void QECorr::calculateDMC(Int v) {
	Double MCtime;
	Double coldPhi_1 = coldPhi(1);
	Double coldPsi_1 = coldPsi(1);
	if (printPointLogOn || pointLogStreamIsSet) {
		MCtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating DMC:" << std::endl;
	}
	Double currentDMC = FDmhalf(coldPhi_1);
	if (printPointLogOn || pointLogStreamIsSet) {
		*pointLOG << "I_{1/2}(\\phi_{C}(1)) = ";
		printer.printSciDouble(*pointLOG, currentDMC, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "DMC_{TF} = (0.5*TC)^{1/2}/(3*\\pi)*(0.5*\\I_{1/2}(phi_{C}(1)) + \\psi_{C}(1)) = ";
	}
	currentDMC = sqrt(0.5*coldT())/3.0/M_PI
        	  * (0.5*currentDMC + coldPsi_1);
	(*(data["DMC"]))[v].setValue(currentDMC);
	calculated["DMC"] = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentDMC, precision, 20, left);
		*pointLOG << std::endl;
		MCtime = pointTimer.getElapsedTimeInMilliSec() - MCtime;
		*pointLOG << "Elapsed time = " << MCtime << " ms" << std::endl;
	}
}