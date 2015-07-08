#pragma once
#include "../lib/numeric/hdr/ODE/ODEsolver.h"
#include "../lib/numeric/hdr/ODE/ODEdata.h"
#include "../lib/numeric/hdr/ODE/steppers/ODEstepperPD853.h"
#include "../lib/numeric/hdr/numericTypes.h"
#include "../lib/numeric/hdr/mathUtils.h"
#include "../lib/numeric/hdr/SF/FermiDirac.h"
#include "../hdr/ThermodynamicFunction.h"
#include "../hdr/FTTFpotential.h"
#include "../hdr/Yfunction.h"
#include "../hdr/Units.h"
#include "../hdr/Timer.h"
#include "../hdr/Printer.h"
#include "../hdr/stringUtils.h"
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>

#define DEFAULT_INPUT_FILE "in/FTTFinput.dat" 
#define DEFAULT_OUTPUT_FILE "out/FTTFoutput.dat" 
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
class FTTF {
public:
	/**
	* @brief A constructor of Thomas-Fermi model.
	*/
    FTTF(Double _Z = 1.0, Double _Mass = 1.0);
    ~FTTF();
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
	* @details The calculation goes through the solving of the ODE with parameters 
	* @f$ a = 2^{7/6}3^{2/3}\pi^{-5/3}T^{1/2}V^{2/3} @f$,
	* @f$ b = 3\sqrt{2}\pi^{-2}VT^{5/2} @f$,
	* @f$ E_0 = -0.76874512422 @f$. Thomas-Fermi energy is solution of the ODE:
	* @f{eqnarray*}{
	*  \frac{d^2\phi}{dx^2} &=& a x I_{1/2}
	*	  \left(
	*		\frac{\phi(x)}{x}
	*	  \right), \\
	*  E'(x) &=& b x^2 I_{3/2}
	*				\left(
	*				  \frac{\phi(x)}{x}
	*				\right), \\
	* \phi(1) &=& \phi'(1), \\
	* E(1) &=& \frac{2\sqrt{2}VT^{5/2}}{\pi^2}I_{3/2}(\phi(1)) - E_0,
	* @f}
	* Finally, @f$ E = E(0) @f$.
	*/
    struct rhsFTTFenergy {
		Double a;
		Double b;
		static const Int dim = 3;
		FermiDirac<Half> FDhalf;
		FermiDirac<ThreeHalf> FD3half;
		rhsFTTFenergy(Double _a = 0, Double _b = 0) : a(_a), b(_b) {}
		void updateParameters(const Double _a, const Double _b) { a = _a; b = _b; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
			dydx[0] = y[1];
    		if (x > 0) {
    			dydx[1] = a*x*FDhalf(y[0]/x);
    			dydx[2] = b*FD3half(y[0]/x)*x*x;
    		}
    		else { 
    			dydx[1] = dydx[2] = 1e+10;
    		}
		}
	};
	/**
	* @brief Problem for calculating entropy.
	* @details The calculation goes through the solving of the ODE with parameters
	* @f$ a = 2^{7/6}3^{2/3}\pi^{-5/3}T^{1/2}V^{2/3} @f$,
	* @f$ b = 7\sqrt{2}\pi^{-2}VT^{3/2} @f$,
	* Thomas-Fermi entropy is solution of the ODE:
	* @f{eqnarray*}{
	*  \frac{d^2\phi}{dx^2} &=& a x I_{1/2}
	*	  \left(
	*		\frac{\phi(x)}{x}
	*	  \right), \\
	*  S'(x) &=& b x^2 I_{3/2}
	*				\left(
	*				  \frac{\phi(x)}{x}
	*				\right), \\
	* \phi(1) &=& \phi'(1), \\
	* S(1) &=& \frac{4\sqrt{2}VT^{3/2}}{\pi^2}I_{3/2}(\phi(1)) - \phi'(0),
	* @f}
	* Finally, @f$ S = S(0) @f$.
	*/
	struct rhsFTTFentropy {
		Double a;
		Double b;
		static const Int dim = 3;
		FermiDirac<Half> FDhalf;
		FermiDirac<ThreeHalf> FD3half;
		rhsFTTFentropy(Double _a = 0, Double _b = 0) : a(_a), b(_b) {}
		void updateParameters(const Double _a, const Double _b) { a = _a; b = _b; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
			dydx[0] = y[1];
    		if (x > 0) {
    			dydx[1] = a*x*FDhalf(y[0]/x);
    			dydx[2] = b*FD3half(y[0]/x)*x*x;
    		}
    		else { 
    			dydx[1] = dydx[2] = 1e+10;
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
	void  calculate(std::string var, Int v, Int t);
	void    calculateP(Int v, Int t); 
	void   calculatePT(Int v, Int t); 
	void   calculatePC(Int v); 
	void    calculateE(Int v, Int t); 
	void   calculateET(Int v, Int t); 
	void   calculateEC(Int v); 
	void    calculateS(Int v, Int t); 
	void   calculateST(Int v, Int t); 
	void   calculateSC(Int v); 
	void    calculateM(Int v, Int t); 
	void   calculateMT(Int v, Int t); 
	void   calculateMC(Int v); 

	void calculateAll(); // uses fast calculation when prepare output
	void performOutput();
	void clearData();
	
    // ode solvers
    rhsFTTFenergy rhsEnergy;
    ODEsolver<ODEstepperPD853<rhsFTTFenergy> > energySolver;
    ODEdata energyData;

    rhsFTTFentropy rhsEntropy;
    ODEsolver<ODEstepperPD853<rhsFTTFentropy> > entropySolver;
    ODEdata entropyData;

    rhsFTTFenergy rhsColdEnergy;
    ODEsolver<ODEstepperPD853<rhsFTTFenergy> > coldEnergySolver;
    ODEdata coldEnergyData;

    rhsFTTFentropy rhsColdEntropy;
    ODEsolver<ODEstepperPD853<rhsFTTFentropy> > coldEntropySolver;
    ODEdata coldEntropyData;

    //transform quantities
    void transformDtoV();
    void transformCtoV();
    void transformVtoD();
    void transformVtoC();

    FermiDirac<ThreeHalf> FD3half;
    PhysQ coldT;

    FTTFpotential phi;
    FTTFpotential coldPhi;

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

FTTF::FTTF(Double _Z, Double _Mass) : 
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
	varName[4] = "PC"; varName[5] = "EC"; varName[6] = "SC"; varName[7] = "MC";
	coldVarBound = 8;
	varName[8] = "P"; varName[9] = "E"; varName[10] = "S"; varName[11] = "M";
	varName[12] = "PT"; varName[13] = "ET"; varName[14] = "ST"; varName[15] = "MT";
	
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

FTTF::~FTTF() {
	clearData();
	setPrintMainLogOff();
	setPrintPointLogOff();
}

void FTTF::calculate(std::string quantity, Int v = -1, Int t = -1) {
	if (v < 0 && t < 0) {
		calculateData(quantity);
	}
	if (v >= 0 && t < 0) {
			 if (!quantity.compare("PC")) calculatePC(v);
		else if (!quantity.compare("EC")) calculateEC(v);
		else if (!quantity.compare("SC")) calculateSC(v);
		else if (!quantity.compare("MC")) calculateMC(v);
	}
	if (v >= 0 && t >= 0) {
			 if (!quantity.compare("P"))  calculateP(v, t);
		else if (!quantity.compare("E"))  calculateE(v, t);
		else if (!quantity.compare("S"))  calculateS(v, t);
		else if (!quantity.compare("M"))  calculateM(v, t);
		else if (!quantity.compare("PT")) calculatePT(v, t);
		else if (!quantity.compare("ET")) calculateET(v, t);
		else if (!quantity.compare("ST")) calculateST(v, t);
		else if (!quantity.compare("MT")) calculateMT(v, t);
	}
}

void FTTF::setTolerance(const Double _eps) { 
	eps = _eps;
	precision = static_cast<int>(-log10(eps));
	if (printMainLogOn || mainLogStreamIsSet) {
		*mainLOG << "FTTF accepted new tolerance, eps = ";
		printer.printSciDouble(*mainLOG, eps, precision);
		*mainLOG << std::endl;
		phi.setLogStream(mainLOG);
		coldPhi.setLogStream(mainLOG);
	}
	phi.setTolerance(eps);
	if (printMainLogOn || mainLogStreamIsSet) *mainLOG << "Cold ";
	coldPhi.setTolerance(eps);
	energySolver.SetTolerance(0.0, eps/10);
	coldEnergySolver.SetTolerance(0.0, eps/10);
	entropySolver.SetTolerance(0.0, eps/10);
	coldEntropySolver.SetTolerance(0.0, eps/10);
	if (printMainLogOn || mainLogStreamIsSet) {
		phi.clearLogStream();
		coldPhi.clearLogStream();	
	}
}

void FTTF::setMainLogStream(std::ofstream* _mainLOG) {
    if (printMainLogOn) setPrintMainLogOff();
    mainLOG = _mainLOG;
    mainTimer.start();
    mainLogStreamIsSet = true;
}

void FTTF::setPrintMainLogOn() {
    if (!mainLogStreamIsSet) {
        mainLOG = new std::ofstream;
        mainTimer.start();
        std::stringstream filename;
        filename << "log/log_FTTF_Z(";
        filename << Z << ")_M("; 
        filename << Mass << ")_";
        filename << mainTimer.getCurrentDatetime();
        filename << ".txt";
        mainLOG->open(filename.str().c_str(), std::ios::out);
        *mainLOG << "FTTF log is started" << std::endl;
        *mainLOG << "Z = " << Z << std::endl;
        *mainLOG << "Atomic Mass = " << Mass << std::endl;
        printMainLogOn = true;
    }
}

void FTTF::setPrintMainLogOff() {
    if (!mainLogStreamIsSet) {
        if (printMainLogOn) {
            printMainLogOn = false;
            if (mainLOG->is_open()) mainLOG->close();
            delete mainLOG;
            mainLOG = NULL;
        }
    }
}

void FTTF::setPointLogStream(std::ofstream* _pointLOG) {
    if (printPointLogOn) setPrintPointLogOff();
    pointLOG = _pointLOG;
    pointLogStreamIsSet = true;
}

void FTTF::setPrintPointLogOn() {
    if (!pointLogStreamIsSet) {
        pointLOG = new std::ofstream;
        printPointLogOn = true;
    }
}

void FTTF::setPrintPointLogOff() {
    if (!pointLogStreamIsSet) {
        if (printPointLogOn) {
            printPointLogOn = false;
            if (pointLOG->is_open()) pointLOG->close();
            delete pointLOG;
            pointLOG = NULL;
        }
    }
}

void FTTF::setShowProgressOn() {
	showProgress = true;
}

void FTTF::setShowProgressOff() {
	showProgress = false;
}

void FTTF::clearData() {
	for (Int i = 0; i < varNum; ++i) {
		if (data[varName[i]] != NULL) delete data[varName[i]];
		if (i >= inputVarBound) { calculated[varName[i]] = false; }
	}
}

void FTTF::setParameters(std::string volRange, std::string tempRange) {
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
		*mainLOG << "FTTF model accepts new parameters:" << std::endl;//
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

void FTTF::transformDtoV() {
	PhysQvec* D = data["D"];
	PhysQvec* V = new PhysQvec(Vsize);
	Double currentV;
	for (Int i = 0; i < Vsize; ++i) {
		currentV = Mass/(Avogadro*(aVol*1e+6))/(*D)[i](gOverCmc);
		(*V)[i].setValue(currentV);
	}
	data["V"] = V;
}

void FTTF::transformCtoV() {
	PhysQvec* C = data["C"];
	PhysQvec* V = new PhysQvec(Vsize);
	Double currentV;
	for (Int i = 0; i < Vsize; ++i) {
		currentV = 1.0/(*C)[i]();
		(*V)[i].setValue(currentV);
	}
	data["V"] = V;
}

void FTTF::transformVtoD() {
	PhysQvec* V = data["V"];
	PhysQvec* D = new PhysQvec(Vsize);
	Double currentD;
	for (Int i = 0; i < Vsize; ++i) {
		currentD = Mass/(Avogadro*(aVol*1e+6))/(*V)[i]();
		(*D)[i].setValue(currentD, gOverCmc);
	}
	data["D"] = D;
}

void FTTF::transformVtoC() {
	PhysQvec* V = data["V"];
	PhysQvec* C = new PhysQvec(Vsize);
	Double currentC;
	for (Int i = 0; i < Vsize; ++i) {
		currentC = 1.0/(*V)[i]();
		(*C)[i].setValue(currentC);
	}
	data["C"] = C;
}

void FTTF::calculateData(std::string inputString) {
	Double localTime;
	/****************************LOG Section*******************/
	if (printMainLogOn || mainLogStreamIsSet) {               //
		localTime = mainTimer.getElapsedTimeInMilliSec();     //
		*mainLOG << "FTTF model accepts physical";            //
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

void FTTF::performOutput() {
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

void FTTF::printOutput(const char* filename) {
	OUT = new std::ofstream;
    OUT->open(filename, std::ios::out);
	performOutput();
    OUT->close();
    delete OUT;
}

void FTTF::printOutput(std::string filename) {
	printOutput(filename.c_str());
}

void FTTF::calculateAll() {
	Double localTime = 0;
	Double pointTime = 0;
	Int progress = 0;
	// if (showProgress) {
	// 	std::cout << "progress:          ";
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
		printer.printString(*mainLOG, "Phi(1)", 20, left);           //
		printer.printString(*mainLOG, "DPhi(0)", 20, left);          //
		printer.printString(*mainLOG, "coldPhi(1)", 20, left);       //
		printer.printString(*mainLOG, "coldDPhi(0)", 20, left);      //
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
            if (isLogVol)  filename << "log/log_FTTF(V=" << (*data["V"])[v]();         //
        	if (isLogMdns) filename << "log/log_FTTF(D=" << (*data["D"])[v](gOverCmc); //
    		if (isLogVdns) filename << "log/log_FTTF(C=" << (*data["C"])[v]();         //
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
			coldPhi.setLogStream(pointLOG);                                            //
		}                                                                              //
		/*******************************************************************************/
		coldPhi.setParameters((*(data["V"]))[v], coldT, Z);
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
            	if (isLogVol)  filename << "log/log_FTTF(V=" << (*data["V"])[v]();         //
        		if (isLogMdns) filename << "log/log_FTTF(D=" << (*data["D"])[v](gOverCmc); //
    			if (isLogVdns) filename << "log/log_FTTF(C=" << (*data["C"])[v]();         //
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
				phi.setLogStream(pointLOG);                                                //
			}                                                                              //
			/*******************************************************************************/
			phi.setParameters((*(data["V"]))[v], (*(data["T"]))[t], Z);
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
				*pointLOG << "FTTF calculation for point is finished. Elapsed time: ";         //
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

PhysQvec& FTTF::operator[] (std::string quantity) { 
	if (!calculated[quantity]) calculate(quantity);
	return *(data[quantity]);
}

/////////////////////////////////////////////////////////////////////////////////////////////////

void FTTF::calculateP(Int v, Int t) {
	Double Ptime;
	Double phi_1 = phi(1);
	if (printPointLogOn || pointLogStreamIsSet) {
		Ptime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating P:" << std::endl;
		*pointLOG << "I_{3/2}(\\phi(1)) = ";
	}
	Double currentP = FD3half(phi_1);;
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentP, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "P_{TF} = (2T)^{5/2}/(6 \\pi^2)*I_{3/2}(\\phi(1)) = ";
	}
	currentP = pow(2*(*(data["T"]))[t](), 5.0/2.0)/6.0/M_PI/M_PI*currentP;
	(*(data["P"]))[v*Tsize + t].setValue(currentP);
	calculated["P"] = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentP, precision, 20, left);
		*pointLOG << std::endl;
		Ptime = pointTimer.getElapsedTimeInMilliSec() - Ptime;
		*pointLOG << "Elapsed time = " << Ptime << " ms" << std::endl;
	}
}
void FTTF::calculatePT(Int v, Int t) {
	if (!calculated["P"])  calculateP(v, t);
	if (!calculated["PC"]) calculatePC(v);
	Double PTtime;
	if (printPointLogOn || pointLogStreamIsSet) {
		PTtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating PT:" << std::endl;
		*pointLOG << "PT = P - PC = ";
	}
	Double currentPT;
	currentPT = (*(data["P"]))[v*Tsize + t]() - (*(data["PC"]))[v]();	
	(*(data["PT"]))[v*Tsize + t].setValue(currentPT);
	calculated["PT"] = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentPT, precision, 20, left);
		*pointLOG << std::endl;
		PTtime = pointTimer.getElapsedTimeInMilliSec() - PTtime;
		*pointLOG << "Elapsed time = " << PTtime << " ms" << std::endl;
	}
}
void FTTF::calculatePC(Int v) {
	Double PCtime;
	Double coldPhi_1 = coldPhi(1);
	if (printPointLogOn || pointLogStreamIsSet) {
		PCtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating PC:" << std::endl;
		*pointLOG << "T_{C}[Hartree] = ";
		printer.printSciDouble(*pointLOG, coldT(), precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "I_{3/2}(\\phi_{C}(1)) = ";
	}
	Double currentPC = FD3half(coldPhi_1);
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentPC, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "PC_{TF} = (2T_{C})^{5/2}/(6 \\pi^2)*I_{3/2}(\\phi_{C}(1)) = ";
	}
	currentPC = pow(2*coldT(), 5.0/2.0)/6.0/M_PI/M_PI*currentPC;
	(*(data["PC"]))[v].setValue(currentPC);
	calculated["PC"] = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentPC, precision, 20, left);
		*pointLOG << std::endl;
		PCtime = pointTimer.getElapsedTimeInMilliSec() - PCtime;
		*pointLOG << "Elapsed time = " << PCtime << " ms" << std::endl;
	}
}
void FTTF::calculateE(Int v, Int t) {
	static const Double E0 = 0.76874512422;
	Double Etime;
	Double currentE;
    DoubleVec startE(0.0, rhsFTTFenergy::dim);
    Double a, b;
    Double V1, T1;
    Double xFrom = 1.0;
    Double xTo = 0.0;
    Double phi_1 = phi(1);

    if (printPointLogOn || pointLogStreamIsSet) {
		Etime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating E:" << std::endl;
		*pointLOG << "Transform parameters " << std::endl;
	}

	V1 = (*(data["V"]))[v]()*Z;
	T1 = (*(data["T"]))[t]()*pow(Z, -4.0/3.0);

    a =   pow(2.0, 7.0/6.0)
        * pow(3.0, 2.0/3.0)
        * pow(M_PI, -5.0/3.0)
        * sqrt(T1)*pow(V1, 2.0/3.0);
    b = 3.0*sqrt(2.0)*V1*pow(T1, 5.0/2.0)/M_PI/M_PI;

	if (printPointLogOn || pointLogStreamIsSet) {
		*pointLOG << "V_1 = V_Z*Z = ";
		printer.printSciDouble(*pointLOG, V1, precision, 20, left);
		*pointLOG << std::endl << "T_1 = T_Z*Z^{-4/3}";
		printer.printSciDouble(*pointLOG, T1, precision, 20, left);
		*pointLOG << std::endl << "E_0 = 0.76874512422" << std::endl;
		*pointLOG << "a = 2^{7/6}*3^{2/3}*\\pi^{-5/3}*\\T_1^{1/2}*V_1^{2/3} = ";
		printer.printSciDouble(*pointLOG, a, precision, 20, left);
		*pointLOG << std::endl << "b = 3*2^{1/2}*V_1*T_1^{5/2}*\\pi^{-2} = ";
		printer.printSciDouble(*pointLOG, b, precision, 20, left);
		*pointLOG << std::endl;
	}

    startE[0] = phi_1;
    startE[1] = phi_1;
    startE[2] = 2.0*sqrt(2.0)*V1*pow(T1, 5.0/2.0)/M_PI/M_PI
                * FD3half(phi_1) + E0;

    if (printPointLogOn || pointLogStreamIsSet) {
    	*pointLOG << "Initial vector for ODE" << std::endl;
    	*pointLOG << "y_0(1) = \\phi(1)";
    	printer.printSciDouble(*pointLOG, startE[0], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_1(1) = \\phi(1)";
    	printer.printSciDouble(*pointLOG, startE[1], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_2(1) = 2^{3/2}*V_{1}*T_{1}^{5/2}*I_{3/2}(\\phi(1)) + E0 = ";
    	printer.printSciDouble(*pointLOG, startE[2], precision, 20, left);
    	*pointLOG << std::endl;
    }

	rhsEnergy.updateParameters(a, b);
	energySolver.Integrate(rhsEnergy, startE, xFrom, xTo);
	currentE = startE[2];

	(*(data["E"]))[v*Tsize + t].setValue(currentE*pow(Z, 7.0/3.0));

	calculated["E"]  = true;

	if (printPointLogOn || pointLogStreamIsSet) {
		Int size = energyData.Count();
        printer.printString((*pointLOG), "x");
        printer.printString((*pointLOG), "phi(x)");
        printer.printString((*pointLOG), "dphi(x)");
        printer.printString((*pointLOG), "E(x)");
        (*pointLOG) << std::endl;
        for (Int i = 0; i < size; ++i) {
            printer.printSciDouble((*pointLOG), energyData.xSave[i], precision);
            printer.printSciDouble((*pointLOG), energyData.ySave[0][i], precision);
            printer.printSciDouble((*pointLOG), energyData.ySave[1][i], precision);
            printer.printSciDouble((*pointLOG), energyData.ySave[2][i], precision);
            (*pointLOG) << std::endl;
        }
        *pointLOG << "Calculated energy E_1 = E(0) = ";
        printer.printSciDouble((*pointLOG), currentE, precision);
        *pointLOG << std::endl;

		*pointLOG << "Transform to E_Z = E_1*Z^{7/3} = ";
		printer.printSciDouble((*pointLOG), currentE*pow(Z, 7.0/3.0), precision);
        *pointLOG << std::endl;

        Etime = pointTimer.getElapsedTimeInMilliSec() - Etime;
        *pointLOG << "Elapsed time = " << Etime << " ms" << std::endl;
    }

}
void FTTF::calculateET(Int v, Int t) {
	if (!calculated["E"]) calculateE(v, t);
	if (!calculated["EC"]) calculateEC(v);
	Double ETtime;
	if (printPointLogOn || pointLogStreamIsSet) {
		ETtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating ET:" << std::endl;
		*pointLOG << "ET = E - EC = ";
	}
	Double currentET;
	currentET = (*(data["E"]))[v*Tsize + t]() - (*(data["EC"]))[v]();
	(*(data["ET"]))[v*Tsize + t].setValue(currentET);
	calculated["ET"] = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentET, precision, 20, left);
		*pointLOG << std::endl;
		ETtime = pointTimer.getElapsedTimeInMilliSec() - ETtime;
		*pointLOG << "Elapsed time = " << ETtime << " ms" << std::endl;
	}
}
void FTTF::calculateEC(Int v) {
	static const Double E0 = 0.76874512422;
	Double ECtime;
	Double currentEC;
	DoubleVec startEC(0.0, rhsFTTFenergy::dim);
	Double a, b;
	Double V1, T1;
	Double xFrom = 1.0;
	Double xTo = 0.0;
	Double coldPhi_1 = coldPhi(1);

	if (printPointLogOn || pointLogStreamIsSet) {
		ECtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating EC:" << std::endl;
		*pointLOG << "Transform parameters " << std::endl;
	}
	
	T1 = coldT()*pow(Z, -4.0/3.0);
	V1 = (*(data["V"]))[v]()*Z;
	
	a =   pow(2.0, 7.0/6.0)
		* pow(3.0, 2.0/3.0)
		* pow(M_PI, -5.0/3.0)
		* sqrt(T1)*pow(V1, 2.0/3.0);
	b = 3.0*sqrt(2.0)*V1*pow(T1, 5.0/2.0)/M_PI/M_PI;

	if (printPointLogOn || pointLogStreamIsSet) {
		*pointLOG << "V_1 = V_Z*Z = ";
		printer.printSciDouble(*pointLOG, V1, precision, 20, left);
		*pointLOG << std::endl << "TC_1 = TC_Z*Z^{-4/3}";
		printer.printSciDouble(*pointLOG, T1, precision, 20, left);
		*pointLOG << std::endl << "E_0 = 0.76874512422" << std::endl;
		*pointLOG << "a = 2^{7/6}*3^{2/3}*\\pi^{-5/3}*\\TC_1^{1/2}*V_1^{2/3} = ";
		printer.printSciDouble(*pointLOG, a, precision, 20, left);
		*pointLOG << std::endl << "b = 3*2^{1/2}*V_1*TC_1^{5/2}*\\pi^{-2} = ";
		printer.printSciDouble(*pointLOG, b, precision, 20, left);
		*pointLOG << std::endl;
	}

	startEC[0] = coldPhi_1;
	startEC[1] = coldPhi_1;
	startEC[2] = 2.0*sqrt(2.0)*V1*pow(T1, 5.0/2.0)/M_PI/M_PI
                * FD3half(coldPhi(1)) + E0;

    if (printPointLogOn || pointLogStreamIsSet) {
    	*pointLOG << "Initial vector for ODE" << std::endl;
    	*pointLOG << "y_0(1) = \\phi_{C}(1)";
    	printer.printSciDouble(*pointLOG, startEC[0], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_1(1) = \\phi_{C}(1)";
    	printer.printSciDouble(*pointLOG, startEC[1], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_2(1) = 2^{3/2}*V_{1}*TC_{1}^{5/2}*I_{3/2}(\\phi(1)) + E0 = ";
    	printer.printSciDouble(*pointLOG, startEC[2], precision, 20, left);
    	*pointLOG << std::endl;
    }

	rhsColdEnergy.updateParameters(a, b);

	coldEnergySolver.Integrate(rhsColdEnergy, startEC, xFrom, xTo);
	currentEC = startEC[2];
	(*(data["EC"]))[v].setValue(currentEC*pow(Z, 7.0/3.0));

	calculated["EC"] = true;

	if (printPointLogOn || pointLogStreamIsSet) {
		Int size = coldEnergyData.Count();
        printer.printString((*pointLOG), "x");
        printer.printString((*pointLOG), "phi_{C}(x)");
        printer.printString((*pointLOG), "dphi_{C}(x)");
        printer.printString((*pointLOG), "EC(x)");
        (*pointLOG) << std::endl;
        for (Int i = 0; i < size; ++i) {
            printer.printSciDouble((*pointLOG), coldEnergyData.xSave[i], precision);
            printer.printSciDouble((*pointLOG), coldEnergyData.ySave[0][i], precision);
            printer.printSciDouble((*pointLOG), coldEnergyData.ySave[1][i], precision);
            printer.printSciDouble((*pointLOG), coldEnergyData.ySave[2][i], precision);
            (*pointLOG) << std::endl;
        }
        *pointLOG << "Calculated energy EC_1 = EC(0) = ";
        printer.printSciDouble((*pointLOG), currentEC, precision);
        *pointLOG << std::endl;

		*pointLOG << "Transform to EC_Z = EC_1*Z^{7/3} = ";
		printer.printSciDouble((*pointLOG), currentEC*pow(Z, 7.0/3.0), precision);
        *pointLOG << std::endl;

        ECtime = pointTimer.getElapsedTimeInMilliSec() - ECtime;
        *pointLOG << "Elapsed time = " << ECtime << " ms" << std::endl;
    }
}
void FTTF::calculateS(Int v, Int t) {
	Double currentS;
	Double Stime;
    DoubleVec startS(0.0, rhsFTTFentropy::dim);
    Double a, b;
    Double V1, T1;
    Double xFrom = 1.0;
    Double xTo = 0.0;
    Double phi_1 = phi(1);

    if (printPointLogOn || pointLogStreamIsSet) {
		Stime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating S:" << std::endl;
		*pointLOG << "Transform parameters " << std::endl;
	}

	V1 = (*(data["V"]))[v]()*Z;
	T1 = (*(data["T"]))[t]()*pow(Z, -4.0/3.0);

    a =   pow(2.0, 7.0/6.0)
        * pow(3.0, 2.0/3.0)
        * pow(M_PI, -5.0/3.0)
        * sqrt(T1)*pow(V1, 2.0/3.0);
    b = 7.0*sqrt(2.0*T1)*V1*T1/M_PI/M_PI;

	if (printPointLogOn || pointLogStreamIsSet) {
		*pointLOG << "V_1 = V_Z*Z = ";
		printer.printSciDouble(*pointLOG, V1, precision, 20, left);
		*pointLOG << std::endl << "T_1 = T_Z*Z^{-4/3}";
		printer.printSciDouble(*pointLOG, T1, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "a = 2^{7/6}*3^{2/3}*\\pi^{-5/3}*\\T_1^{1/2}*V_1^{2/3} = ";
		printer.printSciDouble(*pointLOG, a, precision, 20, left);
		*pointLOG << std::endl << "b = 3*2^{1/2}*V_1*T_1^{5/2}*\\pi^{-2} = ";
		printer.printSciDouble(*pointLOG, b, precision, 20, left);
		*pointLOG << std::endl;
	}
	
	startS[0] = phi_1;
    startS[1] = phi_1;
	startS[2] = 4.0*sqrt(2.0*T1)*V1*T1/M_PI/M_PI
    		    * FD3half(phi_1)
	            - phi.derivative(0);

	if (printPointLogOn || pointLogStreamIsSet) {
    	*pointLOG << "Initial vector for ODE" << std::endl;
    	*pointLOG << "y_0(1) = \\phi(1)";
    	printer.printSciDouble(*pointLOG, startS[0], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_1(1) = \\phi(1)";
    	printer.printSciDouble(*pointLOG, startS[1], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_2(1) = 4*2^{1/2}*V_{1}*T_{1}^{3/2}*I_{3/2}(\\phi(1)) - \\phi'(0) = ";
    	printer.printSciDouble(*pointLOG, startS[2], precision, 20, left);
    	*pointLOG << std::endl;
    }

	rhsEntropy.updateParameters(a, b);

	entropySolver.Integrate(rhsEntropy, startS, xFrom, xTo);
	currentS = startS[2];
	(*(data["S"]))[v*Tsize + t].setValue(currentS*Z);
    	
	calculated["S"]  = true;

	if (printPointLogOn || pointLogStreamIsSet) {
		Int size = entropyData.Count();
        printer.printString((*pointLOG), "x");
        printer.printString((*pointLOG), "phi(x)");
        printer.printString((*pointLOG), "dphi(x)");
        printer.printString((*pointLOG), "S(x)");
        (*pointLOG) << std::endl;
        for (Int i = 0; i < size; ++i) {
            printer.printSciDouble((*pointLOG), entropyData.xSave[i], precision);
            printer.printSciDouble((*pointLOG), entropyData.ySave[0][i], precision);
            printer.printSciDouble((*pointLOG), entropyData.ySave[1][i], precision);
            printer.printSciDouble((*pointLOG), entropyData.ySave[2][i], precision);
            (*pointLOG) << std::endl;
        }
        *pointLOG << "Calculated entropy S_1 = S(0) = ";
        printer.printSciDouble((*pointLOG), currentS, precision);
        *pointLOG << std::endl;

		*pointLOG << "Transform to S_Z = S_1*Z = ";
		printer.printSciDouble((*pointLOG), currentS*Z, precision);
        *pointLOG << std::endl;

        Stime = pointTimer.getElapsedTimeInMilliSec() - Stime;
        *pointLOG << "Elapsed time = " << Stime << " ms" << std::endl;
    }

}
void FTTF::calculateST(Int v, Int t) {
	if (!calculated["S"]) calculateS(v, t);
	if (!calculated["SC"]) calculateSC(v);
	Double STtime;
	if (printPointLogOn || pointLogStreamIsSet) {
		STtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating ST:" << std::endl;
		*pointLOG << "ST = S - SC = ";
	}
	Double currentST;
	currentST = (*(data["S"]))[v*Tsize + t]() - (*(data["SC"]))[v]();
	(*(data["ST"]))[v*Tsize + t].setValue(currentST);
	
	calculated["ST"] = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentST, precision, 20, left);
		*pointLOG << std::endl;
		STtime = pointTimer.getElapsedTimeInMilliSec() - STtime;
		*pointLOG << "Elapsed time = " << STtime << " ms" << std::endl;
	}
}
void FTTF::calculateSC(Int v) {
	Double SCtime;
	Double currentSC;
	DoubleVec startSC(0.0, rhsFTTFentropy::dim);
	Double a, b;
	Double V1, T1;
	Double xFrom = 1.0;
	Double xTo = 0.0;
	Double coldPhi_1 = coldPhi(1);

	if (printPointLogOn || pointLogStreamIsSet) {
		SCtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating SC:" << std::endl;
		*pointLOG << "Transform parameters " << std::endl;
	}
	
	T1 = coldT()*pow(Z, -4.0/3.0);
	V1 = (*(data["V"]))[v]()*Z;
	a =   pow(2.0, 7.0/6.0)
	    * pow(3.0, 2.0/3.0)
	    * pow(M_PI, -5.0/3.0)
	    * sqrt(T1)*pow(V1, 2.0/3.0);
	b = 7.0*sqrt(2.0*T1)*V1*T1/M_PI/M_PI;
	
	if (printPointLogOn || pointLogStreamIsSet) {
		*pointLOG << "V_1 = V_Z*Z = ";
		printer.printSciDouble(*pointLOG, V1, precision, 20, left);
		*pointLOG << std::endl << "TC_1 = TC_Z*Z^{-4/3}";
		printer.printSciDouble(*pointLOG, T1, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "a = 2^{7/6}*3^{2/3}*\\pi^{-5/3}*\\TC_1^{1/2}*V_1^{2/3} = ";
		printer.printSciDouble(*pointLOG, a, precision, 20, left);
		*pointLOG << std::endl << "b = 3*2^{1/2}*V_1*TC_1^{5/2}*\\pi^{-2} = ";
		printer.printSciDouble(*pointLOG, b, precision, 20, left);
		*pointLOG << std::endl;
	}
	
	startSC[0] = coldPhi_1;
	startSC[1] = coldPhi_1;
	startSC[2] = 4.0*sqrt(2.0*T1)*V1*T1/M_PI/M_PI
				* FD3half(coldPhi_1)
				- coldPhi.derivative(0);

	if (printPointLogOn || pointLogStreamIsSet) {
    	*pointLOG << "Initial vector for ODE" << std::endl;
    	*pointLOG << "y_0(1) = \\phi_{C}(1)";
    	printer.printSciDouble(*pointLOG, startSC[0], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_1(1) = \\phi_{C}(1)";
    	printer.printSciDouble(*pointLOG, startSC[1], precision, 20, left);
    	*pointLOG << std::endl;
    	*pointLOG << "y_2(1) = 4*2^{1/2}*V_{1}*T_{1}^{3/2}*I_{3/2}(\\phi_{C}(1)) - \\phi'_{C}(0) = ";
    	printer.printSciDouble(*pointLOG, startSC[2], precision, 20, left);
    	*pointLOG << std::endl;
    }

	rhsColdEntropy.updateParameters(a, b);

	coldEntropySolver.Integrate(rhsColdEntropy, startSC, xFrom, xTo);
	currentSC = startSC[2];
	(*(data["SC"]))[v].setValue(currentSC*Z);

	calculated["SC"] = true;

	if (printPointLogOn || pointLogStreamIsSet) {
		Int size = coldEntropyData.Count();
        printer.printString((*pointLOG), "x");
        printer.printString((*pointLOG), "phi_C(x)");
        printer.printString((*pointLOG), "dphi_C(x)");
        printer.printString((*pointLOG), "SC(x)");
        (*pointLOG) << std::endl;
        for (Int i = 0; i < size; ++i) {
            printer.printSciDouble((*pointLOG), coldEntropyData.xSave[i], precision);
            printer.printSciDouble((*pointLOG), coldEntropyData.ySave[0][i], precision);
            printer.printSciDouble((*pointLOG), coldEntropyData.ySave[1][i], precision);
            printer.printSciDouble((*pointLOG), coldEntropyData.ySave[2][i], precision);
            (*pointLOG) << std::endl;
        }
        *pointLOG << "Calculated entropy SC_1 = SC(0) = ";
        printer.printSciDouble((*pointLOG), currentSC, precision);
        *pointLOG << std::endl;

		*pointLOG << "Transform to SC_Z = SC_1*Z = ";
		printer.printSciDouble((*pointLOG), currentSC*Z, precision);
        *pointLOG << std::endl;

        SCtime = pointTimer.getElapsedTimeInMilliSec() - SCtime;
        *pointLOG << "Elapsed time = " << SCtime << " ms" << std::endl;
    }
}
void FTTF::calculateM(Int v, Int t) {
	Double Mtime;
	Double phi_1 = phi(1);
	if (printPointLogOn || pointLogStreamIsSet) {
		Mtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating M:" << std::endl;
	}
	Double currentM;
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentM, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "M_{TF} = T*\\phi(1) = ";
	}
	currentM = (*(data["T"]))[t]()*phi_1;
	(*(data["M"]))[v*Tsize + t].setValue(currentM);
	calculated["M"]  = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentM, precision, 20, left);
		*pointLOG << std::endl;
		Mtime = pointTimer.getElapsedTimeInMilliSec() - Mtime;
		*pointLOG << "Elapsed time = " << Mtime << " ms" << std::endl;
	}
}
void FTTF::calculateMT(Int v, Int t) {
	if (!calculated["M"]) calculateM(v, t);
	if (!calculated["MC"]) calculateMC(v);
	Double MTtime;
	if (printPointLogOn || pointLogStreamIsSet) {
		MTtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating MT:" << std::endl;
		*pointLOG << "MT = M - MC = ";
	}
	
	Double currentMT;
	currentMT = (*(data["M"]))[v*Tsize + t]() - (*(data["MC"]))[v]();
	(*(data["MT"]))[v*Tsize + t].setValue(currentMT);
	calculated["MT"] = true;

	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentMT, precision, 20, left);
		*pointLOG << std::endl;
		MTtime = pointTimer.getElapsedTimeInMilliSec() - MTtime;
		*pointLOG << "Elapsed time = " << MTtime << " ms" << std::endl;
	}
}
void FTTF::calculateMC(Int v) {
	Double MCtime;
	Double coldPhi_1 = coldPhi(1);
	if (printPointLogOn || pointLogStreamIsSet) {
		MCtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating MC:" << std::endl;
	}
	Double currentMC;
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentMC, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "MC_{TF} = TC*\\phi_{C}(1) = ";
	}
	currentMC = coldT()*coldPhi_1;
	(*(data["MC"]))[v].setValue(currentMC);
	calculated["MC"] = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentMC, precision, 20, left);
		*pointLOG << std::endl;
		MCtime = pointTimer.getElapsedTimeInMilliSec() - MCtime;
		*pointLOG << "Elapsed time = " << MCtime << " ms" << std::endl;
	}
}