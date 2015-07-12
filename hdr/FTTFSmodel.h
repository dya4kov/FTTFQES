#pragma once
#include "../lib/numeric/hdr/ODE/ODEsolver.h"
#include "../lib/numeric/hdr/ODE/ODEdata.h"
#include "../lib/numeric/hdr/ODE/steppers/ODEstepperPD853.h"
#include "../lib/numeric/hdr/numericTypes.h"
#include "../lib/numeric/hdr/mathUtils.h"
#include "../lib/numeric/hdr/SF/FermiDirac.h"
#include "../hdr/ThermodynamicFunction.h"
#include "../hdr/FTTFpotential.h"
#include "../hdr/FTTFelectronicStates.h"
#include "../hdr/Units.h"
#include "../hdr/stringUtils.h"
#include "../hdr/Timer.h"
#include "../hdr/Printer.h"
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>

#define DEFAULT_INPUT_FILE "FTTFSinput.dat"
#define DEFAULT_OUTPUT_FILE "FTTFSoutput.dat"

class FTTFSmodel {
public:
	/**
	* @brief A constructor of Thomas-Fermi model with shell corrections.
	*/
    FTTFSmodel(Double _Z = 1.0, Double _Mass = 1.0);
    /**
	* @brief A constructor of Thomas-Fermi model with shell corrections.
	*/
    ~FTTFSmodel();
    /**
	* @brief Set tolerance eps for the further calculations. Default is @f$ 10^{-6} @f$.
	*/
    void setTolerance(const Double eps);
    /**
	* @brief Set maximum n quantum number for calculating energy levels. Default is 6.
	*/
    void setLevelsNumber(const Int n);
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
	struct RHSEnergyInt {
		static const Int dim = 3;
		FermiDirac<Mhalf> FDmhalf;
		FermiDirac<Half> FDhalf;
		RHSEnergyInt(Double _a = 0) : a(_a) {}
		void updateParameter(const Double _a) { a = _a; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
    		dydx[0] = y[1];
		    if (x > 0) {
		        dydx[1] = a*x*FDhalf(y[0] / x);
		        dydx[2] = -x*y[0]*FDmhalf(y[0] / x);
		    }
		    else {
		        dydx[1] = 1e+10;
		        dydx[2] = 0;
		    }
		}
		private:
			Double a;
	};

	struct RHSCPInt {
		
		static const Int dim = 3;
		FermiDirac<Mhalf> FDmhalf;
		FermiDirac<Half> FDhalf;
		RHSCPInt(Double _a = 0) : a(_a) {}
		void updateParameter(const Double _a) { a = _a; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
    		dydx[0] = y[1];
		    if (x > 0) {
		        dydx[1] = a*x*FDhalf(y[0] / x);
		        dydx[2] = -x*x*FDmhalf(y[0] / x);
		    }
		    else {
		        dydx[1] = 1e+10;
		        dydx[2] = 0;
		    }
		}
	private:
		Double a;
	};
private:
	Double Z;
	Double Mass;
	Double eps;
	//input data
	std::map<std::string, Unit> unit;
	std::map<std::string, Scaling> scale;
	std::map<std::string, bool> need;
	std::map<std::string, bool> calculated;
	std::map<std::string, PhysQvec*> data;
	// output data:
	// shell corrections
	static const Int varNum = 12;
	std::string varName[varNum];
	Int inputVars;
	Int inputVarBound;
	Int Vsize;
	Int Tsize;
	//calculate single point
	void  calculate(std::string var, Int v, Int t);
	void calculateLDP(Int v, Int t); 
	void calculateLDE(Int v, Int t); 
	void calculateLDM(Int v, Int t); 

	void calculateNLDP(Int v, Int t); 
	void calculateNLDE(Int v, Int t); 
	void calculateNLDM(Int v, Int t); 

	void calculateCPInt(Int v, Int t);
	void calculateEnergyInt(Int v, Int t);

	void calculateAll(); // uses fast calculation when prepare output
	void performOutput();
	void clearData();

    // ode solvers
    RHSEnergyInt rhsEnergyInt;
    ODEsolver<ODEstepperPD853<RHSEnergyInt> > energyIntSolver;
    ODEdata energyIntData;

    RHSCPInt rhsCPInt;
    ODEsolver<ODEstepperPD853<RHSCPInt> > CPIntSolver;
    ODEdata CPIntData;

    //transform quantities
    void transformDtoV();
    void transformCtoV();
    void transformVtoD();
    void transformVtoC();

    FermiDirac<ThreeHalf> FD3half;
    FermiDirac<Half> FDhalf;
    FermiDirac<Mhalf> FDmhalf;

    FTTFpotential phi;
    FTTFelectronicStates eStates;

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

FTTFSmodel::FTTFSmodel(Double _Z, Double _Mass) : 
	Z(_Z), Mass(_Mass), eps(1e-6),
	energyIntData(-1), // save all steps
	energyIntSolver(1e-6, 0.0),
	rhsEnergyInt(0.0),
	CPIntData(-1), // save all steps
	CPIntSolver(1e-6, 0.0),
	rhsCPInt(0.0),
	eStates(6) {
		varName[0] = "T"; varName[1] = "V"; varName[2] = "D"; varName[3] = "C";  
		inputVarBound = 4; 
		varName[4] = "LDP"; varName[5] = "LDE"; varName[6] = "LDM";
		varName[7] = "NLDP"; varName[8] = "NLDE"; varName[9] = "NLDM";
		varName[10] = "Eint"; varName[11] = "CPInt";
		for (Int i = 0; i < varNum; ++i) {
			data[varName[i]] = NULL;
			unit[varName[i]] = Atomic;
			scale[varName[i]] = lin;
			need[varName[i]] = false;
			if (i >= inputVarBound) calculated[varName[i]] = false;
		}
		precision = static_cast<int>(-log10(eps));
		CPIntSolver.SetOutput(energyIntData);
		energyIntSolver.SetOutput(CPIntData);
		eStates.setEnergyRange(-4e+3, 1e+2);
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

FTTFSmodel::~FTTFSmodel() {
	clearData();
	setPrintMainLogOff();
	setPrintPointLogOff();
}

void FTTFSmodel::clearData() {
	for (Int i = 0; i < varNum; ++i) {
		if (data[varName[i]] != NULL) delete data[varName[i]];
		if (i >= inputVarBound) { calculated[varName[i]] = false; }
	}
}

void FTTFSmodel::setMainLogStream(std::ofstream* _mainLOG) {
    if (printMainLogOn) setPrintMainLogOff();
    mainLOG = _mainLOG;
    mainTimer.start();
    mainLogStreamIsSet = true;
}

void FTTFSmodel::setPrintMainLogOn() {
    if (!mainLogStreamIsSet) {
        mainLOG = new std::ofstream;
        mainTimer.start();
        std::stringstream filename;
        filename << "log/log_FTTFS_Z(";
        filename << Z << ")_M("; 
        filename << Mass << ")_";
        filename << mainTimer.getCurrentDatetime();
        filename << ".txt";
        mainLOG->open(filename.str().c_str(), std::ios::out);
        *mainLOG << "FTTFS log is started" << std::endl;
        *mainLOG << "Z = " << Z << std::endl;
        *mainLOG << "Atomic Mass = " << Mass << std::endl;
        printMainLogOn = true;
    }
}

void FTTFSmodel::setPrintMainLogOff() {
    if (!mainLogStreamIsSet) {
        if (printMainLogOn) {
        	*mainLOG << "Log stops. Final time: ";
        	*mainLOG << mainTimer.getElapsedTimeInMilliSec();
        	*mainLOG << " [ms]" << std::endl;
        	mainTimer.stop();
            printMainLogOn = false;
            if (mainLOG->is_open()) mainLOG->close();
            delete mainLOG;
            mainLOG = NULL;
        }
    }
}

void FTTFSmodel::setPointLogStream(std::ofstream* _pointLOG) {
    if (printPointLogOn) setPrintPointLogOff();
    pointLOG = _pointLOG;
    pointLogStreamIsSet = true;
}

void FTTFSmodel::setPrintPointLogOn() {
    if (!pointLogStreamIsSet) {
        pointLOG = new std::ofstream;
        printPointLogOn = true;
    }
}

void FTTFSmodel::setPrintPointLogOff() {
    if (!pointLogStreamIsSet) {
        if (printPointLogOn) {
            printPointLogOn = false;
            if (pointLOG->is_open()) pointLOG->close();
            delete pointLOG;
            pointLOG = NULL;
        }
    }
}

void FTTFSmodel::setShowProgressOn() {
	showProgress = true;
}

void FTTFSmodel::setShowProgressOff() {
	showProgress = false;
}

void FTTFSmodel::setTolerance(const Double _eps) { 
	eps = _eps;
	precision = static_cast<int>(-log10(eps));
	if (printMainLogOn || mainLogStreamIsSet) {
		*mainLOG << "FTTFS accepted new tolerance, eps = ";
		printer.printSciDouble(*mainLOG, eps, precision);
		*mainLOG << std::endl;
		phi.setLogStream(mainLOG);
	}
	phi.setTolerance(eps);
	eStates.setTolerance(eps);
	CPIntSolver.SetTolerance(0.0, eps);
	energyIntSolver.SetTolerance(0.0, eps);
	if (printMainLogOn || mainLogStreamIsSet) {
		phi.clearLogStream();
	}
}

void FTTFSmodel::calculate(std::string quantity, Int v = -1, Int t = -1) {
	if (v < 0 && t < 0) {
		calculateData(quantity);
	}
	if (v >= 0 && t >= 0) {
			 if (!quantity.compare("LDP"))   calculateLDP(v, t);
		else if (!quantity.compare("LDE"))   calculateLDE(v, t);
		else if (!quantity.compare("LDM"))   calculateLDM(v, t);
		else if (!quantity.compare("NLDP"))  calculateNLDP(v, t);
		else if (!quantity.compare("NLDE"))  calculateNLDE(v, t);
		else if (!quantity.compare("NLDM"))  calculateNLDM(v, t);
		else if (!quantity.compare("CPInt")) calculateCPInt(v, t);
		else if (!quantity.compare("Eint"))  calculateEnergyInt(v, t);
	}
}

void FTTFSmodel::setParameters(std::string volRange, std::string tempRange) {
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
	/****************************LOG Section****************************/
	if (printMainLogOn || mainLogStreamIsSet) {                        //
		*mainLOG << "FTTFS model accepts new parameters:" << std::endl;//
		if (isVolume) *mainLOG << "Volume: V[";                        //
		if (isMdns) *mainLOG << "Mass density: D[";                    //
		if (isVdns) *mainLOG << "Number density: C[";                  //
		*mainLOG << unitToString(volUnit) << ","                       //
			<< scaleToString(volScale);                                //
		*mainLOG << "](";                                              //
		*mainLOG << volRangeStart << ",";                              //
		*mainLOG << volRangeEnd   << ",";                              //
		*mainLOG << Vsize         << ")";                              //
		*mainLOG << std::endl;                                         //
		*mainLOG << "Temperature: T[";                                 //
		*mainLOG << unitToString(tempUnit) << ","                      //
			<< scaleToString(tempScale);                               //
		*mainLOG << "](";                                              //
		*mainLOG << tempRangeStart << ",";                             //
		*mainLOG << tempRangeEnd   << ",";                             //
		*mainLOG << Tsize          << ")";                             //
		*mainLOG << std::endl;                                         //
		*mainLOG << "elapsed time: " <<                                //
			mainTimer.getElapsedTimeInMilliSec() - localTime;          //
		*mainLOG << " ms" << std::endl;                                //
		*mainLOG << "Allocating memory for input parameters, ";        //
		localTime = mainTimer.getElapsedTimeInMilliSec();              //
	}                                                                  //
	/*******************************************************************/
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

void FTTFSmodel::transformDtoV() {
	PhysQvec* D = data["D"];
	PhysQvec* V = new PhysQvec(Vsize);
	Double currentV;
	for (Int i = 0; i < Vsize; ++i) {
		currentV = Mass/(Avogadro*(aVol*1e+6))/(*D)[i](gOverCmc);
		(*V)[i].setValue(currentV);
	}
	data["V"] = V;
}

void FTTFSmodel::transformCtoV() {
	PhysQvec* C = data["C"];
	PhysQvec* V = new PhysQvec(Vsize);
	Double currentV;
	for (Int i = 0; i < Vsize; ++i) {
		currentV = 1.0/(*C)[i]();
		(*V)[i].setValue(currentV);
	}
	data["V"] = V;
}

void FTTFSmodel::transformVtoD() {
	PhysQvec* V = data["V"];
	PhysQvec* D = new PhysQvec(Vsize);
	Double currentD;
	for (Int i = 0; i < Vsize; ++i) {
		currentD = Mass/(Avogadro*(aVol*1e+6))/(*V)[i]();
		(*D)[i].setValue(currentD, gOverCmc);
	}
	data["D"] = D;
}

void FTTFSmodel::transformVtoC() {
	PhysQvec* V = data["V"];
	PhysQvec* C = new PhysQvec(Vsize);
	Double currentC;
	for (Int i = 0; i < Vsize; ++i) {
		currentC = 1.0/(*V)[i]();
		(*C)[i].setValue(currentC);
	}
	data["C"] = C;
}

void FTTFSmodel::calculateData(std::string inputString) {
	Double localTime;
	/****************************LOG Section*******************/
	if (printMainLogOn || mainLogStreamIsSet) {               //
		localTime = mainTimer.getElapsedTimeInMilliSec();     //
		*mainLOG << "FTTFS model accepts physical";           //
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
    for (Int i = inputVarBound; i < varNum; ++i) {
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

void FTTFSmodel::performOutput() {
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

void FTTFSmodel::printOutput(const char* filename) {
	OUT = new std::ofstream;
    OUT->open(filename, std::ios::out);
	performOutput();
    OUT->close();
    delete OUT;
}

void FTTFSmodel::printOutput(std::string filename) {
	printOutput(filename.c_str());
}

void FTTFSmodel::calculateAll() {
	Double localTime = 0;
	Double pointTime = 0;
	Int progress = 0;
	if (showProgress) {
		std::cout << "progress:          ";
	}
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
		/*******************************************************************************/
		for (Int t = 0; t < Tsize; ++t) {
			/****************************LOG Section****************************************/
			if (printMainLogOn || mainLogStreamIsSet)  {                        //
				localTime = mainTimer.getElapsedTimeInMilliSec();                          //
				pointTime = 0;                                                             //
			}                                                                              //
			if (printPointLogOn) {                                                //
            	std::stringstream filename;                                                //
            	if (isLogVol)  filename << "log/log_FTTFS(V=" << (*data["V"])[v]();         //
        		if (isLogMdns) filename << "log/log_FTTFS(D=" << (*data["D"])[v](gOverCmc); //
    			if (isLogVdns) filename << "log/log_FTTFS(C=" << (*data["C"])[v]();         //
				filename << ", T=" << (*data["T"])[t]() << ")_";                           //
            	filename << pointTimer.getCurrentDatetime();                               //
            	filename << ".txt";                                                        //
            	if (pointLOG->is_open()) pointLOG->close();                                //
            	pointLOG->open(filename.str().c_str(), std::ios::out);                     //
			}                                                                              //
			if (printPointLogOn || pointLogStreamIsSet) {                       //
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
				eStates.setMainLogStream(pointLOG);
			}                                                                              //
			/*******************************************************************************/
			phi.setParameters((*(data["V"]))[v], (*(data["T"]))[t], Z);
			phi(1);
			eStates.setParameters((*(data["V"]))[v], (*(data["T"]))[t], Z);
			eStates.calculateEnergyLevels();
			if (showProgress) {
				std::stringstream ss;
				int progress = (100*(v*Tsize + t + 1))/(Vsize*Tsize);
				ss << progress;
				std::cout << "\r";
				std::cout << "progress: ";
				std::cout << ss.str();
				std::cout << "\%";
			}	
			for (Int i = inputVarBound; i < varNum; ++i) {                    
				if (need[varName[i]]) calculate(varName[i], v, t);
			}
			for (Int i = inputVarBound; i < varNum; ++i) {
				calculated[varName[i]] = false;
			}
			/****************************LOG Section********************************************/
			if (printPointLogOn || pointLogStreamIsSet) {                                      //
				*pointLOG << "FTTFS calculation for point is finished. Elapsed time: ";         //
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
						if (varName[i][0] != 'T' && i < inputVarBound)                                              //
							printer.printSciDouble(                                            //
								*mainLOG,                                                      //
								(*(data[varName[i]]))[v](unit[varName[i]], scale[varName[i]]), //
								precision, 20, left);                                          //
						if (varName[i][0] == 'T' && i < inputVarBound)                                              //
							printer.printSciDouble(                                            //
								*mainLOG,                                                      //
								(*(data[varName[i]]))[t](unit[varName[i]], scale[varName[i]]), //
								precision, 20, left);                                          //
						if (i >= inputVarBound)                                                 //
							printer.printSciDouble(                                            //
								*mainLOG,                                                      //
								(*(data[varName[i]]))[v*Tsize + t](unit[varName[i]],           //
								scale[varName[i]]),                                            //
								precision, 20, left);                                          //
					}                                                                          //
				}                                                                              //
				printer.printSciDouble(*mainLOG, phi(1), precision, 20, left);                 //
				printer.printSciDouble(*mainLOG, phi.derivative(0), precision, 20, left);      //
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

PhysQvec& FTTFSmodel::operator[] (std::string quantity) { 
	if (!calculated[quantity]) calculate(quantity);
	return *(data[quantity]);
}

/////////////////////////////////////////////////////////////////////////////////////////////////

void FTTFSmodel::calculateLDP(Int v, Int t) {
	if (!calculated["LDM"]) calculateLDM(v, t);
	Double LPtime;
	Double phi_1;;
    Double density;
    Double currentLP;
    Double currentLM;
    Double Temp;
	if (printPointLogOn || pointLogStreamIsSet) {
		LPtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating linear shell correction to P:" << std::endl;
		*pointLOG << "I_{1/2}(\\phi(1)) = ";
	}
	phi_1 = phi(1);
    currentLP = FDhalf(phi_1);
	Temp = (*(data["T"]))[t]();
    density = pow(Temp, 1.5)*sqrt(2.0) / M_PI / M_PI * currentLP;
    if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentLP, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "\\rho_{TF} = (2T)^{3/2}/(2 \\pi^2)*I_{1/2}(\\phi(1)) = ";
		printer.printSciDouble(*pointLOG, density, precision, 20, left);
		*pointLOG << std::endl;
	}
	// correction to chemical potential
	currentLM = (*(data["LDM"]))[v*Tsize + t]();
	// \delta P = \rho_{TF} * \delta \mu
	currentLP = density*currentLM;
    (*(data["LDP"]))[v*Tsize + t].setValue(currentLP);
	calculated["LDP"] = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		*pointLOG << "\\delta \\mu_{sh} = ";
		printer.printSciDouble(*pointLOG, currentLM, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "\\delta \\P_{sh} = ";
		printer.printSciDouble(*pointLOG, currentLP, precision, 20, left);
		*pointLOG << std::endl;
		LPtime = pointTimer.getElapsedTimeInMilliSec() - LPtime;
		*pointLOG << "Elapsed time = " << LPtime << " ms" << std::endl;
	}
}

void FTTFSmodel::calculateLDE(Int v, Int t) {
	if (!calculated["LDM"]) calculateLDM(v, t);
	if (!calculated["Eint"]) calculateEnergyInt(v, t);
	Double LEtime;
	Double currentLE;
	Double currentLM;
	Double currentEint;
	if (printPointLogOn || pointLogStreamIsSet) {
		LEtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating linear shell correction to E:" << std::endl;
	}
	currentLM = (*(data["LDM"]))[v*Tsize + t]();
	currentEint = (*(data["Eint"]))[v*Tsize + t]();
	currentLE = (1.5*Z - currentEint)*currentLM;
	(*(data["LDE"]))[v*Tsize + t].setValue(currentLE);
	calculated["LDE"]  = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		*pointLOG << "\\delta \\mu_{sh} = ";
		printer.printSciDouble(*pointLOG, currentLM, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "\\int = ";
		printer.printSciDouble(*pointLOG, currentEint, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "\\delta \\E_{sh} = (3/2*Z - \\int)*\\delta \\mu_{sh} = ";
		printer.printSciDouble(*pointLOG, currentLE, precision, 20, left);
		*pointLOG << std::endl;
		LEtime = pointTimer.getElapsedTimeInMilliSec() - LEtime;
		*pointLOG << "Elapsed time = " << LEtime << " ms" << std::endl;
	}
}

void FTTFSmodel::calculateLDM(Int v, Int t) {
    Double currentDN;
    currentDN = eStates.DNlow();
	if (!calculated["CPInt"]) calculateCPInt(v, t);
    Double currentLM;
    Double currentCPInt;
    Double LMtime;
    if (printPointLogOn || pointLogStreamIsSet) {
		LMtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating linear shell correction to \\mu:" << std::endl;
	}
    currentCPInt = (*(data["CPInt"]))[v*Tsize + t]();
    currentLM = -currentDN/currentCPInt;
    (*(data["LDM"]))[v*Tsize + t].setValue(currentLM);
	calculated["LDM"]  = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		*pointLOG << "\\int = ";
		printer.printSciDouble(*pointLOG, currentCPInt, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "\\delta \\mu_{sh} = ";
		printer.printSciDouble(*pointLOG, currentLM, precision, 20, left);
		*pointLOG << std::endl;
		LMtime = pointTimer.getElapsedTimeInMilliSec() - LMtime;
		*pointLOG << "Elapsed time = " << LMtime << " ms" << std::endl;
	}	
}

void FTTFSmodel::calculateCPInt(Int v, Int t) {
	Double CPIntTime;
	if (printPointLogOn || pointLogStreamIsSet) {
		CPIntTime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating CPInt:" << std::endl;
	}
	Double currentCPInt;
	DoubleVec startCPInt(0.0, RHSCPInt::dim);
	Double xFrom = 1.0;
	Double xTo = 0.0;
	Double Vol = (*(data["V"]))[v]();
	Double Temp = (*(data["T"]))[t]();

	startCPInt[0] = phi(1);
	startCPInt[1] = startCPInt[0];
	startCPInt[2] = 0.0;

	CPIntSolver.Integrate(rhsCPInt, startCPInt, xFrom, xTo);
	currentCPInt = 3*Vol/ 4.0 / M_PI * 2 / M_PI*sqrt(2 * Temp)*startCPInt[2];
	(*(data["CPInt"]))[v*Tsize + t].setValue(currentCPInt);

	calculated["CPInt"]  = true;

	if (printPointLogOn || pointLogStreamIsSet) {
		Int size = CPIntData.Count();
        printer.printString((*pointLOG), "x");
        printer.printString((*pointLOG), "phi(x)");
        printer.printString((*pointLOG), "dphi(x)");
        printer.printString((*pointLOG), "CPInt(x)");
        (*pointLOG) << std::endl;
        for (Int i = 0; i < size; ++i) {
            printer.printSciDouble((*pointLOG), CPIntData.xSave[i], precision);
            printer.printSciDouble((*pointLOG), CPIntData.ySave[0][i], precision);
            printer.printSciDouble((*pointLOG), CPIntData.ySave[1][i], precision);
            printer.printSciDouble((*pointLOG), CPIntData.ySave[2][i], precision);
            (*pointLOG) << std::endl;
        }
        *pointLOG << "Calculated CPInt = ";
        printer.printSciDouble((*pointLOG), currentCPInt, precision);
        *pointLOG << std::endl;

        CPIntTime = pointTimer.getElapsedTimeInMilliSec() - CPIntTime;
        *pointLOG << "Elapsed time = " << CPIntTime << " ms" << std::endl;
    }
}

void FTTFSmodel::calculateEnergyInt(Int v, Int t) {
	Double EnergyIntTime;
	if (printPointLogOn || pointLogStreamIsSet) {
		EnergyIntTime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating EnergyInt:" << std::endl;
	}
	Double currentEnergyInt;
	DoubleVec startEnergyInt(0.0, RHSEnergyInt::dim);
	Double xFrom = 1.0;
	Double xTo = 0.0;
	Double Vol = (*(data["V"]))[v]();
	Double Temp = (*(data["T"]))[t]();
	
	startEnergyInt[0] = phi(1);
	startEnergyInt[1] = startEnergyInt[0];
	startEnergyInt[2] = 0.0;
	
	energyIntSolver.Integrate(rhsEnergyInt, startEnergyInt, xFrom, xTo);

	currentEnergyInt = 3*Vol / 4.0 / M_PI *pow(2.0 * Temp, 1.5) / M_PI*startEnergyInt[2];
	(*(data["Eint"]))[v*Tsize + t].setValue(currentEnergyInt);
	
	calculated["Eint"]  = true;

	if (printPointLogOn || pointLogStreamIsSet) {
		Int size = energyIntData.Count();
		*pointLOG << "Number of integration steps: " << size << std::endl;
        printer.printString((*pointLOG), "x");
        printer.printString((*pointLOG), "phi(x)");
        printer.printString((*pointLOG), "dphi(x)");
        printer.printString((*pointLOG), "EnergyInt(x)");
        (*pointLOG) << std::endl;
        for (Int i = 0; i < size; ++i) {
            printer.printSciDouble((*pointLOG), energyIntData.xSave[i], precision);
            printer.printSciDouble((*pointLOG), energyIntData.ySave[0][i], precision);
            printer.printSciDouble((*pointLOG), energyIntData.ySave[1][i], precision);
            printer.printSciDouble((*pointLOG), energyIntData.ySave[2][i], precision);
            (*pointLOG) << std::endl;
        }
        *pointLOG << "Calculated EnergyInt = ";
        printer.printSciDouble((*pointLOG), currentEnergyInt, precision);
        *pointLOG << std::endl;

        EnergyIntTime = pointTimer.getElapsedTimeInMilliSec() - EnergyIntTime;
        *pointLOG << "Elapsed time = " << EnergyIntTime << " ms" << std::endl;
    }
}

void FTTFSmodel::calculateNLDP(Int v, Int t) {
    if (!calculated["NLDM"]) calculateNLDM(v, t);
	Double NLPtime;
	Double phi_1;;
    Double density;
    Double currentNLP;
    Double currentNLM;
	if (printPointLogOn || pointLogStreamIsSet) {
		NLPtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating nonlinear shell correction to P:" << std::endl;
		*pointLOG << "I_{1/2}(\\phi(1)) = ";
	}
	Double Temp = (*(data["T"]))[t]();
	phi_1 = phi(1);
    currentNLP = FDhalf(phi_1);
    density = pow(Temp, 1.5)*sqrt(2.0) / M_PI / M_PI * currentNLP;
    if (printPointLogOn || pointLogStreamIsSet) {
		printer.printSciDouble(*pointLOG, currentNLP, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "\\rho_{TF} = (2T)^{3/2}/(2 \\pi^2)*I_{1/2}(\\phi(1)) = ";
		printer.printSciDouble(*pointLOG, density, precision, 20, left);
		*pointLOG << std::endl;
	}
	// correction to chemical potential
	currentNLM = (*(data["NLDM"]))[v*Tsize + t]();
	// \delta P = \rho_{TF} * \delta \mu
	currentNLP = density*currentNLM;
    (*(data["NLDP"]))[v*Tsize + t].setValue(currentNLP);
	calculated["NLDP"] = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		*pointLOG << "\\delta \\mu_{sh} = ";
		printer.printSciDouble(*pointLOG, currentNLM, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "\\delta \\P_{sh} = ";
		printer.printSciDouble(*pointLOG, currentNLP, precision, 20, left);
		*pointLOG << std::endl;
		NLPtime = pointTimer.getElapsedTimeInMilliSec() - NLPtime;
		*pointLOG << "Elapsed time = " << NLPtime << " ms" << std::endl;
	}
}

void FTTFSmodel::calculateNLDE(Int v, Int t) {
	if (!calculated["NLDM"]) calculateNLDM(v, t);
	if (!calculated["Eint"]) calculateEnergyInt(v, t);
	Double NLEtime;
	Double currentNLE;
	Double currentNLM;
	Double currentEint;
	if (printPointLogOn || pointLogStreamIsSet) {
		NLEtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating nonlinear shell correction to E:" << std::endl;
	}
	currentNLM = (*(data["NLDM"]))[v*Tsize + t]();
	currentEint = (*(data["Eint"]))[v*Tsize + t]();
	currentNLE = (1.5*Z - currentEint)*currentNLM;
	(*(data["NLDE"]))[v*Tsize + t].setValue(currentNLE);
	calculated["NLDE"]  = true;
	if (printPointLogOn || pointLogStreamIsSet) {
		*pointLOG << "\\delta \\mu_{sh} = ";
		printer.printSciDouble(*pointLOG, currentNLM, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "\\int = ";
		printer.printSciDouble(*pointLOG, currentEint, precision, 20, left);
		*pointLOG << std::endl;
		*pointLOG << "\\delta \\E_{sh} = (3/2*Z - \\int)*\\delta \\mu_{sh} = ";
		printer.printSciDouble(*pointLOG, currentNLE, precision, 20, left);
		*pointLOG << std::endl;
		NLEtime = pointTimer.getElapsedTimeInMilliSec() - NLEtime;
		*pointLOG << "Elapsed time = " << NLEtime << " ms" << std::endl;
	}
}

void FTTFSmodel::calculateNLDM(Int v, Int t) {
	Double currentN;
	Double MTF = phi(1);
	Double exactN = Z;
	Double leftM = -fabs(MTF * 1.1);
	Double rightM = fabs(MTF * 1.2);
	Double centerM = 0;
	Double currentNLM;
	Double NLMtime;
    Double Temp = (*(data["T"]))[t]();
	

	Int nStep = 0;
	Double err;

	if (printPointLogOn || pointLogStreamIsSet) {
		NLMtime = pointTimer.getElapsedTimeInMilliSec();
		*pointLOG << "Calculating nonlinear shell correction to M:" << std::endl;
		printer.printString(*pointLOG, "nStep", 20, left);
		printer.printString(*pointLOG, "leftMsh/T", 20, left);
		printer.printString(*pointLOG, "rightMsh/T", 20, left);
		printer.printString(*pointLOG, "centerMsh/T", 20, left);
		printer.printString(*pointLOG, "Nstates", 20, left);
		printer.printString(*pointLOG, "error", 20, left);
		printer.printString(*pointLOG, "time[ms]", 20, left);
		*pointLOG << std::endl;
	}
	err = abs(leftM - rightM) / abs(leftM + rightM);
	centerM = 0.5*(rightM + leftM);
	// disable eStates log writing for good table appearance
	eStates.clearMainLogStream();
	while (err > eps*10) {
		eStates.setPhiShift(centerM);
		currentN = eStates.N();
		// std::cout << "DM = " << 0.5*(leftM + rightM) << std::endl;
		// std::cout << "currentN = " << currentN << std::endl;
		if (printPointLogOn || pointLogStreamIsSet) {
			printer.printInt(*pointLOG, nStep, 20, left);
			printer.printSciDouble(*pointLOG, leftM, precision, 20, left);
			printer.printSciDouble(*pointLOG, rightM, precision, 20, left);
			printer.printSciDouble(*pointLOG, centerM, precision, 20, left);
			printer.printSciDouble(*pointLOG, currentN, precision, 20, left);
		}
		if (abs(currentN - exactN) < 1e-20) {
			leftM = rightM = 0;
			break;
		}
		if (currentN - exactN > 0) {
			rightM -= (rightM - leftM)*0.5;
		}
		else {
			leftM += (rightM - leftM)*0.5;
		}
		err = abs(leftM - rightM) / abs(leftM + rightM);
		if (printPointLogOn || pointLogStreamIsSet) {
			NLMtime = pointTimer.getElapsedTimeInMilliSec() - NLMtime;
			printer.printSciDouble(*pointLOG, err, precision, 20, left);
			printer.printSciDouble(*pointLOG, NLMtime, precision, 20, left);
			*pointLOG << std::endl;
		}
		centerM = 0.5*(rightM + leftM);
	}

	currentNLM = centerM*Temp;

	(*(data["NLDM"]))[v*Tsize + t].setValue(currentNLM);
	calculated["NLDM"] = true;
	// enable eStates main log stream
	eStates.setMainLogStream(pointLOG);
	if (printPointLogOn || pointLogStreamIsSet) {
		*pointLOG << "Finally calculated nonlinear Msh: " << std::endl;
		printer.printSciDouble(*pointLOG, currentNLM, precision, 20, left);
		*pointLOG << std::endl;
	}
}