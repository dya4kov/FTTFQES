#pragma once
#include "../lib/numeric/hdr/ODE/ODEsolver.h"
#include "../lib/numeric/hdr/ODE/ODEdata.h"
#include "../lib/numeric/hdr/ODE/steppers/ODEstepperPD853.h"
#include "../lib/numeric/hdr/numericTypes.h"
#include "../lib/numeric/hdr/mathUtils.h"
#include "../lib/numeric/hdr/SF/FermiDirac.h"
#include "../hdr/ThermodynamicFunction.h"
#include "../hdr/FTTFpotential.h"
#include "../hdr/Units.h"
#include "../hdr/stringUtils.h"
#include <fstream>
#include <iostream>
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
class FTTFmodel {
public:
	/**
	* @brief A constructor of Thomas-Fermi model.
	*/
    FTTFmodel(Double _Z = 1.0, Double _Mass = 1.0);
    /**
	* @brief Set tolerance eps for the further calculations. Default is @f$ 10^{-6} @f$.
	*/
    void setTolerance(const Double eps);
    /**
	* @brief Set temperature and volume vectors for the further calculations.
	*/
    void setParameters(const VolVec &V, const TempVec &T);
    /**
	* @brief Set temperature and density vectors for the further calculations.
	*/
    void setParameters(const DensVec &D, const TempVec &T);
    /**
	* @brief Set temperature and concentration vectors for the further calculations.
	*/
    void setParameters(const ConVec &C, const TempVec &T);
    /**
	* @brief Calculate output data as defined in inputString.
	* @details Example of string: P[GPa] PT[GPa] E[eV] S M[eV] PV/T Q ...
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
	* - QEoS[Units, Scaling] - caloric equation of state @f$ E - 3/2PV @f$, units: eV, Atomic;
	* - TEoS - thermal equation of state;
	*/
    void prepareOutput(std::string inputString);
    /**
	* @brief Print output data in file as defined in inputString.
	*/
    void printOutput(const char* filename);

	PressVec&   getP(); 
	PressVec&   getPT();
	PressVec&   getPC();
	EnergyVec&  getE();
	EnergyVec&  getET();
	EnergyVec&  getEC();
	EntropyVec& getS();
	EntropyVec& getST();
	EntropyVec& getSC();
	ChemPotVec& getM();
	ChemPotVec& getMT();
	ChemPotVec& getMC();
	EnergyVec&  getQEoS();
	DoubleVec&  getTEoS();
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
	//input data
	TempVec* T; Unit unitT; Scaling scaleT; bool needT;
	VolVec*  V; Unit unitV; Scaling scaleV; bool needV;
	DensVec* D; Unit unitD; Scaling scaleD; bool needD;
	ConVec*  C; Unit unitC; Scaling scaleC; bool needC;

	// output data
	PressVec*     P; Unit  unitP;   Scaling  scaleP;   bool    calculatedP; bool    needP;
	PressVec*    PT; Unit unitPT;   Scaling scalePT;   bool   calculatedPT; bool   needPT; 
	PressVec*    PC; Unit unitPC;   Scaling scalePC;   bool   calculatedPC; bool   needPC; 
	EnergyVec*    E; Unit  unitE;   Scaling  scaleE;   bool    calculatedE; bool    needE; 
	EnergyVec*   ET; Unit unitET;   Scaling scaleET;   bool   calculatedET; bool   needET; 
	EnergyVec*   EC; Unit unitEC;   Scaling scaleEC;   bool   calculatedEC; bool   needEC; 
	EntropyVec*   S; Unit  unitS;   Scaling  scaleS;   bool    calculatedS; bool    needS; 
	EntropyVec*  ST; Unit unitST;   Scaling scaleST;   bool   calculatedST; bool   needST; 
	EntropyVec*  SC; Unit unitSC;   Scaling scaleSC;   bool   calculatedSC; bool   needSC; 
	ChemPotVec*   M; Unit  unitM;   Scaling  scaleM;   bool    calculatedM; bool    needM; 
	ChemPotVec*  MT; Unit unitMT;   Scaling scaleMT;   bool   calculatedMT; bool   needMT; 
	ChemPotVec*  MC; Unit unitMC;   Scaling scaleMC;   bool   calculatedMC; bool   needMC;
	EnergyVec* QEoS; Unit unitQEoS; Scaling scaleQEoS; bool calculatedQEoS; bool needQEoS; 
	DoubleVec* TEoS;                                   bool calculatedTEoS; bool needTEoS;

	// calculate full vector |      calculate single point
	void        calculateP();  void    calculateP(Int v, Int t); 
	void       calculatePT();  void   calculatePT(Int v, Int t); 
	void       calculatePC();  void   calculatePC(Int v); 
	void        calculateE();  void    calculateE(Int v, Int t); 
	void       calculateET();  void   calculateET(Int v, Int t); 
	void       calculateEC();  void   calculateEC(Int v); 
	void        calculateS();  void    calculateS(Int v, Int t); 
	void       calculateST();  void   calculateST(Int v, Int t); 
	void       calculateSC();  void   calculateSC(Int v); 
	void        calculateM();  void    calculateM(Int v, Int t); 
	void       calculateMT();  void   calculateMT(Int v, Int t); 
	void       calculateMC();  void   calculateMC(Int v); 
	void     calculateQEoS();  void calculateQEoS(Int v, Int t); 
	void     calculateTEoS();  void calculateTEoS(Int v, Int t); 

	void calculateAll(); // uses fast calculation when prepare output
	void printString(std::ofstream &out, std::string data);
	void printValue(std::ofstream &out, Double value, Int digitsToPrint);

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
    FTTFpotential phi;

    FTTFpotential coldPhi;
    Temperature coldT;
};

FTTFmodel::FTTFmodel(Double _Z, Double _Mass) : 
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
		T = NULL;    unitT = Atomic;    scaleT = lin;    needT = false;                         
		V = NULL;    unitV = Atomic;    scaleV = lin;    needV = false;                         
		D = NULL;    unitD = gOverCmc;  scaleD = lin;    needD = false;                         
		C = NULL;    unitC = Atomic;    scaleC = lin;    needC = false;

		P =  NULL;    unitP = Atomic;    scaleP = lin;    needP = false;    calculatedP = false;    
		PT = NULL;   unitPT = Atomic;   scalePT = lin;   needPT = false;   calculatedPT = false;   
		PC = NULL;   unitPC = Atomic;   scalePC = lin;   needPC = false;   calculatedPC = false;   
		E =  NULL;    unitE = Atomic;    scaleE = lin;    needE = false;    calculatedE = false;    
		ET = NULL;   unitET = Atomic;   scaleET = lin;   needET = false;   calculatedET = false;   
		EC = NULL;   unitEC = Atomic;   scaleEC = lin;   needEC = false;   calculatedEC = false;   
		S =  NULL;    unitS = Atomic;    scaleS = lin;    needS = false;    calculatedS = false;    
		ST = NULL;   unitST = Atomic;   scaleST = lin;   needST = false;   calculatedST = false;   
		SC = NULL;   unitSC = Atomic;   scaleSC = lin;   needSC = false;   calculatedSC = false;   
		M =  NULL;    unitM = Atomic;    scaleM = lin;    needM = false;    calculatedM = false;    
		MT = NULL;   unitMT = Atomic;   scaleMT = lin;   needMT = false;   calculatedMT = false;   
		MC = NULL;   unitMC = Atomic;   scaleMC = lin;   needMC = false;   calculatedMC = false;   

		QEoS = NULL; unitQEoS = Atomic; scaleQEoS = lin; needQEoS = false; calculatedQEoS = false; 
		TEoS = NULL;									 needTEoS = false;
	}

void FTTFmodel::setTolerance(const Double _eps) { 
	eps = _eps; 
	phi.setTolerance(eps);
	coldPhi.setTolerance(eps);
}

void FTTFmodel::setParameters(const VolVec &_V, const TempVec &_T) {
	if (V != NULL) delete V;
	if (T != NULL) delete T;
	V = new VolVec(_V);
	T = new TempVec(_T);
	transformVtoC();
	transformVtoD();
	needV = true;
	needT = true;
	calculatedP  = false; if (P != NULL) delete P;
	calculatedPT = false; if (PT != NULL) delete PT;
	calculatedPC = false; if (PC != NULL) delete PC;
	calculatedE  = false; if (E != NULL) delete E;
	calculatedET = false; if (ET != NULL) delete ET;
	calculatedEC = false; if (EC != NULL) delete EC;
	calculatedS  = false; if (S != NULL) delete S;
	calculatedST = false; if (ST != NULL) delete ST;
	calculatedSC = false; if (SC != NULL) delete SC;
	calculatedM  = false; if (M != NULL) delete M;
	calculatedMT = false; if (MT != NULL) delete MT;
	calculatedMC = false; if (MC != NULL) delete MC;
	calculatedQEoS = false; if (QEoS != NULL) delete QEoS;
	calculatedTEoS = false; if (TEoS != NULL) delete TEoS;
}

void FTTFmodel::setParameters(const DensVec &_D, const TempVec &_T) {
	if (D != NULL) delete D;
	if (T != NULL) delete T;
	D = new DensVec(_D);
	T = new TempVec(_T);
	needD = true;
	needT = true;
	transformDtoV();
	transformVtoC();
	calculatedP  = false; if (P != NULL) delete P;
	calculatedPT = false; if (PT != NULL) delete PT;
	calculatedPC = false; if (PC != NULL) delete PC;
	calculatedE  = false; if (E != NULL) delete E;
	calculatedET = false; if (ET != NULL) delete ET;
	calculatedEC = false; if (EC != NULL) delete EC;
	calculatedS  = false; if (S != NULL) delete S;
	calculatedST = false; if (ST != NULL) delete ST;
	calculatedSC = false; if (SC != NULL) delete SC;
	calculatedM  = false; if (M != NULL) delete M;
	calculatedMT = false; if (MT != NULL) delete MT;
	calculatedMC = false; if (MC != NULL) delete MC;
	calculatedQEoS = false; if (QEoS != NULL) delete QEoS;
	calculatedTEoS = false; if (TEoS != NULL) delete TEoS;
}

void FTTFmodel::setParameters(const ConVec &_C, const TempVec &_T) {
	if (C != NULL) delete C;
	if (T != NULL) delete T;
	C = new ConVec(_C);
	T = new TempVec(_T);
	needC = true;
	needT = true;
	transformCtoV();
	transformVtoD();
	calculatedP  = false; if (P != NULL) delete P;
	calculatedPT = false; if (PT != NULL) delete PT;
	calculatedPC = false; if (PC != NULL) delete PC;
	calculatedE  = false; if (E != NULL) delete E;
	calculatedET = false; if (ET != NULL) delete ET;
	calculatedEC = false; if (EC != NULL) delete EC;
	calculatedS  = false; if (S != NULL) delete S;
	calculatedST = false; if (ST != NULL) delete ST;
	calculatedSC = false; if (SC != NULL) delete SC;
	calculatedM  = false; if (M != NULL) delete M;
	calculatedMT = false; if (MT != NULL) delete MT;
	calculatedMC = false; if (MC != NULL) delete MC;
	calculatedQEoS = false; if (QEoS != NULL) delete QEoS;
	calculatedTEoS = false; if (TEoS != NULL) delete TEoS;
}

void FTTFmodel::transformDtoV() {
	Int n = D->size();
	V = new VolVec(n);
	Double currentV;
	for (Int i = 0; i < n; ++i) {
		currentV = Mass/(Avogadro*1.4818e-25)/(*D)[i](gOverCmc);
		(*V)[i].setValue(currentV);
	}
}

void FTTFmodel::transformCtoV() {
	Int n = C->size();
	V = new VolVec(n);
	Double currentV;
	for (Int i = 0; i < n; ++i) {
		currentV = 1.0/(*C)[i]();
		(*V)[i].setValue(currentV);
	}
}

void FTTFmodel::transformVtoD() {
	Int n = V->size();
	D = new DensVec(n);
	Double currentD;
	for (Int i = 0; i < n; ++i) {
		currentD = Mass/(Avogadro*1.4818e-25)/(*V)[i]();
		(*D)[i].setValue(currentD, gOverCmc);
	}
}

void FTTFmodel::transformVtoC() {
	Int n = V->size();
	C = new ConVec(n);
	Double currentC;
	for (Int i = 0; i < n; ++i) {
		currentC = 1.0/(*V)[i]();
		(*C)[i].setValue(currentC);
	}
}

void FTTFmodel::prepareOutput(std::string inputString) {
	const char startUnit = '[';
    const char endUnit = ']';
    const char delim = ' ';
    const char unitScaleDelim = ',';
    // split inputString into tokens
    std::vector<std::string> tokens = split(inputString, delim);
    // analyze tokens
    std::vector<std::string>::iterator itoken;
    for (itoken = tokens.begin(); itoken != tokens.end(); ++itoken) {
    	std::string quantity, unit, scale;
    	quantity = *itoken;
    	if (!isalpha(*(quantity.end() - 1))) {
			unit = split(quantity, startUnit)[1];
			unit = std::string(unit.begin(), unit.end() - 1);
			quantity = split(quantity, startUnit)[0];
			if (split(unit, unitScaleDelim).size() > 1){
				scale = split(unit, unitScaleDelim)[1];
				unit = split(unit, unitScaleDelim)[0];
			}
    	}
    	if (!quantity.compare("V"))    { needV    = true; if (!unit.empty()) unitV = stringToUnit(unit);    if (!scale.empty()) scaleV = stringToScale(scale); }
    	if (!quantity.compare("D"))    { needD    = true; if (!unit.empty()) unitD = stringToUnit(unit);    if (!scale.empty()) scaleD = stringToScale(scale); }
 		if (!quantity.compare("C"))    { needC    = true; if (!unit.empty()) unitC = stringToUnit(unit);    if (!scale.empty()) scaleC = stringToScale(scale); }
    	if (!quantity.compare("T"))    { needT    = true; if (!unit.empty()) unitT = stringToUnit(unit);    if (!scale.empty()) scaleT = stringToScale(scale); }
    	if (!quantity.compare("P"))    { needP    = true; if (!unit.empty()) unitP = stringToUnit(unit);    if (!scale.empty()) scaleP = stringToScale(scale); }
		if (!quantity.compare("PT"))   { needPT   = true; if (!unit.empty()) unitPT = stringToUnit(unit);   if (!scale.empty()) scalePT = stringToScale(scale); }
		if (!quantity.compare("PC"))   { needPC   = true; if (!unit.empty()) unitPC = stringToUnit(unit);   if (!scale.empty()) scalePC = stringToScale(scale); }
		if (!quantity.compare("E"))    { needE    = true; if (!unit.empty()) unitE = stringToUnit(unit);    if (!scale.empty()) scaleE = stringToScale(scale); }
		if (!quantity.compare("ET"))   { needET   = true; if (!unit.empty()) unitET = stringToUnit(unit);   if (!scale.empty()) scaleET = stringToScale(scale); }
		if (!quantity.compare("EC"))   { needEC   = true; if (!unit.empty()) unitEC = stringToUnit(unit);   if (!scale.empty()) scaleEC = stringToScale(scale); }
		if (!quantity.compare("S"))    { needS    = true; if (!unit.empty()) unitS = stringToUnit(unit);    if (!scale.empty()) scaleS = stringToScale(scale); }
		if (!quantity.compare("ST"))   { needST   = true; if (!unit.empty()) unitST = stringToUnit(unit);   if (!scale.empty()) scaleST = stringToScale(scale); }
		if (!quantity.compare("SC"))   { needSC   = true; if (!unit.empty()) unitSC = stringToUnit(unit);   if (!scale.empty()) scaleSC = stringToScale(scale); }
		if (!quantity.compare("M"))    { needM    = true; if (!unit.empty()) unitM = stringToUnit(unit);    if (!scale.empty()) scaleM = stringToScale(scale); }
		if (!quantity.compare("MT"))   { needMT   = true; if (!unit.empty()) unitMT = stringToUnit(unit);   if (!scale.empty()) scaleMT = stringToScale(scale); }
		if (!quantity.compare("MC"))   { needMC   = true; if (!unit.empty()) unitMC = stringToUnit(unit);   if (!scale.empty()) scaleMC = stringToScale(scale); }
		if (!quantity.compare("QEoS")) { needQEoS = true; if (!unit.empty()) unitQEoS = stringToUnit(unit); if (!scale.empty()) scaleQEoS = stringToScale(scale); }
		if (!quantity.compare("TEoS")) { needTEoS = true; }
    }
    calculateAll();
}

void FTTFmodel::printString(std::ofstream &out, std::string data) {
	out.setf(std::ios::left);
    out.width(20);
    out.fill(' ');
    out << data;
}

void FTTFmodel::printValue(std::ofstream &out, Double value, Int digitsToPrint) {
	out.setf(std::ios::scientific | std::ios::left | std::ios::showpos);
    out.width(20);
    out.precision(digitsToPrint);
    out.fill(' ');
    out << value;
}

void FTTFmodel::printOutput(const char* filename) {
	std::ofstream out;
    Int digitsToPrint = static_cast<int>(-log10(eps));
    Int Vsize = V->size();
    Int Tsize = T->size();
    out.open(filename, std::ios::out);
    // print first line:
    std::string s;
    if (needV)    { s = "V[";    s += unitToString(unitV);    s += ","; s += scaleToString(scaleV);    s += "]"; printString(out, s); }
	if (needD)    { s = "D[";    s += unitToString(unitD);    s += ","; s += scaleToString(scaleD);    s += "]"; printString(out, s); }
	if (needC)    { s = "C[";    s += unitToString(unitC);    s += ","; s += scaleToString(scaleC);    s += "]"; printString(out, s); }
	if (needT)    { s = "T[";    s += unitToString(unitT);    s += ","; s += scaleToString(scaleT);    s += "]"; printString(out, s); }
	if (needP)    { s = "P[";    s += unitToString(unitP);    s += ","; s += scaleToString(scaleP);    s += "]"; printString(out, s); }
	if (needPT)   { s = "PT[";   s += unitToString(unitPT);   s += ","; s += scaleToString(scalePT);   s += "]"; printString(out, s); }
	if (needPC)   { s = "PC[";   s += unitToString(unitPC);   s += ","; s += scaleToString(scalePC);   s += "]"; printString(out, s); }
	if (needE)    { s = "E[";    s += unitToString(unitE);    s += ","; s += scaleToString(scaleE);    s += "]"; printString(out, s); }
	if (needET)   { s = "ET[";   s += unitToString(unitET);   s += ","; s += scaleToString(scaleET);   s += "]"; printString(out, s); }
	if (needEC)   { s = "EC[";   s += unitToString(unitEC);   s += ","; s += scaleToString(scaleEC);   s += "]"; printString(out, s); }
	if (needS)    { s = "S[";    s += unitToString(unitS);    s += ","; s += scaleToString(scaleS);    s += "]"; printString(out, s); }
	if (needST)   { s = "ST[";   s += unitToString(unitST);   s += ","; s += scaleToString(scaleST);   s += "]"; printString(out, s); }
	if (needSC)   { s = "SC[";   s += unitToString(unitSC);   s += ","; s += scaleToString(scaleSC);   s += "]"; printString(out, s); }
	if (needM)    { s = "M[";    s += unitToString(unitM);    s += ","; s += scaleToString(scaleM);    s += "]"; printString(out, s); }
	if (needMT)   { s = "MT[";   s += unitToString(unitMT);   s += ","; s += scaleToString(scaleMT);   s += "]"; printString(out, s); }
	if (needMC)   { s = "MC[";   s += unitToString(unitMC);   s += ","; s += scaleToString(scaleMC);   s += "]"; printString(out, s); }
	if (needQEoS) { s = "QEoS["; s += unitToString(unitQEoS); s += ","; s += scaleToString(scaleQEoS); s += "]"; printString(out, s); }
	if (needTEoS) { s = "TEoS"; printString(out, s); }
	out << std::endl;
	// print data
    for (Int v = 0; v < Vsize; ++v) {
	    for (Int t = 0; t < Tsize; ++t) {
	    	if (needV)    printValue(out, (*V)[v](unitV, scaleV), digitsToPrint);
	    	if (needD)    printValue(out, (*D)[v](unitD, scaleD), digitsToPrint);
	    	if (needC)    printValue(out, (*C)[v](unitC, scaleC), digitsToPrint);
	    	if (needT)    printValue(out, (*T)[t](unitT, scaleT), digitsToPrint);
			if (needP)    printValue(out, (*P)[v*Tsize + t](unitP, scaleP), digitsToPrint);
			if (needPT)   printValue(out, (*PT)[v*Tsize + t](unitPT, scalePT), digitsToPrint);
			if (needPC)   printValue(out, (*PC)[v](unitPC, scalePC), digitsToPrint);
			if (needE)    printValue(out, (*E)[v*Tsize + t](unitE, scaleE), digitsToPrint);
			if (needET)   printValue(out, (*ET)[v*Tsize + t](unitET, scaleET), digitsToPrint);
			if (needEC)   printValue(out, (*EC)[v](unitEC, scaleEC), digitsToPrint);
			if (needS)    printValue(out, (*S)[v*Tsize + t](unitS, scaleS), digitsToPrint);
			if (needST)   printValue(out, (*ST)[v*Tsize + t](unitST, scaleST), digitsToPrint);
			if (needSC)   printValue(out, (*SC)[v](unitSC, scaleSC), digitsToPrint);
			if (needM)    printValue(out, (*M)[v*Tsize + t](unitM, scaleM), digitsToPrint);
			if (needMT)   printValue(out, (*MT)[v*Tsize + t](unitMT, scaleMT), digitsToPrint);
			if (needMC)   printValue(out, (*MC)[v](unitMC, scaleMC), digitsToPrint);
			if (needQEoS) printValue(out, (*QEoS)[v*Tsize + t](unitQEoS, scaleQEoS), digitsToPrint);
			if (needTEoS) printValue(out, (*TEoS)[v*Tsize + t], digitsToPrint);
			out << std::endl;
	    }
    }
}

void FTTFmodel::calculateAll() {
	Int Vsize = V->size();
	Int Tsize = T->size();
    // generate necessary arrays
	if ( needP) P = new PressVec(Vsize*Tsize);
	if (needPT) {
		if (!needP) P  = new PressVec(Vsize*Tsize); 
		if (!needPC) PC = new PressVec(Vsize);
		PT = new PressVec(Vsize*Tsize);
	}
	if (needPC) PC = new PressVec(Vsize);
	if ( needE) E = new EnergyVec(Vsize*Tsize);
	if (needET) {
		if ( !needE) E  = new EnergyVec(Vsize*Tsize); 
		if (!needEC) EC = new EnergyVec(Vsize);	
		ET = new EnergyVec(Vsize*Tsize);
	}
	if (needEC) EC = new EnergyVec(Vsize);
	if ( needS) S = new EntropyVec(Vsize*Tsize);
	if (needST) {
		if ( !needS) S  = new EntropyVec(Vsize*Tsize); 
		if (!needSC) SC = new EntropyVec(Vsize);	
		ST = new EntropyVec(Vsize*Tsize);
	}
	if (needSC) SC = new EntropyVec(Vsize);
	if ( needM) M = new ChemPotVec(Vsize*Tsize);
	if (needMT) {
		if ( !needM) M  = new ChemPotVec(Vsize*Tsize); 
		if (!needMC) MC = new ChemPotVec(Vsize);	
		MT = new ChemPotVec(Vsize*Tsize);
	}
	if (needMC) MC = new ChemPotVec(Vsize);
	if (needQEoS) {
		if (!needP) P  = new PressVec(Vsize*Tsize); 
		if (!needE) E  = new EnergyVec(Vsize*Tsize); 
		QEoS = new EnergyVec(Vsize*Tsize);
	}
	if (needTEoS) {
		if (!needP) P  = new PressVec(Vsize*Tsize); 
		TEoS = new DoubleVec(Vsize*Tsize);
	}

	energySolver.SetOutput(energyData);
	energySolver.SetTolerance(0.0, eps/10);
	coldEnergySolver.SetOutput(coldEnergyData);
	coldEnergySolver.SetTolerance(0.0, eps/10);
	entropySolver.SetOutput(entropyData);
	entropySolver.SetTolerance(0.0, eps/10);
	coldEntropySolver.SetOutput(coldEnergyData);
	coldEntropySolver.SetTolerance(0.0, eps/10);

	// calculate potential only one time per point
	for (Int v = 0; v < Vsize; ++v) {
		calculatedPC = false;
		calculatedEC = false;
		calculatedSC = false;
		calculatedMC = false;
		coldPhi.setParameters((*V)[v](), coldT(), Z);
		if (needPC) calculatePC(v); 
		if (needEC) calculateEC(v); 
		if (needSC) calculateSC(v); 
		if (needMC) calculateMC(v); 
		for (Int t = 0; t < Tsize; ++t) {
			phi.setParameters((*V)[v](), (*T)[t](), Z);
			calculatedP = false;
			calculatedPT = false;
			calculatedE = false;
			calculatedET = false;
			calculatedS = false;
			calculatedST = false;
			calculatedM = false;
			calculatedMT = false;
			calculatedQEoS = false;
			calculatedTEoS = false;
			if (needP) calculateP(v, t); 
			if (needPT) calculatePT(v, t); 
			if (needE) calculateE(v, t); 
			if (needET) calculateET(v, t); 
			if (needS) calculateS(v, t); 
			if (needST) calculateST(v, t); 
			if (needM) calculateM(v, t); 
			if (needMT) calculateMT(v, t); 
			if (needQEoS) calculateQEoS(v, t); 
			if (needTEoS) calculateTEoS(v, t); 
		}
	}
}

PressVec& FTTFmodel::getP()    { if (!calculatedP)  calculateP();  return *P; }
PressVec& FTTFmodel::getPT()   { if (!calculatedPT) calculatePT(); return *PT; }
PressVec& FTTFmodel::getPC()   { if (!calculatedPC) calculatePC(); return *PC; }

EnergyVec& FTTFmodel::getE()   { if (!calculatedE)  calculateE();  return *E; }
EnergyVec& FTTFmodel::getET()  { if (!calculatedET) calculateET(); return *ET; }
EnergyVec& FTTFmodel::getEC()  { if (!calculatedEC) calculateEC(); return *EC; }

EntropyVec& FTTFmodel::getS()  { if (!calculatedS)  calculateS();  return *S; }
EntropyVec& FTTFmodel::getST() { if (!calculatedST) calculateST(); return *ST; }
EntropyVec& FTTFmodel::getSC() { if (!calculatedSC) calculateSC(); return *SC; }

ChemPotVec& FTTFmodel::getM()  { if (!calculatedM)  calculateM();  return *M; }
ChemPotVec& FTTFmodel::getMT() { if (!calculatedMT) calculateMT(); return *MT; }
ChemPotVec& FTTFmodel::getMC() { if (!calculatedMC) calculateMC(); return *MC; }

void FTTFmodel::calculateP() {
	Int Vsize = V->size();
	Int Tsize = T->size();
	Int Psize = Vsize*Tsize;
	P = new PressVec(Psize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			phi.setParameters((*V)[v], (*T)[t], Z);
			calculateP(v, t);
		}
	}
	calculatedP  = true;
}
void FTTFmodel::calculatePT() {
	if (!calculatedP) calculateP();
	if (!calculatedPC) calculatePC();
	Int PTsize = P->size();
	Int Vsize = V->size();
	Int Tsize = T->size();
	PT = new PressVec(PTsize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			calculatePT(v, t);
		}
	}
	calculatedPT = true;
}
void FTTFmodel::calculatePC() {
	Int Vsize = V->size();
	PC = new PressVec(Vsize);
	for (Int v = 0; v < Vsize; ++v) {
		coldPhi.setParameters((*V)[v], coldT, Z);
		calculatePC(v);
	}
	calculatedPC = true;
}
void FTTFmodel::calculateE() {
	Int Vsize = V->size();
	Int Tsize = T->size();
    E = new EnergyVec(Vsize*Tsize);

	energySolver.SetOutput(energyData);
	energySolver.SetTolerance(0.0, eps/10);

    for (Int v = 0; v < Vsize; ++v) {
    	for (Int t = 0; t < Tsize; ++t) {
    		phi.setParameters((*V)[v], (*T)[t], Z);
    		calculateE(v, t);
    	}
    }
	calculatedE  = true;
}
void FTTFmodel::calculateET() {
	if (!calculatedE) calculateE();
	if (!calculatedEC) calculateEC();
	Int ETsize = E->size();
	Int Vsize = V->size();
	Int Tsize = T->size();
	ET = new EnergyVec(ETsize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			calculateET(v, t);
		}
	}
	calculatedET = true;
}
void FTTFmodel::calculateEC() {
	Int Vsize = V->size();
	EC = new EnergyVec(Vsize);

	coldEnergySolver.SetOutput(coldEnergyData);
	coldEnergySolver.SetTolerance(0.0, eps/10);

	for (Int v = 0; v < Vsize; ++v) {
		coldPhi.setParameters((*V)[v], coldT, Z);
		calculateEC(v);
	}
	calculatedEC = true;
}
void FTTFmodel::calculateS() {
	Int Vsize = V->size();
	Int Tsize = T->size();
    S = new EntropyVec(Vsize*Tsize);

	entropySolver.SetOutput(entropyData);
	entropySolver.SetTolerance(0.0, eps/10);

    for (Int v = 0; v < Vsize; ++v) {
    	for (Int t = 0; t < Tsize; ++t) {
    		phi.setParameters((*V)[v], (*T)[t], Z);
    		calculateS(v, t);
    	}
    }
	calculatedS  = true;
}
void FTTFmodel::calculateST() {
	if (!calculatedS) calculateS();
	if (!calculatedSC) calculateSC();
	Int STsize = S->size();
	Int Vsize = V->size();
	Int Tsize = T->size();
	ST = new EntropyVec(STsize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			calculateST(v, t);
		}
	}
	calculatedST = true;
}
void FTTFmodel::calculateSC() {
	Int Vsize = V->size();
	SC = new EntropyVec(Vsize);

	coldEntropySolver.SetOutput(coldEntropyData);
	coldEntropySolver.SetTolerance(0.0, eps/10);

	for (Int v = 0; v < Vsize; ++v) {
		coldPhi.setParameters((*V)[v], coldT, Z);
		calculateSC(v);
	}
	calculatedSC = true;
}
void FTTFmodel::calculateM() {
	Int Vsize = V->size();
	Int Tsize = T->size();
	Int Msize = Vsize*Tsize;
	M = new ChemPotVec(Msize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			phi.setParameters((*V)[v], (*T)[t], Z);
			calculateM(v, t);
		}
	}
	calculatedM  = true;
}
void FTTFmodel::calculateMT() {
	if (!calculatedM) calculateM();
	if (!calculatedMC) calculateMC();
	Int MTsize = M->size();
	Int Vsize = V->size();
	Int Tsize = T->size();
	MT = new ChemPotVec(MTsize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			calculateMT(v, t);
		}
	}
	calculatedMT = true;
}
void FTTFmodel::calculateMC() {
	Int Vsize = V->size();
	MC = new ChemPotVec(Vsize);
	for (Int v = 0; v < Vsize; ++v) {
		coldPhi.setParameters((*V)[v], coldT, Z);
		calculateMC(v);
	}
	calculatedMC = true;
}
void FTTFmodel::calculateQEoS() {
	if (!calculatedE) calculateE();
	if (!calculatedP) calculateP();
	Int Qsize = E->size();
	Int Vsize = V->size();
	Int Tsize = T->size();
	QEoS = new EnergyVec(Qsize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			calculateQEoS(v, t);
		}
	}
	calculatedQEoS = true;
}
void FTTFmodel::calculateTEoS() {
	if (!calculatedP) calculateP();
	Int TEoSsize = P->size();
	Int Vsize = V->size();
	Int Tsize = T->size();
	TEoS = new DoubleVec(TEoSsize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			calculateTEoS(v, t);
		}
	}
	calculatedTEoS = true;
}
///////////////////////////////////////////////////////////////////////

void FTTFmodel::calculateP(Int v, Int t) {
	Double currentP;
	Int Tsize = T->size();
	currentP = pow(2*(*T)[t](), 5.0/2.0)/6.0/M_PI/M_PI
    *FD3half(phi.valueAt_1());
	(*P)[v*Tsize + t].setValue(currentP);
	calculatedP  = true;
}
void FTTFmodel::calculatePT(Int v, Int t) {
	if (!calculatedP) calculateP(v, t);
	if (!calculatedPC) calculatePC(v);
	Double currentPT;
	Int Tsize = T->size();
	currentPT = (*P)[v*Tsize + t]() - (*PC)[v]();
	(*PT)[v*Tsize + t].setValue(currentPT);
	calculatedPT = true;
}
void FTTFmodel::calculatePC(Int v) {
	Double currentPC;
	currentPC = pow(2*coldT(), 5.0/2.0)/6.0/M_PI/M_PI
    *FD3half(coldPhi.valueAt_1());
	(*PC)[v].setValue(currentPC);
	calculatedPC = true;
}
void FTTFmodel::calculateE(Int v, Int t) {
	Double currentE;
    DoubleVec startE(0.0, rhsFTTFenergy::dim);
    Double a, b;
    Double V1, T1;
    Double xFrom = 1.0;
    Double xTo = 0.0;
    Int Tsize = T->size();

	V1 = (*V)[v]()*Z;
	T1 = (*T)[t]()*pow(Z, -4.0/3.0);

    startE[0] = phi.valueAt_1();
    startE[1] = startE[0];
    startE[2] = 2.0*sqrt(2.0)*V1*pow(T1, 5.0/2.0)/M_PI/M_PI
                * FD3half(phi.valueAt_1()) + 0.76874512422;

    a =   pow(2.0, 7.0/6.0)
        * pow(3.0, 2.0/3.0)
        * pow(M_PI, -5.0/3.0)
        * sqrt(T1)*pow(V1, 2.0/3.0);
    b = 3.0*sqrt(2.0)*V1*pow(T1, 5.0/2.0)/M_PI/M_PI;
	rhsEnergy.updateParameters(a, b);

	energySolver.Integrate(rhsEnergy, startE, xFrom, xTo);
	currentE = startE[2];
	(*E)[v*Tsize + t].setValue(currentE*pow(Z, 7.0/3.0));

	calculatedE  = true;
}
void FTTFmodel::calculateET(Int v, Int t) {
	if (!calculatedE) calculateE(v, t);
	if (!calculatedEC) calculateEC(v);
	Double currentET;
	Int Tsize = T->size();
	currentET = (*E)[v*Tsize + t]() - (*EC)[v]();
	(*ET)[v*Tsize + t].setValue(currentET);
	calculatedET = true;
}
void FTTFmodel::calculateEC(Int v) {
	Double currentEC;
	DoubleVec startE(0.0, rhsFTTFenergy::dim);
	Double a, b;
	Double V1, T1;
	Double xFrom = 1.0;
	Double xTo = 0.0;
	
	T1 = coldT()*pow(Z, -4.0/3.0);
	V1 = (*V)[v]()*Z;
	
	startE[0] = coldPhi.valueAt_1();
	startE[1] = startE[0];
	startE[2] = 2.0*sqrt(2.0)*V1*pow(T1, 5.0/2.0)/M_PI/M_PI
                * FD3half(coldPhi.valueAt_1()) + 0.76874512422;

	a =   pow(2.0, 7.0/6.0)
		* pow(3.0, 2.0/3.0)
		* pow(M_PI, -5.0/3.0)
		* sqrt(T1)*pow(V1, 2.0/3.0);
	b = 3.0*sqrt(2.0)*V1*pow(T1, 5.0/2.0)/M_PI/M_PI;
	rhsColdEnergy.updateParameters(a, b);

	coldEnergySolver.Integrate(rhsColdEnergy, startE, xFrom, xTo);
	currentEC = startE[2];
	(*EC)[v].setValue(currentEC*pow(Z, 7.0/3.0));

	calculatedEC = true;
}
void FTTFmodel::calculateS(Int v, Int t) {
	Double currentS;
    DoubleVec startS(0.0, rhsFTTFentropy::dim);
    Double a, b;
    Double V1, T1;
    Double xFrom = 1.0;
    Double xTo = 0.0;
    Int Tsize = T->size();

	V1 = (*V)[v]()*Z;
	T1 = (*T)[t]()*pow(Z, -4.0/3.0);
	
	startS[0] = phi.valueAt_1();
    startS[1] = startS[0];
	startS[2] = 4.0*sqrt(2.0*T1)*V1*T1/M_PI/M_PI
    		    * FD3half(phi.valueAt_1())
	            - phi.derivativeAt_0();

    a =   pow(2.0, 7.0/6.0)
        * pow(3.0, 2.0/3.0)
        * pow(M_PI, -5.0/3.0)
        * sqrt(T1)*pow(V1, 2.0/3.0);
    b = 7.0*sqrt(2.0*T1)*V1*T1/M_PI/M_PI;
	rhsEntropy.updateParameters(a, b);

	entropySolver.Integrate(rhsEntropy, startS, xFrom, xTo);
	currentS = startS[2];
	(*S)[v*Tsize + t].setValue(currentS*Z);
    	
	calculatedS  = true;
}
void FTTFmodel::calculateST(Int v, Int t) {
	if (!calculatedS) calculateS(v, t);
	if (!calculatedSC) calculateSC(v);
	Int Tsize = T->size();
	Double currentST;
	currentST = (*S)[v*Tsize + t]() - (*SC)[v]();
	(*ST)[v*Tsize + t].setValue(currentST);
	
	calculatedST = true;
}
void FTTFmodel::calculateSC(Int v) {
	Double currentSC;
	DoubleVec startS(0.0, rhsFTTFentropy::dim);
	Double a, b;
	Double V1, T1;
	Double xFrom = 1.0;
	Double xTo = 0.0;
	
	T1 = coldT()*pow(Z, -4.0/3.0);
	V1 = (*V)[v]()*Z;
	
	startS[0] = coldPhi.valueAt_1();
	startS[1] = startS[0];
	startS[2] = 4.0*sqrt(2.0*T1)*V1*T1/M_PI/M_PI
				* FD3half(coldPhi.valueAt_1())
				- coldPhi.derivativeAt_0();

	a =   pow(2.0, 7.0/6.0)
	    * pow(3.0, 2.0/3.0)
	    * pow(M_PI, -5.0/3.0)
	    * sqrt(T1)*pow(V1, 2.0/3.0);
	b = 7.0*sqrt(2.0*T1)*V1*T1/M_PI/M_PI;
	rhsColdEntropy.updateParameters(a, b);

	coldEntropySolver.Integrate(rhsColdEntropy, startS, xFrom, xTo);
	currentSC = startS[2];
	(*SC)[v].setValue(currentSC*Z);

	calculatedSC = true;
}
void FTTFmodel::calculateM(Int v, Int t) {
	Double currentM;
	currentM = (*T)[t]()*phi.valueAt_1();
	Int Tsize = T->size();
	(*M)[v*Tsize + t].setValue(currentM);
	calculatedM  = true;
}
void FTTFmodel::calculateMT(Int v, Int t) {
	if (!calculatedM) calculateM(v, t);
	if (!calculatedMC) calculateMC(v);
	Double currentMT;
	Int Tsize = T->size();
	currentMT = (*M)[v*Tsize + t]() - (*MC)[v]();
	(*MT)[v*Tsize + t].setValue(currentMT);
	calculatedMT = true;
}
void FTTFmodel::calculateMC(Int v) {
	Double currentMC;
	currentMC = coldT()*coldPhi.valueAt_1();
	(*MC)[v].setValue(currentMC);
	calculatedMC = true;
}
void FTTFmodel::calculateQEoS(Int v, Int t) {
	if (!calculatedE) calculateE(v, t);
	if (!calculatedP) calculateP(v, t);
	Double currentQ;
	Int Tsize = T->size();
	currentQ = (*E)[v*Tsize + t]() - (*P)[v*Tsize + t]()*(*V)[v]();
	(*QEoS)[v*Tsize + t].setValue(currentQ);
	calculatedQEoS = true;
}
void FTTFmodel::calculateTEoS(Int v, Int t) {
	if (!calculatedP) calculateP(v, t);
	Int Tsize = T->size();
	Double currentTEoS;
	currentTEoS = (*P)[v*Tsize + t]()*(*V)[v]()/(*T)[t]();
	(*TEoS)[v*Tsize + t] = currentTEoS;
	calculatedTEoS = true;
}