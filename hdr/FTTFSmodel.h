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
#include "../hdr/FTTFSmodel.h"
#include "../hdr/Units.h"
#include "../hdr/stringUtils.h"
#include <fstream>
#include <iostream>

class FTTFSmodel {
public:
	/**
	* @brief A constructor of Thomas-Fermi model with shell corrections.
	*/
    FTTFSmodel(Double _Z = 1.0, Double _Mass = 1.0);
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
	* - M[Units, Scaling] - full FTTF chemical potential, units: eV, Atomic;
	* - MT[Units, Scaling] - thermal FTTF chemcial potential;
	* - MC[Units, Scaling] - cold FTTF chemical potential;
	*/
    void prepareOutput(std::string inputString);
    /**
	* @brief Print output data in file as defined in inputString.
	*/
    void printOutput(const char* filename);

	PressVec&   getLP();
	EnergyVec&  getLE();
	ChemPotVec& getLM();

	PressVec&   getNLP();
	EnergyVec&  getNLE();
	ChemPotVec& getNLM();

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
		rhsEnergyInt(Double _a = 0) : a(_a) {}
		void updateParameter(const Double _a) { a = _a; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
    		dydx[0] = y[1];
		    if (x > 0) {
		        dydx[1] = a*x*FDhalf(y[0] / x);
		        dydx[2] = -x*y[0] * FDmhalf(y[0] / x);
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
		rhsCPInt(Double _a = 0) : a(_a) {}
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
	TempVec* T; Unit unitT; Scaling scaleT; bool needT;
	VolVec*  V; Unit unitV; Scaling scaleV; bool needV;
	DensVec* D; Unit unitD; Scaling scaleD; bool needD;
	ConVec*  C; Unit unitC; Scaling scaleC; bool needC;
	// output data:
	// shell corrections
	PressVec*     LP;  Unit  unitLP; Scaling scaleLP;  bool calculatedLP;  bool needLP;
	EnergyVec*    LE;  Unit  unitLE; Scaling scaleLE;  bool calculatedLE;  bool needLE; 
	ChemPotVec*   LM;  Unit  unitLM; Scaling scaleLM;  bool calculatedLM;  bool needLM; 

	PressVec*     NLP; Unit unitNLP; Scaling scaleNLP; bool calculatedNLP; bool needNLP;
	EnergyVec*    NLE; Unit unitNLE; Scaling scaleNLE; bool calculatedNLE; bool needNLE; 
	ChemPotVec*   NLM; Unit unitNLM; Scaling scaleNLM; bool calculatedNLM; bool needNLM; 

	DoubleVec* CPInt;     bool calculatedCPInt;
	DoubleVec* EnergyInt; bool calculatedEnergyInt;

	// calculate full vector   |  calculate single point
	void         calculateLP(); void calculateLP(Int v, Int t); 
	void         calculateLE(); void calculateLE(Int v, Int t); 
	void         calculateLM(); void calculateLM(Int v, Int t); 

	void        calculateNLP(); void calculateNLP(Int v, Int t); 
	void        calculateNLE(); void calculateNLE(Int v, Int t); 
	void        calculateNLM(); void calculateNLM(Int v, Int t); 

	void      calculateCPInt(); void calculateCPInt(Int v, Int t);
	void  calculateEnergyInt(); void calculateEnergyInt(Int v, Int t);

	void calculateAll(); // uses fast calculation when prepare output
	void printString(std::ofstream &out, std::string data);
	void printValue(std::ofstream &out, Double value, Int digitsToPrint);

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
};

FTTFSmodel::FTTFSmodel(Double _Z, Double _Mass) : 
	Z(_Z), Mass(_Mass), eps(1e-6),
	energyIntData(-1), // save all steps
	energyIntSolver(1e-6, 0.0),
	rhsEnergyInt(0.0),
	CPIntData(-1), // save all steps
	CPIntSolver(1e-6, 0.0),
	rhsCPInt(0.0),
	eStates(7) {
		T = NULL;    unitT = Atomic;    scaleT = lin;    needT = false;
		V = NULL;    unitV = Atomic;    scaleV = lin;    needV = false;
		D = NULL;    unitD = gOverCmc;  scaleD = lin;    needD = false;
		C = NULL;    unitC = Atomic;    scaleC = lin;    needC = false;
		// FTTF
		LP  = NULL;  unitLP = Atomic;  scaleLP = lin;  needLP = false;  calculatedLP = false;
		LE  = NULL;  unitLE = Atomic;  scaleLE = lin;  needLE = false;  calculatedLE = false;
		LM  = NULL;  unitLM = Atomic;  scaleLM = lin;  needLM = false;  calculatedLM = false;

		NLP  = NULL; unitNLP = Atomic; scaleNLP = lin; needNLP = false; calculatedNLP = false;
		NLE  = NULL; unitNLE = Atomic; scaleNLE = lin; needNLE = false; calculatedNLE = false;
		NLM  = NULL; unitNLM = Atomic; scaleNLM = lin; needNLM = false; calculatedNLM = false;

		CPInt = NULL;     calculatedCPInt = false;
		EnergyInt = NULL; calculatedEnergyInt = false;
}

void FTTFSmodel::setParameters(const VolVec &_V, const TempVec &_T) {
	if (V != NULL) delete V;
	if (T != NULL) delete T;
	V = new VolVec(_V);
	T = new TempVec(_T);
	transformVtoC();
	transformVtoD();
	needV = true;
	needT = true;
	calculatedLP  = false; if (LP != NULL) delete LP;
	calculatedLE  = false; if (LE != NULL) delete LE;
	calculatedLM  = false; if (LM != NULL) delete LM;
	calculatedNLP  = false; if (NLP != NULL) delete NLP;
	calculatedNLE  = false; if (NLE != NULL) delete NLE;
	calculatedNLM  = false; if (NLM != NULL) delete NLM;
	calculatedCPInt = false; if (CPInt != NULL) delete CPInt;
	calculatedEnergyInt = false; if (EnergyInt != NULL) delete EnergyInt;
}

void FTTFSmodel::setParameters(const DensVec &_D, const TempVec &_T) {
	if (D != NULL) delete D;
	if (T != NULL) delete T;
	D = new DensVec(_D);
	T = new TempVec(_T);
	needD = true;
	needT = true;
	transformDtoV();
	transformVtoC();
	calculatedLP  = false; if (LP != NULL) delete LP;
	calculatedLE  = false; if (LE != NULL) delete LE;
	calculatedLM  = false; if (LM != NULL) delete LM;
	calculatedNLP  = false; if (NLP != NULL) delete NLP;
	calculatedNLE  = false; if (NLE != NULL) delete NLE;
	calculatedNLM  = false; if (NLM != NULL) delete NLM;
	calculatedCPInt = false; if (CPInt != NULL) delete CPInt;
	calculatedEnergyInt = false; if (EnergyInt != NULL) delete EnergyInt;
}

void FTTFSmodel::setParameters(const ConVec &_C, const TempVec &_T) {
	if (C != NULL) delete C;
	if (T != NULL) delete T;
	C = new ConVec(_C);
	T = new TempVec(_T);
	needC = true;
	needT = true;
	transformCtoV();
	transformVtoD();
	calculatedLP  = false; if (LP != NULL) delete LP;
	calculatedLE  = false; if (LE != NULL) delete LE;
	calculatedLM  = false; if (LM != NULL) delete LM;
	calculatedNLP  = false; if (NLP != NULL) delete NLP;
	calculatedNLE  = false; if (NLE != NULL) delete NLE;
	calculatedNLM  = false; if (NLM != NULL) delete NLM;
	calculatedCPInt = false; if (CPInt != NULL) delete CPInt;
	calculatedEnergyInt = false; if (EnergyInt != NULL) delete EnergyInt;
}

void FTTFSmodel::transformDtoV() {
	Int n = D->size();
	V = new VolVec(n);
	Double currentV;
	for (Int i = 0; i < n; ++i) {
		currentV = Mass/(Avogadro*1.4818e-25)/(*D)[i](gOverCmc);
		(*V)[i].setValue(currentV);
	}
}

void FTTFSmodel::transformCtoV() {
	Int n = C->size();
	V = new VolVec(n);
	Double currentV;
	for (Int i = 0; i < n; ++i) {
		currentV = 1.0/(*C)[i]();
		(*V)[i].setValue(currentV);
	}
}

void FTTFSmodel::transformVtoD() {
	Int n = V->size();
	D = new DensVec(n);
	Double currentD;
	for (Int i = 0; i < n; ++i) {
		currentD = Mass/(Avogadro*1.4818e-25)/(*V)[i]();
		(*D)[i].setValue(currentD, gOverCmc);
	}
}

void FTTFSmodel::transformVtoC() {
	Int n = V->size();
	C = new ConVec(n);
	Double currentC;
	for (Int i = 0; i < n; ++i) {
		currentC = 1.0/(*V)[i]();
		(*C)[i].setValue(currentC);
	}
}

void FTTFSmodel::prepareOutput(std::string inputString) {
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
    	if (!quantity.compare("V"))   { needV   = true; if (!unit.empty()) unitV   = stringToUnit(unit); if (!scale.empty()) scaleV   = stringToScale(scale); }
    	if (!quantity.compare("D"))   { needD   = true; if (!unit.empty()) unitD   = stringToUnit(unit); if (!scale.empty()) scaleD   = stringToScale(scale); }
 		if (!quantity.compare("C"))   { needC   = true; if (!unit.empty()) unitC   = stringToUnit(unit); if (!scale.empty()) scaleC   = stringToScale(scale); }
    	if (!quantity.compare("T"))   { needT   = true; if (!unit.empty()) unitT   = stringToUnit(unit); if (!scale.empty()) scaleT   = stringToScale(scale); }
    	if (!quantity.compare("NP"))  { needLP  = true; if (!unit.empty()) unitLP  = stringToUnit(unit); if (!scale.empty()) scaleLP  = stringToScale(scale); }
		if (!quantity.compare("NE"))  { needLE  = true; if (!unit.empty()) unitLE  = stringToUnit(unit); if (!scale.empty()) scaleLE  = stringToScale(scale); }
		if (!quantity.compare("NM"))  { needLM  = true; if (!unit.empty()) unitLM  = stringToUnit(unit); if (!scale.empty()) scaleLM  = stringToScale(scale); }
		if (!quantity.compare("NLP")) { needNLP = true; if (!unit.empty()) unitNLP = stringToUnit(unit); if (!scale.empty()) scaleNLP = stringToScale(scale); }
		if (!quantity.compare("NLE")) { needNLE = true; if (!unit.empty()) unitNLE = stringToUnit(unit); if (!scale.empty()) scaleNLE = stringToScale(scale); }
		if (!quantity.compare("NLM")) { needNLM = true; if (!unit.empty()) unitNLM = stringToUnit(unit); if (!scale.empty()) scaleNLM = stringToScale(scale); }
    }
    calculateAll();
}

void FTTFSmodel::printString(std::ofstream &out, std::string data) {
	out.setf(std::ios::left);
    out.width(20);
    out.fill(' ');
    out << data;
}

void FTTFSmodel::printValue(std::ofstream &out, Double value, Int digitsToPrint) {
	out.setf(std::ios::scientific | std::ios::left | std::ios::showpos);
    out.width(20);
    out.precision(digitsToPrint);
    out.fill(' ');
    out << value;
}

void FTTFSmodel::printOutput(const char* filename) {
	std::ofstream out;
    Int digitsToPrint = static_cast<int>(-log10(eps));
    Int Vsize = V->size();
    Int Tsize = T->size();
    out.open(filename, std::ios::out);
    // print first line:
    std::string s;
    if (needV)   { s = "V[";   s += unitToString(unitV);   s += ","; s += scaleToString(scaleV);   s += "]"; printString(out, s); }
	if (needD)   { s = "D[";   s += unitToString(unitD);   s += ","; s += scaleToString(scaleD);   s += "]"; printString(out, s); }
	if (needC)   { s = "C[";   s += unitToString(unitC);   s += ","; s += scaleToString(scaleC);   s += "]"; printString(out, s); }
	if (needT)   { s = "T[";   s += unitToString(unitT);   s += ","; s += scaleToString(scaleT);   s += "]"; printString(out, s); }
	if (needLP)  { s = "LP[";  s += unitToString(unitLP);  s += ","; s += scaleToString(scaleLP);  s += "]"; printString(out, s); }
	if (needLE)  { s = "LE[";  s += unitToString(unitLE);  s += ","; s += scaleToString(scaleLE);  s += "]"; printString(out, s); }
	if (needLM)  { s = "LM[";  s += unitToString(unitLM);  s += ","; s += scaleToString(scaleLM);  s += "]"; printString(out, s); }
	if (needNLP) { s = "NLP["; s += unitToString(unitNLP); s += ","; s += scaleToString(scaleNLP); s += "]"; printString(out, s); }
	if (needNLE) { s = "NLE["; s += unitToString(unitNLE); s += ","; s += scaleToString(scaleNLE); s += "]"; printString(out, s); }
	if (needNLM) { s = "NLM["; s += unitToString(unitNLM); s += ","; s += scaleToString(scaleNLM); s += "]"; printString(out, s); }
	out << std::endl;
	// print data
    for (Int v = 0; v < Vsize; ++v) {
	    for (Int t = 0; t < Tsize; ++t) {
	    	if (needV)   printValue(out, (*V  )[v]          (unitV,   scaleV),   digitsToPrint);
	    	if (needD)   printValue(out, (*D  )[v]          (unitD,   scaleD),   digitsToPrint);
	    	if (needC)   printValue(out, (*C  )[v]          (unitC,   scaleC),   digitsToPrint);
	    	if (needT)   printValue(out, (*T  )[t]          (unitT,   scaleT),   digitsToPrint);
			if (needLP)  printValue(out, (*LP )[v*Tsize + t](unitLP,  scaleLP),  digitsToPrint);
			if (needLE)  printValue(out, (*LE )[v*Tsize + t](unitLE,  scaleLE),  digitsToPrint);
			if (needLM)  printValue(out, (*LM )[v*Tsize + t](unitLM,  scaleLM),  digitsToPrint);
			if (needNLP) printValue(out, (*NLP)[v*Tsize + t](unitNLP, scaleNLP), digitsToPrint);
			if (needNLE) printValue(out, (*NLE)[v*Tsize + t](unitNLE, scaleNLE), digitsToPrint);
			if (needNLM) printValue(out, (*NLM)[v*Tsize + t](unitNLM, scaleNLM), digitsToPrint);
			out << std::endl;
	    }
    }
}

void FTTFSmodel::calculateAll() {
	Int Vsize = V->size();
	Int Tsize = T->size();
    // generate necessary arrays
	if ( needLP) LP  = new PressVec(Vsize*Tsize);
	if ( needLE) LE  = new EnergyVec(Vsize*Tsize);
	LM  = new ChemPotVec(Vsize*Tsize);

	if (needNLP) NLP = new PressVec(Vsize*Tsize);
	if (needNLE) NLE = new EnergyVec(Vsize*Tsize);
	NLM = new ChemPotVec(Vsize*Tsize);

	CPInt = new DoubleVec(Vsize*Tsize);
	EnergyInt = new DoubleVec(Vsize*Tsize);

	energyIntSolver.SetOutput(energyData);
	energyIntSolver.SetTolerance(0.0, eps/10);
	CPIntSolver.SetOutput(energyData);
	CPIntSolver.SetTolerance(0.0, eps/10);

	// calculate potential only one time per point
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			phi.setParameters((*V)[v](), (*T)[t](), Z);
			eStates.setParameters((*V)[v](), (*T)[t](), Z);
			calculatedCPInt = false;
			calculatedEnergyInt = false;
			calculatedLP  = false;
			calculatedLE  = false;
			calculatedLM  = false;
			calculatedNLP = false;
			calculatedNLE = false;
			calculatedNLM = false;
			calculateCPInt(v, t);
			calculateEnergyInt(v, t);
			calculateLM(v, t); 
			calculateNLM(v, t); 
			if (needLP)  calculateLP(v, t); 
			if (needLE)  calculateLE(v, t); 
			if (needNLP) calculateNLP(v, t); 
			if (needNLE) calculateNLE(v, t); 
		}
	}
}

PressVec&   FTTFSmodel::getLP()  { if (!calculatedLP)  calculateLP();  return *LP; }
EnergyVec&  FTTFSmodel::getLE()  { if (!calculatedLE)  calculateLE();  return *LE; }
ChemPotVec& FTTFSmodel::getLM()  { if (!calculatedLM)  calculateLM();  return *LM; }

PressVec&   FTTFSmodel::getNLP() { if (!calculatedNLP) calculateNLP(); return *NLP; }
EnergyVec&  FTTFSmodel::getNLE() { if (!calculatedNLE) calculateNLE(); return *NLE; }
ChemPotVec& FTTFSmodel::getNLM() { if (!calculatedNLM) calculateNLM(); return *NLM; }

void FTTFSmodel::calculateLP() {
	if (!calculatedLM) calculateLM();
	Int Vsize = V->size();
	Int Tsize = T->size();
	Int Psize = Vsize*Tsize;
	LP = new PressVec(Psize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			phi.setParameters((*V)[v], (*T)[t], Z);
			eStates.setParameters((*V)[v](), (*T)[t](), Z);
			calculateLP(v, t);
		}
	}
	calculatedLP  = true;
}

void FTTFSmodel::calculateLE() {
	if (!calculatedLM) calculateLM();
	if (!calculatedEnergyInt) calculateEnergyInt();
	Int Vsize = V->size();
	Int Tsize = T->size();
    LE = new EnergyVec(Vsize*Tsize);

	energySolver.SetOutput(energyData);
	energySolver.SetTolerance(0.0, eps/10);

    for (Int v = 0; v < Vsize; ++v) {
    	for (Int t = 0; t < Tsize; ++t) {
    		phi.setParameters((*V)[v], (*T)[t], Z);
    		eStates.setParameters((*V)[v](), (*T)[t](), Z);
    		calculateLE(v, t);
    	}
    }
	calculatedLE  = true;
}

void FTTFSmodel::calculateLM() {
	if (!calculatedCPInt) calculateCPInt();
	Int Vsize = V->size();
	Int Tsize = T->size();
	Int Msize = Vsize*Tsize;
	LM = new ChemPotVec(Msize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			phi.setParameters((*V)[v], (*T)[t], Z);
			eStates.setParameters((*V)[v](), (*T)[t](), Z);
			calculateLM(v, t);
		}
	}
	calculatedLM  = true;
}

void FTTFSmodel::calculateCPInt() {
	Int Vsize = V->size();
	Int Tsize = T->size();
	CPInt = new DoubleVec(Vsize*Tsize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			phi.setParameters((*V)[v], (*T)[t], Z);
			eStates.setParameters((*V)[v](), (*T)[t](), Z);
			calculateCPInt(v, t);
		}
	}
	calculatedCPInt  = true;	
}

void FTTFSmodel::calculateEnergyInt() {
	Int Vsize = V->size();
	Int Tsize = T->size();
	EnergyInt = new DoubleVec(Vsize*Tsize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			phi.setParameters((*V)[v], (*T)[t], Z);
			eStates.setParameters((*V)[v](), (*T)[t](), Z);
			calculateEnergyInt(v, t);
		}
	}
	calculatedEnergyInt  = true;	
}

void FTTFSmodel::calculateNLP() {
	Int Vsize = V->size();
	Int Tsize = T->size();
	Int Psize = Vsize*Tsize;
	NLP = new PressVec(Psize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			phi.setParameters((*V)[v], (*T)[t], Z);
			eStates.setParameters((*V)[v](), (*T)[t](), Z);
			calculateNLP(v, t);
		}
	}
	calculatedNLP  = true;
}

void FTTFSmodel::calculateNLE() {
	Int Vsize = V->size();
	Int Tsize = T->size();
    NLE = new EnergyVec(Vsize*Tsize);

	energySolver.SetOutput(energyData);
	energySolver.SetTolerance(0.0, eps/10);

    for (Int v = 0; v < Vsize; ++v) {
    	for (Int t = 0; t < Tsize; ++t) {
    		phi.setParameters((*V)[v], (*T)[t], Z);
    		eStates.setParameters((*V)[v](), (*T)[t](), Z);
    		calculateNLE(v, t);
    	}
    }
	calculatedNLE  = true;
}

void FTTFSmodel::calculateNLM() {
	Int Vsize = V->size();
	Int Tsize = T->size();
	Int Msize = Vsize*Tsize;
	NLM = new ChemPotVec(Msize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			phi.setParameters((*V)[v], (*T)[t], Z);
			eStates.setParameters((*V)[v](), (*T)[t](), Z);
			calculateNLM(v, t);
		}
	}
	calculatedNLM  = true;
}

///////////////////////////////////////////////////////////////////////

void FTTFSmodel::calculateLP(Int v, Int t) {
    Double density;
	Int Tsize = T->size();
    density = pow((*T)[t](), 1.5)*sqrt(2.0) / M_PI / M_PI
            * FDhalf(phi.valueAt_1());
    (*LP)[v*Tsize + t].setValue(density*(*LM)[v*Tsize + t]());
	calculatedLP = true;
}

void FTTFSmodel::calculateLE(Int v, Int t) {
	Int Tsize = T->size();
	Double value;
    value = (*EnergyInt)[v*Tsize + t];
    value = (1.5*Z - value)*(*LM)[v*Tsize + t]();
	(*LE)[v*Tsize + t].setValue(value);
	calculatedLE  = true;
}

void FTTFSmodel::calculateLM(Int v, Int t) {
    Double value;
    Int Tsize = T->size();
    value = eStates.DNlow() / (*CPint)[v*Tsize + t] * (*T)[t]();
    (*LM)[v*Tsize + t].setValue(value);
	calculatedLM  = true;
}

void FTTFSmodel::calculateCPInt(Int v, Int t) {
	Double currentCPInt;
    DoubleVec startCPInt(0.0, rhsCPInt::dim);
    Double xFrom = 1.0;
    Double xTo = 0.0;
    Int Tsize = T->size();

    startCPInt[0] = phi.valueAt_1();
    startCPInt[1] = startCPInt[0];
    startCPInt[2] = 0.0;

	CPIntSolver.Integrate(rhsCPInt, startCPInt, xFrom, xTo);
	currentCPInt = 3*(*V)[v]() / 4.0 / M_PI * 2 / M_PI*sqrt(2 * (*T)[t]())*startCPInt[2];
	(*CPInt)[v*Tsize + t].setValue(currentCPInt);

	calculatedCPInt  = true;
}

void FTTFSmodel::calculateEnergyInt(Int v, Int t) {
	Double currentEnergyInt;
	DoubleVec startEnergyInt(0.0, rhsEnergyInt::dim);
	Double xFrom = 1.0;
	Double xTo = 0.0;
	Int Tsize = T->size();
	
	startEnergyInt[0] = phi.valueAt_1();
	startEnergyInt[1] = startEnergyInt[0];
	startEnergyInt[2] = 0.0;
	
	EnergyIntSolver.Integrate(rhsEnergyInt, startEnergyInt, xFrom, xTo);
	currentEnergyInt = 3*(*V)[v]() / 4.0 / M_PI *pow(2.0 * (*T)[t](), 1.5) / M_PI*startEnergyInt[2];
	(*EnergyInt)[v*Tsize + t].setValue(currentEnergyInt);
	
	calculatedEnergyInt  = true;
}

void FTTFSmodel::calculateNLP(Int v, Int t) {
    Double density;
	Int Tsize = T->size();
    density = pow((*T)[t](), 1.5)*sqrt(2.0) / M_PI / M_PI
            * FDhalf(phi.valueAt_1());
    (*LP)[v*Tsize + t].setValue(density*(*NLM)[v*Tsize + t]());
	calculatedNLP = true;
}

void FTTFSmodel::calculateNLE(Int v, Int t) {
	Int Tsize = T->size();
	Double value;
    value = (*EnergyInt)[v*Tsize + t];
    value = (1.5*Z - value)*(*NLM)[v*Tsize + t]();
	(*NLE)[v*Tsize + t].setValue(value);
	calculatedNLE  = true;
}

void FTTFSmodel::calculateNLM(Int v, Int t) {
	Double currentN;
    Double exactN = Z;
    Int Tsize = T->size();
    Double MTF = phi.valueAt_1();
    Double leftM = -fabs(MTF * 2);
    Double rightM = fabs(MTF * 3);

    std::cout << "BE = " << boundEnergy << std::endl;
    std::cout << "MTF = " << MTF << std::endl;

    while (abs(leftM - rightM) / abs(leftM + rightM) > eps*10) {
    	eStates.setPhiShift(0.5*(rightM + leftM));
        currentN = eStates.N();
        std::cout << "DM = " << 0.5*(leftM + rightM) << std::endl;
        std::cout << "currentN = " << currentN << std::endl;
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
    }

    (*NLM)[v*Tsize + t].setValue((rightM + leftM)*0.5*(*T)[v*Tsize + t]());
	calculatedNLM = true;
}