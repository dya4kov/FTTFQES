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
#include "../hdr/FTTFmodel.h"
#include "../hdr/Yfunction.h"
#include "../hdr/Units.h"
#include "../hdr/stringUtils.h"
#include <fstream>
#include <iostream>
/**
* @brief This class implements interface for FTTFQE model.
* @details The main formulas for calculation of corrections to thermodynamic quantities
* are presented below and expressed throgh the FTTF potential and its correction.
* - the formula for the quantum and exchange correction to pressure: 
* @f[
*	\Delta P = \frac{3T^2}{\pi^3}I_{1/2}(\phi(1))\psi(1) + Y(\phi(1)).
* @f]
* - the formula for the quantum and exchange correction to chemical potential:
* @f[
*	\Delta \mu = \frac{1}{3\pi}\sqrt{\frac{T}{2}}
*    \left[\frac12 I_{-1/2}(\phi(1)) + \psi(1)\right].
* @f]
* - The methods for calculating energy and entropy corrections are presented inside 
    model class
* - Thermal parts are caculated in the following way:
* @f[
*	\Delta P_T = \Delta P - \left. \Delta P \right|_{T = 0}.
* @f]
* @f[
* \Delta E_T = \Delta E - \left.\Delta E\right|_{T = 0}.
* @f]
* @f[
*	\Delta S_T = \Delta S - \left.\Delta S\right|_{T = 0}.
* @f]
* @f[
*   \Delta \mu_T = \Delta \mu - \left.\Delta \mu\right|_{T = 0}.
* @f]
*/
class FTTFQEmodel {
public:
	/**
	* @brief A constructor of Thomas-Fermi model.
	*/
    FTTFQEmodel(Double _Z = 1.0, Double _Mass = 1.0);
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

    PressVec&   getP();  PressVec&   getDP();  PressVec&   getP0(); 
	PressVec&   getPT(); PressVec&   getDPT(); PressVec&   getPT0();
	PressVec&   getPC(); PressVec&   getDPC(); PressVec&   getPC0();
	EnergyVec&  getE();  EnergyVec&  getDE();  EnergyVec&  getE0();
	EnergyVec&  getET(); EnergyVec&  getDET(); EnergyVec&  getET0();
	EnergyVec&  getEC(); EnergyVec&  getDEC(); EnergyVec&  getEC0();
	EntropyVec& getS();  EntropyVec& getDS();  EntropyVec& getS0();
	EntropyVec& getST(); EntropyVec& getDST(); EntropyVec& getST0();
	EntropyVec& getSC(); EntropyVec& getDSC(); EntropyVec& getSC0();
	ChemPotVec& getM();  ChemPotVec& getDM();  ChemPotVec& getM0();
	ChemPotVec& getMT(); ChemPotVec& getDMT(); ChemPotVec& getMT0();
	ChemPotVec& getMC(); ChemPotVec& getDMC(); ChemPotVec& getMC0();

	EnergyVec&  getQEoS();
	DoubleVec&  getTEoS();
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
    struct rhsFTTFQEenergy {
		Double a;
		Double b;
		static const Int dim = 4;
		FermiDirac<Mhalf> FDmhalf;
		FermiDirac<Half> FDhalf;
		Yfunction Y;
		rhsFTTFQEenergy(Double _a = 0, Double _b = 0) : a(_a), b(_b) {}
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
	struct rhsFTTFQEentropy {
		Double a;
		Double b;
		static const Int dim = 4;
		FermiDirac<Half> FDhalf;
		FermiDirac<Mhalf> FDmhalf;
		Yfunction Y;
		rhsFTTFQEentropy(Double _a = 0, Double _b = 0) : a(_a), b(_b) {}
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
	//input data
	TempVec* T; Unit unitT; Scaling scaleT; bool needT;
	VolVec*  V; Unit unitV; Scaling scaleV; bool needV;
	DensVec* D; Unit unitD; Scaling scaleD; bool needD;
	ConVec*  C; Unit unitC; Scaling scaleC; bool needC;
	// output data:
	//FTTF model
	PressVec*     P0; Unit  unitP0; Scaling  scaleP0; bool  calculatedP0; bool  needP0;
	PressVec*    PT0; Unit unitPT0; Scaling scalePT0; bool calculatedPT0; bool needPT0; 
	PressVec*    PC0; Unit unitPC0; Scaling scalePC0; bool calculatedPC0; bool needPC0; 
	EnergyVec*    E0; Unit  unitE0; Scaling  scaleE0; bool  calculatedE0; bool  needE0; 
	EnergyVec*   ET0; Unit unitET0; Scaling scaleET0; bool calculatedET0; bool needET0; 
	EnergyVec*   EC0; Unit unitEC0; Scaling scaleEC0; bool calculatedEC0; bool needEC0; 
	EntropyVec*   S0; Unit  unitS0; Scaling  scaleS0; bool  calculatedS0; bool  needS0; 
	EntropyVec*  ST0; Unit unitST0; Scaling scaleST0; bool calculatedST0; bool needST0; 
	EntropyVec*  SC0; Unit unitSC0; Scaling scaleSC0; bool calculatedSC0; bool needSC0; 
	ChemPotVec*   M0; Unit  unitM0; Scaling  scaleM0; bool  calculatedM0; bool  needM0; 
	ChemPotVec*  MT0; Unit unitMT0; Scaling scaleMT0; bool calculatedMT0; bool needMT0; 
	ChemPotVec*  MC0; Unit unitMC0; Scaling scaleMC0; bool calculatedMC0; bool needMC0;
	// QE corrections only
	PressVec*     DP; Unit  unitDP; Scaling  scaleDP; bool  calculatedDP; bool  needDP;
	PressVec*    DPT; Unit unitDPT; Scaling scaleDPT; bool calculatedDPT; bool needDPT; 
	PressVec*    DPC; Unit unitDPC; Scaling scaleDPC; bool calculatedDPC; bool needDPC; 
	EnergyVec*    DE; Unit  unitDE; Scaling  scaleDE; bool  calculatedDE; bool  needDE; 
	EnergyVec*   DET; Unit unitDET; Scaling scaleDET; bool calculatedDET; bool needDET; 
	EnergyVec*   DEC; Unit unitDEC; Scaling scaleDEC; bool calculatedDEC; bool needDEC; 
	EntropyVec*   DS; Unit  unitDS; Scaling  scaleDS; bool  calculatedDS; bool  needDS; 
	EntropyVec*  DST; Unit unitDST; Scaling scaleDST; bool calculatedDST; bool needDST; 
	EntropyVec*  DSC; Unit unitDSC; Scaling scaleDSC; bool calculatedDSC; bool needDSC; 
	ChemPotVec*   DM; Unit  unitDM; Scaling  scaleDM; bool  calculatedDM; bool  needDM; 
	ChemPotVec*  DMT; Unit unitDMT; Scaling scaleDMT; bool calculatedDMT; bool needDMT; 
	ChemPotVec*  DMC; Unit unitDMC; Scaling scaleDMC; bool calculatedDMC; bool needDMC;
	// FTTF + QE
	PressVec*     P;  Unit  unitP;  Scaling   scaleP; bool   calculatedP; bool   needP;
	PressVec*    PT;  Unit unitPT;  Scaling  scalePT; bool  calculatedPT; bool  needPT; 
	PressVec*    PC;  Unit unitPC;  Scaling  scalePC; bool  calculatedPC; bool  needPC; 
	EnergyVec*    E;  Unit  unitE;  Scaling   scaleE; bool   calculatedE; bool   needE; 
	EnergyVec*   ET;  Unit unitET;  Scaling  scaleET; bool  calculatedET; bool  needET; 
	EnergyVec*   EC;  Unit unitEC;  Scaling  scaleEC; bool  calculatedEC; bool  needEC; 
	EntropyVec*   S;  Unit  unitS;  Scaling   scaleS; bool   calculatedS; bool   needS; 
	EntropyVec*  ST;  Unit unitST;  Scaling  scaleST; bool  calculatedST; bool  needST; 
	EntropyVec*  SC;  Unit unitSC;  Scaling  scaleSC; bool  calculatedSC; bool  needSC; 
	ChemPotVec*   M;  Unit  unitM;  Scaling   scaleM; bool   calculatedM; bool   needM; 
	ChemPotVec*  MT;  Unit unitMT;  Scaling  scaleMT; bool  calculatedMT; bool  needMT; 
	ChemPotVec*  MC;  Unit unitMC;  Scaling  scaleMC; bool  calculatedMC; bool  needMC;
	
	EnergyVec* QEoS; Unit unitQEoS; Scaling scaleQEoS; bool calculatedQEoS; bool needQEoS; 
	DoubleVec* TEoS;                                   bool calculatedTEoS; bool needTEoS;

	// calculate full vector |      calculate single point
	void        calculateDP();  void    calculateDP(Int v, Int t); 
	void       calculateDPT();  void   calculateDPT(Int v, Int t); 
	void       calculateDPC();  void   calculateDPC(Int v); 
	void        calculateDE();  void    calculateDE(Int v, Int t); 
	void       calculateDET();  void   calculateDET(Int v, Int t); 
	void       calculateDEC();  void   calculateDEC(Int v); 
	void        calculateDS();  void    calculateDS(Int v, Int t); 
	void       calculateDST();  void   calculateDST(Int v, Int t); 
	void       calculateDSC();  void   calculateDSC(Int v); 
	void        calculateDM();  void    calculateDM(Int v, Int t); 
	void       calculateDMT();  void   calculateDMT(Int v, Int t); 
	void       calculateDMC();  void   calculateDMC(Int v); 
	
	void         calculateP();  void    calculateP(Int v, Int t); 
	void        calculatePT();  void   calculatePT(Int v, Int t); 
	void        calculatePC();  void   calculatePC(Int v); 
	void         calculateE();  void    calculateE(Int v, Int t); 
	void        calculateET();  void   calculateET(Int v, Int t); 
	void        calculateEC();  void   calculateEC(Int v); 
	void         calculateS();  void    calculateS(Int v, Int t); 
	void        calculateST();  void   calculateST(Int v, Int t); 
	void        calculateSC();  void   calculateSC(Int v); 
	void         calculateM();  void    calculateM(Int v, Int t); 
	void        calculateMT();  void   calculateMT(Int v, Int t); 
	void        calculateMC();  void   calculateMC(Int v); 

	void      calculateQEoS();  void calculateQEoS(Int v, Int t); 
	void      calculateTEoS();  void calculateTEoS(Int v, Int t); 

	void calculateAll(); // uses fast calculation when prepare output
	void printString(std::ofstream &out, std::string data);
	void printValue(std::ofstream &out, Double value, Int digitsToPrint);

    // ode solvers
    rhsFTTFQEenergy rhsEnergy;
    ODEsolver<ODEstepperPD853<rhsFTTFQEenergy> > energySolver;
    ODEdata energyData;

    rhsFTTFQEentropy rhsEntropy;
    ODEsolver<ODEstepperPD853<rhsFTTFQEentropy> > entropySolver;
    ODEdata entropyData;

    rhsFTTFQEenergy rhsColdEnergy;
    ODEsolver<ODEstepperPD853<rhsFTTFQEenergy> > coldEnergySolver;
    ODEdata coldEnergyData;

    rhsFTTFQEentropy rhsColdEntropy;
    ODEsolver<ODEstepperPD853<rhsFTTFQEentropy> > coldEntropySolver;
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

    FTTFpotential phi;
    FTTFQEpotential psi;

    FTTFpotential coldPhi;
    FTTFQEpotential coldPsi;
    Temperature coldT;

    FTTFmodel fttf;
};

FTTFQEmodel::FTTFQEmodel(Double _Z, Double _Mass) : 
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
	coldT(1e-6),
	fttf(_Z, _Mass) {
		T = NULL;    unitT = Atomic;    scaleT = lin;    needT = false;
		V = NULL;    unitV = Atomic;    scaleV = lin;    needV = false;
		D = NULL;    unitD = gOverCmc;  scaleD = lin;    needD = false;
		C = NULL;    unitC = Atomic;    scaleC = lin;    needC = false;
		// FTTF
		P0  = NULL;  unitP0 = Atomic;  scaleP0 = lin;  needP0 = false;  calculatedP0 = false;
		PT0 = NULL; unitPT0 = Atomic; scalePT0 = lin; needPT0 = false; calculatedPT0 = false;
		PC0 = NULL; unitPC0 = Atomic; scalePC0 = lin; needPC0 = false; calculatedPC0 = false;
		E0  = NULL;  unitE0 = Atomic;  scaleE0 = lin;  needE0 = false;  calculatedE0 = false;
		ET0 = NULL; unitET0 = Atomic; scaleET0 = lin; needET0 = false; calculatedET0 = false;
		EC0 = NULL; unitEC0 = Atomic; scaleEC0 = lin; needEC0 = false; calculatedEC0 = false;
		S0  = NULL;  unitS0 = Atomic;  scaleS0 = lin;  needS0 = false;  calculatedS0 = false;
		ST0 = NULL; unitST0 = Atomic; scaleST0 = lin; needST0 = false; calculatedST0 = false;
		SC0 = NULL; unitSC0 = Atomic; scaleSC0 = lin; needSC0 = false; calculatedSC0 = false;
		M0  = NULL;  unitM0 = Atomic;  scaleM0 = lin;  needM0 = false;  calculatedM0 = false;
		MT0 = NULL; unitMT0 = Atomic; scaleMT0 = lin; needMT0 = false; calculatedMT0 = false;
		MC0 = NULL; unitMC0 = Atomic; scaleMC0 = lin; needMC0 = false; calculatedMC0 = false;
		// QE corrections
		DP  = NULL;  unitDP = Atomic;  scaleDP = lin;  needDP = false;  calculatedDP = false;    
		DPT = NULL; unitDPT = Atomic; scaleDPT = lin; needDPT = false; calculatedDPT = false;   
		DPC = NULL; unitDPC = Atomic; scaleDPC = lin; needDPC = false; calculatedDPC = false;   
		DE  = NULL;  unitDE = Atomic;  scaleDE = lin;  needDE = false;  calculatedDE = false;    
		DET = NULL; unitDET = Atomic; scaleDET = lin; needDET = false; calculatedDET = false;   
		DEC = NULL; unitDEC = Atomic; scaleDEC = lin; needDEC = false; calculatedDEC = false;   
		DS  = NULL;  unitDS = Atomic;  scaleDS = lin;  needDS = false;  calculatedDS = false;    
		DST = NULL; unitDST = Atomic; scaleDST = lin; needDST = false; calculatedDST = false;   
		DSC = NULL; unitDSC = Atomic; scaleDSC = lin; needDSC = false; calculatedDSC = false;   
		DM  = NULL;  unitDM = Atomic;  scaleDM = lin;  needDM = false;  calculatedDM = false;    
		DMT = NULL; unitDMT = Atomic; scaleDMT = lin; needDMT = false; calculatedDMT = false;   
		DMC = NULL; unitDMC = Atomic; scaleDMC = lin; needDMC = false; calculatedDMC = false;   
		// FTTF + QE
		P  = NULL;    unitP = Atomic;   scaleP = lin;   needP = false;   calculatedP = false;
		PT = NULL;   unitPT = Atomic;  scalePT = lin;  needPT = false;  calculatedPT = false;
		PC = NULL;   unitPC = Atomic;  scalePC = lin;  needPC = false;  calculatedPC = false;
		E  = NULL;    unitE = Atomic;   scaleE = lin;   needE = false;   calculatedE = false;
		ET = NULL;   unitET = Atomic;  scaleET = lin;  needET = false;  calculatedET = false;
		EC = NULL;   unitEC = Atomic;  scaleEC = lin;  needEC = false;  calculatedEC = false;
		S  = NULL;    unitS = Atomic;   scaleS = lin;   needS = false;   calculatedS = false;
		ST = NULL;   unitST = Atomic;  scaleST = lin;  needST = false;  calculatedST = false;
		SC = NULL;   unitSC = Atomic;  scaleSC = lin;  needSC = false;  calculatedSC = false;
		M  = NULL;    unitM = Atomic;   scaleM = lin;   needM = false;   calculatedM = false;
		MT = NULL;   unitMT = Atomic;  scaleMT = lin;  needMT = false;  calculatedMT = false;
		MC = NULL;   unitMC = Atomic;  scaleMC = lin;  needMC = false;  calculatedMC = false;

		QEoS = NULL; unitQEoS = Atomic; scaleQEoS = lin; needQEoS = false;  calculatedQEoS = false; 
		TEoS = NULL;									 needTEoS = false;
	}

void FTTFQEmodel::setTolerance(const Double _eps) { 
	eps = _eps; 
	fttf.setTolerance(eps);
	phi.setTolerance(eps);
	psi.setTolerance(eps);
	coldPhi.setTolerance(eps);
	coldPsi.setTolerance(eps);
}

void FTTFQEmodel::setParameters(const VolVec &_V, const TempVec &_T) {
	if (V != NULL) delete V;
	if (T != NULL) delete T;
	V = new VolVec(_V);
	T = new TempVec(_T);
	transformVtoC();
	transformVtoD();

	needV = true;
	needT = true;

	calculatedP  = false; calculatedP0  = false; calculatedDP  = false; 
	calculatedPT = false; calculatedPT0 = false; calculatedDPT = false; 
	calculatedPC = false; calculatedPC0 = false; calculatedDPC = false; 
	calculatedE  = false; calculatedE0  = false; calculatedDE  = false; 
	calculatedET = false; calculatedET0 = false; calculatedDET = false; 
	calculatedEC = false; calculatedEC0 = false; calculatedDEC = false; 
	calculatedS  = false; calculatedS0  = false; calculatedDS  = false; 
	calculatedST = false; calculatedST0 = false; calculatedDST = false; 
	calculatedSC = false; calculatedSC0 = false; calculatedDSC = false; 
	calculatedM  = false; calculatedM0  = false; calculatedDM  = false; 
	calculatedMT = false; calculatedMT0 = false; calculatedDMT = false; 
	calculatedMC = false; calculatedMC0 = false; calculatedDMC = false; 

	if (P  != NULL) delete P;  if (P0  != NULL) delete P0;  if (DP  != NULL) delete DP;
	if (PT != NULL) delete PT; if (PT0 != NULL) delete PT0; if (DPT != NULL) delete DPT;
	if (PC != NULL) delete PC; if (PC0 != NULL) delete PC0; if (DPC != NULL) delete DPC;
	if (E  != NULL) delete E;  if (E0  != NULL) delete E0;  if (DE  != NULL) delete DE;
	if (ET != NULL) delete ET; if (ET0 != NULL) delete ET0; if (DET != NULL) delete DET;
	if (EC != NULL) delete EC; if (EC0 != NULL) delete EC0; if (DEC != NULL) delete DEC;
	if (S  != NULL) delete S;  if (S0  != NULL) delete S0;  if (DS  != NULL) delete DS;
	if (ST != NULL) delete ST; if (ST0 != NULL) delete ST0; if (DST != NULL) delete DST;
	if (SC != NULL) delete SC; if (SC0 != NULL) delete SC0; if (DSC != NULL) delete DSC;
	if (M  != NULL) delete M;  if (M0  != NULL) delete M0;  if (DM  != NULL) delete DM;
	if (MT != NULL) delete MT; if (MT0 != NULL) delete MT0; if (DMT != NULL) delete DMT;
	if (MC != NULL) delete MC; if (MC0 != NULL) delete MC0; if (DMC != NULL) delete DMC;

	calculatedQEoS = false; if (QEoS != NULL) delete QEoS;
	calculatedTEoS = false; if (TEoS != NULL) delete TEoS;

	fttf.setParameters(_V, _T);
}

void FTTFQEmodel::setParameters(const DensVec &_D, const TempVec &_T) {
	if (D != NULL) delete D;
	if (T != NULL) delete T;
	D = new DensVec(_D);
	T = new TempVec(_T);
	transformDtoV();
	transformVtoC();

	needD = true;
	needT = true;

	calculatedP  = false; calculatedP0  = false; calculatedDP  = false; 
	calculatedPT = false; calculatedPT0 = false; calculatedDPT = false; 
	calculatedPC = false; calculatedPC0 = false; calculatedDPC = false; 
	calculatedE  = false; calculatedE0  = false; calculatedDE  = false; 
	calculatedET = false; calculatedET0 = false; calculatedDET = false; 
	calculatedEC = false; calculatedEC0 = false; calculatedDEC = false; 
	calculatedS  = false; calculatedS0  = false; calculatedDS  = false; 
	calculatedST = false; calculatedST0 = false; calculatedDST = false; 
	calculatedSC = false; calculatedSC0 = false; calculatedDSC = false; 
	calculatedM  = false; calculatedM0  = false; calculatedDM  = false; 
	calculatedMT = false; calculatedMT0 = false; calculatedDMT = false; 
	calculatedMC = false; calculatedMC0 = false; calculatedDMC = false; 

	if (P  != NULL) delete P;  if (P0  != NULL) delete P0;  if (DP  != NULL) delete DP;
	if (PT != NULL) delete PT; if (PT0 != NULL) delete PT0; if (DPT != NULL) delete DPT;
	if (PC != NULL) delete PC; if (PC0 != NULL) delete PC0; if (DPC != NULL) delete DPC;
	if (E  != NULL) delete E;  if (E0  != NULL) delete E0;  if (DE  != NULL) delete DE;
	if (ET != NULL) delete ET; if (ET0 != NULL) delete ET0; if (DET != NULL) delete DET;
	if (EC != NULL) delete EC; if (EC0 != NULL) delete EC0; if (DEC != NULL) delete DEC;
	if (S  != NULL) delete S;  if (S0  != NULL) delete S0;  if (DS  != NULL) delete DS;
	if (ST != NULL) delete ST; if (ST0 != NULL) delete ST0; if (DST != NULL) delete DST;
	if (SC != NULL) delete SC; if (SC0 != NULL) delete SC0; if (DSC != NULL) delete DSC;
	if (M  != NULL) delete M;  if (M0  != NULL) delete M0;  if (DM  != NULL) delete DM;
	if (MT != NULL) delete MT; if (MT0 != NULL) delete MT0; if (DMT != NULL) delete DMT;
	if (MC != NULL) delete MC; if (MC0 != NULL) delete MC0; if (DMC != NULL) delete DMC;

	calculatedQEoS = false; if (QEoS != NULL) delete QEoS;
	calculatedTEoS = false; if (TEoS != NULL) delete TEoS;

	fttf.setParameters(_D, _T);
}

void FTTFQEmodel::setParameters(const ConVec &_C, const TempVec &_T) {
	if (C != NULL) delete C;
	if (T != NULL) delete T;
	C = new ConVec(_C);
	T = new TempVec(_T);
	transformCtoV();
	transformVtoD();

	needC = true;
	needT = true;

	calculatedP  = false; calculatedP0  = false; calculatedDP  = false; 
	calculatedPT = false; calculatedPT0 = false; calculatedDPT = false; 
	calculatedPC = false; calculatedPC0 = false; calculatedDPC = false; 
	calculatedE  = false; calculatedE0  = false; calculatedDE  = false; 
	calculatedET = false; calculatedET0 = false; calculatedDET = false; 
	calculatedEC = false; calculatedEC0 = false; calculatedDEC = false; 
	calculatedS  = false; calculatedS0  = false; calculatedDS  = false; 
	calculatedST = false; calculatedST0 = false; calculatedDST = false; 
	calculatedSC = false; calculatedSC0 = false; calculatedDSC = false; 
	calculatedM  = false; calculatedM0  = false; calculatedDM  = false; 
	calculatedMT = false; calculatedMT0 = false; calculatedDMT = false; 
	calculatedMC = false; calculatedMC0 = false; calculatedDMC = false; 

	if (P  != NULL) delete P;  if (P0  != NULL) delete P0;  if (DP  != NULL) delete DP;
	if (PT != NULL) delete PT; if (PT0 != NULL) delete PT0; if (DPT != NULL) delete DPT;
	if (PC != NULL) delete PC; if (PC0 != NULL) delete PC0; if (DPC != NULL) delete DPC;
	if (E  != NULL) delete E;  if (E0  != NULL) delete E0;  if (DE  != NULL) delete DE;
	if (ET != NULL) delete ET; if (ET0 != NULL) delete ET0; if (DET != NULL) delete DET;
	if (EC != NULL) delete EC; if (EC0 != NULL) delete EC0; if (DEC != NULL) delete DEC;
	if (S  != NULL) delete S;  if (S0  != NULL) delete S0;  if (DS  != NULL) delete DS;
	if (ST != NULL) delete ST; if (ST0 != NULL) delete ST0; if (DST != NULL) delete DST;
	if (SC != NULL) delete SC; if (SC0 != NULL) delete SC0; if (DSC != NULL) delete DSC;
	if (M  != NULL) delete M;  if (M0  != NULL) delete M0;  if (DM  != NULL) delete DM;
	if (MT != NULL) delete MT; if (MT0 != NULL) delete MT0; if (DMT != NULL) delete DMT;
	if (MC != NULL) delete MC; if (MC0 != NULL) delete MC0; if (DMC != NULL) delete DMC;

	calculatedQEoS = false; if (QEoS != NULL) delete QEoS;
	calculatedTEoS = false; if (TEoS != NULL) delete TEoS;

	fttf.setParameters(_C, _T);
}

void FTTFQEmodel::transformDtoV() {
	Int n = D->size();
	V = new VolVec(n);
	Double currentV;
	for (Int i = 0; i < n; ++i) {
		currentV = Mass/(Avogadro*1.4818e-25)/(*D)[i](gOverCmc);
		(*V)[i].setValue(currentV);
	}
}

void FTTFQEmodel::transformCtoV() {
	Int n = C->size();
	V = new VolVec(n);
	Double currentV;
	for (Int i = 0; i < n; ++i) {
		currentV = 1.0/(*C)[i]();
		(*V)[i].setValue(currentV);
	}
}

void FTTFQEmodel::transformVtoD() {
	Int n = V->size();
	D = new DensVec(n);
	Double currentD;
	for (Int i = 0; i < n; ++i) {
		currentD = Mass/(Avogadro*1.4818e-25)/(*V)[i]();
		(*D)[i].setValue(currentD, gOverCmc);
	}
}

void FTTFQEmodel::transformVtoC() {
	Int n = V->size();
	C = new ConVec(n);
	Double currentC;
	for (Int i = 0; i < n; ++i) {
		currentC = 1.0/(*V)[i]();
		(*C)[i].setValue(currentC);
	}
}

void FTTFQEmodel::prepareOutput(std::string inputString) {
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
    	if (!isalnum(*(quantity.end() - 1))) {
			unit = split(quantity, startUnit)[1];
			unit = std::string(unit.begin(), unit.end() - 1);
			quantity = split(quantity, startUnit)[0];
			if (split(unit, unitScaleDelim).size() > 1){
				scale = split(unit, unitScaleDelim)[1];
				unit = split(unit, unitScaleDelim)[0];
			}
    	}
    	if (!quantity.compare("V")) { needV = true; if (!unit.empty()) unitV = stringToUnit(unit); if (!scale.empty()) scaleV = stringToScale(scale); }
    	if (!quantity.compare("D")) { needD = true; if (!unit.empty()) unitD = stringToUnit(unit); if (!scale.empty()) scaleD = stringToScale(scale); }
 		if (!quantity.compare("C")) { needC = true; if (!unit.empty()) unitC = stringToUnit(unit); if (!scale.empty()) scaleC = stringToScale(scale); }
    	if (!quantity.compare("T")) { needT = true; if (!unit.empty()) unitT = stringToUnit(unit); if (!scale.empty()) scaleT = stringToScale(scale); }

    	if (!quantity.compare("P"))  { needP  = true; if (!unit.empty()) unitP  = stringToUnit(unit); if (!scale.empty()) scaleP =  stringToScale(scale); }
		if (!quantity.compare("PT")) { needPT = true; if (!unit.empty()) unitPT = stringToUnit(unit); if (!scale.empty()) scalePT = stringToScale(scale); }
		if (!quantity.compare("PC")) { needPC = true; if (!unit.empty()) unitPC = stringToUnit(unit); if (!scale.empty()) scalePC = stringToScale(scale); }
		if (!quantity.compare("E"))  { needE  = true; if (!unit.empty()) unitE  = stringToUnit(unit); if (!scale.empty()) scaleE =  stringToScale(scale); }
		if (!quantity.compare("ET")) { needET = true; if (!unit.empty()) unitET = stringToUnit(unit); if (!scale.empty()) scaleET = stringToScale(scale); }
		if (!quantity.compare("EC")) { needEC = true; if (!unit.empty()) unitEC = stringToUnit(unit); if (!scale.empty()) scaleEC = stringToScale(scale); }
		if (!quantity.compare("S"))  { needS  = true; if (!unit.empty()) unitS  = stringToUnit(unit); if (!scale.empty()) scaleS =  stringToScale(scale); }
		if (!quantity.compare("ST")) { needST = true; if (!unit.empty()) unitST = stringToUnit(unit); if (!scale.empty()) scaleST = stringToScale(scale); }
		if (!quantity.compare("SC")) { needSC = true; if (!unit.empty()) unitSC = stringToUnit(unit); if (!scale.empty()) scaleSC = stringToScale(scale); }
		if (!quantity.compare("M"))  { needM  = true; if (!unit.empty()) unitM  = stringToUnit(unit); if (!scale.empty()) scaleM =  stringToScale(scale); }
		if (!quantity.compare("MT")) { needMT = true; if (!unit.empty()) unitMT = stringToUnit(unit); if (!scale.empty()) scaleMT = stringToScale(scale); }
		if (!quantity.compare("MC")) { needMC = true; if (!unit.empty()) unitMC = stringToUnit(unit); if (!scale.empty()) scaleMC = stringToScale(scale); }

		if (!quantity.compare("P0"))  { needP0  = true; if (!unit.empty()) unitP0  = stringToUnit(unit); if (!scale.empty()) scaleP0 =  stringToScale(scale); }
		if (!quantity.compare("PT0")) { needPT0 = true; if (!unit.empty()) unitPT0 = stringToUnit(unit); if (!scale.empty()) scalePT0 = stringToScale(scale); }
		if (!quantity.compare("PC0")) { needPC0 = true; if (!unit.empty()) unitPC0 = stringToUnit(unit); if (!scale.empty()) scalePC0 = stringToScale(scale); }
		if (!quantity.compare("E0"))  { needE0  = true; if (!unit.empty()) unitE0  = stringToUnit(unit); if (!scale.empty()) scaleE0 =  stringToScale(scale); }
		if (!quantity.compare("ET0")) { needET0 = true; if (!unit.empty()) unitET0 = stringToUnit(unit); if (!scale.empty()) scaleET0 = stringToScale(scale); }
		if (!quantity.compare("EC0")) { needEC0 = true; if (!unit.empty()) unitEC0 = stringToUnit(unit); if (!scale.empty()) scaleEC0 = stringToScale(scale); }
		if (!quantity.compare("S0"))  { needS0  = true; if (!unit.empty()) unitS0  = stringToUnit(unit); if (!scale.empty()) scaleS0 =  stringToScale(scale); }
		if (!quantity.compare("ST0")) { needST0 = true; if (!unit.empty()) unitST0 = stringToUnit(unit); if (!scale.empty()) scaleST0 = stringToScale(scale); }
		if (!quantity.compare("SC0")) { needSC0 = true; if (!unit.empty()) unitSC0 = stringToUnit(unit); if (!scale.empty()) scaleSC0 = stringToScale(scale); }
		if (!quantity.compare("M0"))  { needM0  = true; if (!unit.empty()) unitM0  = stringToUnit(unit); if (!scale.empty()) scaleM0 =  stringToScale(scale); }
		if (!quantity.compare("MT0")) { needMT0 = true; if (!unit.empty()) unitMT0 = stringToUnit(unit); if (!scale.empty()) scaleMT0 = stringToScale(scale); }
		if (!quantity.compare("MC0")) { needMC0 = true; if (!unit.empty()) unitMC0 = stringToUnit(unit); if (!scale.empty()) scaleMC0 = stringToScale(scale); }

		if (!quantity.compare("DP"))  { needDP  = true; if (!unit.empty()) unitDP  = stringToUnit(unit); if (!scale.empty()) scaleDP =  stringToScale(scale); }
		if (!quantity.compare("DPT")) { needDPT = true; if (!unit.empty()) unitDPT = stringToUnit(unit); if (!scale.empty()) scaleDPT = stringToScale(scale); }
		if (!quantity.compare("DPC")) { needDPC = true; if (!unit.empty()) unitDPC = stringToUnit(unit); if (!scale.empty()) scaleDPC = stringToScale(scale); }
		if (!quantity.compare("DE"))  { needDE  = true; if (!unit.empty()) unitDE  = stringToUnit(unit); if (!scale.empty()) scaleDE =  stringToScale(scale); }
		if (!quantity.compare("DET")) { needDET = true; if (!unit.empty()) unitDET = stringToUnit(unit); if (!scale.empty()) scaleDET = stringToScale(scale); }
		if (!quantity.compare("DEC")) { needDEC = true; if (!unit.empty()) unitDEC = stringToUnit(unit); if (!scale.empty()) scaleDEC = stringToScale(scale); }
		if (!quantity.compare("DS"))  { needDS  = true; if (!unit.empty()) unitDS  = stringToUnit(unit); if (!scale.empty()) scaleDS =  stringToScale(scale); }
		if (!quantity.compare("DST")) { needDST = true; if (!unit.empty()) unitDST = stringToUnit(unit); if (!scale.empty()) scaleDST = stringToScale(scale); }
		if (!quantity.compare("DSC")) { needDSC = true; if (!unit.empty()) unitDSC = stringToUnit(unit); if (!scale.empty()) scaleDSC = stringToScale(scale); }
		if (!quantity.compare("DM"))  { needDM  = true; if (!unit.empty()) unitDM  = stringToUnit(unit); if (!scale.empty()) scaleDM =  stringToScale(scale); }
		if (!quantity.compare("DMT")) { needDMT = true; if (!unit.empty()) unitDMT = stringToUnit(unit); if (!scale.empty()) scaleDMT = stringToScale(scale); }
		if (!quantity.compare("DMC")) { needDMC = true; if (!unit.empty()) unitDMC = stringToUnit(unit); if (!scale.empty()) scaleDMC = stringToScale(scale); }

		if (!quantity.compare("QEoS")) { needQEoS = true; if (!unit.empty()) unitQEoS = stringToUnit(unit); if (!scale.empty()) scaleQEoS = stringToScale(scale); }
		if (!quantity.compare("TEoS")) { needTEoS = true; }
    }
    calculateAll();
}

void FTTFQEmodel::printString(std::ofstream &out, std::string data) {
	out.setf(std::ios::left);
    out.width(20);
    out.fill(' ');
    out << data;
}

void FTTFQEmodel::printValue(std::ofstream &out, Double value, Int digitsToPrint) {
	out.setf(std::ios::scientific | std::ios::left | std::ios::showpos);
    out.width(20);
    out.precision(digitsToPrint);
    out.fill(' ');
    out << value;
}

void FTTFQEmodel::printOutput(const char* filename) {
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

	if (needP0)   { s = "P0[";   s += unitToString(unitP0);   s += ","; s += scaleToString(scaleP0);   s += "]"; printString(out, s); }
	if (needPT0)  { s = "PT0[";  s += unitToString(unitPT0);  s += ","; s += scaleToString(scalePT0);  s += "]"; printString(out, s); }
	if (needPC0)  { s = "PC0[";  s += unitToString(unitPC0);  s += ","; s += scaleToString(scalePC0);  s += "]"; printString(out, s); }
	if (needE0)   { s = "E0[";   s += unitToString(unitE0);   s += ","; s += scaleToString(scaleE0);   s += "]"; printString(out, s); }
	if (needET0)  { s = "ET0[";  s += unitToString(unitET0);  s += ","; s += scaleToString(scaleET0);  s += "]"; printString(out, s); }
	if (needEC0)  { s = "EC0[";  s += unitToString(unitEC0);  s += ","; s += scaleToString(scaleEC0);  s += "]"; printString(out, s); }
	if (needS0)   { s = "S0[";   s += unitToString(unitS0);   s += ","; s += scaleToString(scaleS0);   s += "]"; printString(out, s); }
	if (needST0)  { s = "ST0[";  s += unitToString(unitST0);  s += ","; s += scaleToString(scaleST0);  s += "]"; printString(out, s); }
	if (needSC0)  { s = "SC0[";  s += unitToString(unitSC0);  s += ","; s += scaleToString(scaleSC0);  s += "]"; printString(out, s); }
	if (needM0)   { s = "M0[";   s += unitToString(unitM0);   s += ","; s += scaleToString(scaleM0);   s += "]"; printString(out, s); }
	if (needMT0)  { s = "MT0[";  s += unitToString(unitMT0);  s += ","; s += scaleToString(scaleMT0);  s += "]"; printString(out, s); }
	if (needMC0)  { s = "MC0[";  s += unitToString(unitMC0);  s += ","; s += scaleToString(scaleMC0);  s += "]"; printString(out, s); }

	if (needDP)   { s = "DP[";   s += unitToString(unitDP);   s += ","; s += scaleToString(scaleDP);   s += "]"; printString(out, s); }
	if (needDPT)  { s = "DPT[";  s += unitToString(unitDPT);  s += ","; s += scaleToString(scaleDPT);  s += "]"; printString(out, s); }
	if (needDPC)  { s = "DPC[";  s += unitToString(unitDPC);  s += ","; s += scaleToString(scaleDPC);  s += "]"; printString(out, s); }
	if (needDE)   { s = "DE[";   s += unitToString(unitDE);   s += ","; s += scaleToString(scaleDE);   s += "]"; printString(out, s); }
	if (needDET)  { s = "DET[";  s += unitToString(unitDET);  s += ","; s += scaleToString(scaleDET);  s += "]"; printString(out, s); }
	if (needDEC)  { s = "DEC[";  s += unitToString(unitDEC);  s += ","; s += scaleToString(scaleDEC);  s += "]"; printString(out, s); }
	if (needDS)   { s = "DS[";   s += unitToString(unitDS);   s += ","; s += scaleToString(scaleDS);   s += "]"; printString(out, s); }
	if (needDST)  { s = "DST[";  s += unitToString(unitDST);  s += ","; s += scaleToString(scaleDST);  s += "]"; printString(out, s); }
	if (needDSC)  { s = "DSC[";  s += unitToString(unitDSC);  s += ","; s += scaleToString(scaleDSC);  s += "]"; printString(out, s); }
	if (needDM)   { s = "DM[";   s += unitToString(unitDM);   s += ","; s += scaleToString(scaleDM);   s += "]"; printString(out, s); }
	if (needDMT)  { s = "DMT[";  s += unitToString(unitDMT);  s += ","; s += scaleToString(scaleDMT);  s += "]"; printString(out, s); }
	if (needDMC)  { s = "DMC[";  s += unitToString(unitDMC);  s += ","; s += scaleToString(scaleDMC);  s += "]"; printString(out, s); }

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

			if (needP)     printValue(out,  (*P )[v*Tsize + t](unitP,   scaleP),   digitsToPrint);
			if (needPT)    printValue(out,  (*PT)[v*Tsize + t](unitPT,  scalePT),  digitsToPrint);
			if (needPC)    printValue(out,  (*PC)[v](          unitPC,  scalePC),  digitsToPrint);
			if (needE)     printValue(out,  (*E )[v*Tsize + t](unitE,   scaleE),   digitsToPrint);
			if (needET)    printValue(out,  (*ET)[v*Tsize + t](unitET,  scaleET),  digitsToPrint);
			if (needEC)    printValue(out,  (*EC)[v](          unitEC,  scaleEC),  digitsToPrint);
			if (needS)     printValue(out,  (*S )[v*Tsize + t](unitS,   scaleS),   digitsToPrint);
			if (needST)    printValue(out,  (*ST)[v*Tsize + t](unitST,  scaleST),  digitsToPrint);
			if (needSC)    printValue(out,  (*SC)[v](          unitSC,  scaleSC),  digitsToPrint);
			if (needM)     printValue(out,  (*M )[v*Tsize + t](unitM,   scaleM ),  digitsToPrint);
			if (needMT)    printValue(out,  (*MT)[v*Tsize + t](unitMT,  scaleMT),  digitsToPrint);
			if (needMC)    printValue(out,  (*MC)[v](          unitMC,  scaleMC),  digitsToPrint);

			if (needP0)    printValue(out, (*P0 )[v*Tsize + t](unitP0,  scaleP0),  digitsToPrint);
			if (needPT0)   printValue(out, (*PT0)[v*Tsize + t](unitPT0, scalePT0), digitsToPrint);
			if (needPC0)   printValue(out, (*PC0)[v](          unitPC0, scalePC0), digitsToPrint);
			if (needE0)    printValue(out, (*E0 )[v*Tsize + t](unitE0,  scaleE0),  digitsToPrint);
			if (needET0)   printValue(out, (*ET0)[v*Tsize + t](unitET0, scaleET0), digitsToPrint);
			if (needEC0)   printValue(out, (*EC0)[v](          unitEC0, scaleEC0), digitsToPrint);
			if (needS0)    printValue(out, (*S0 )[v*Tsize + t](unitS0,  scaleS0),  digitsToPrint);
			if (needST0)   printValue(out, (*ST0)[v*Tsize + t](unitST0, scaleST0), digitsToPrint);
			if (needSC0)   printValue(out, (*SC0)[v](          unitSC0, scaleSC0), digitsToPrint);
			if (needM0)    printValue(out, (*M0 )[v*Tsize + t](unitM0,  scaleM0 ), digitsToPrint);
			if (needMT0)   printValue(out, (*MT0)[v*Tsize + t](unitMT0, scaleMT0), digitsToPrint);
			if (needMC0)   printValue(out, (*MC0)[v](          unitMC0, scaleMC0), digitsToPrint);

			if (needDP)    printValue(out, (*DP )[v*Tsize + t](unitDP,  scaleDP),  digitsToPrint);
			if (needDPT)   printValue(out, (*DPT)[v*Tsize + t](unitDPT, scaleDPT), digitsToPrint);
			if (needDPC)   printValue(out, (*DPC)[v](          unitDPC, scaleDPC), digitsToPrint);
			if (needDE)    printValue(out, (*DE )[v*Tsize + t](unitDE,  scaleDE),  digitsToPrint);
			if (needDET)   printValue(out, (*DET)[v*Tsize + t](unitDET, scaleDET), digitsToPrint);
			if (needDEC)   printValue(out, (*DEC)[v](          unitDEC, scaleDEC), digitsToPrint);
			if (needDS)    printValue(out, (*DS )[v*Tsize + t](unitDS,  scaleDS),  digitsToPrint);
			if (needDST)   printValue(out, (*DST)[v*Tsize + t](unitDST, scaleDST), digitsToPrint);
			if (needDSC)   printValue(out, (*DSC)[v](          unitDSC, scaleDSC), digitsToPrint);
			if (needDM)    printValue(out, (*DM )[v*Tsize + t](unitDM,  scaleDM ), digitsToPrint);
			if (needDMT)   printValue(out, (*DMT)[v*Tsize + t](unitDMT, scaleDMT), digitsToPrint);
			if (needDMC)   printValue(out, (*DMC)[v](          unitDMC, scaleDMC), digitsToPrint);

			if (needQEoS) printValue(out, (*QEoS)[v*Tsize + t](unitQEoS, scaleQEoS), digitsToPrint);
			if (needTEoS) printValue(out, (*TEoS)[v*Tsize + t], digitsToPrint);
			out << std::endl;
	    }
    }
}

PressVec&   FTTFQEmodel::getP()  {if (!calculatedP)  calculateP();  return *P;  }
PressVec&   FTTFQEmodel::getPT() {if (!calculatedPT) calculatePT(); return *PT; }
PressVec&   FTTFQEmodel::getPC() {if (!calculatedPC) calculatePC(); return *PC; }
EnergyVec&  FTTFQEmodel::getE()  {if (!calculatedE)  calculateE();  return *E;  }
EnergyVec&  FTTFQEmodel::getET() {if (!calculatedET) calculateET(); return *ET; }
EnergyVec&  FTTFQEmodel::getEC() {if (!calculatedEC) calculateEC(); return *EC; }
EntropyVec& FTTFQEmodel::getS()  {if (!calculatedS)  calculateS();  return *S;  }
EntropyVec& FTTFQEmodel::getST() {if (!calculatedST) calculateST(); return *ST; }
EntropyVec& FTTFQEmodel::getSC() {if (!calculatedSC) calculateSC(); return *SC; }
ChemPotVec& FTTFQEmodel::getM()  {if (!calculatedM)  calculateM();  return *M;  }
ChemPotVec& FTTFQEmodel::getMT() {if (!calculatedMT) calculateMT(); return *MT; }
ChemPotVec& FTTFQEmodel::getMC() {if (!calculatedMC) calculateMC(); return *MC; }

PressVec&   FTTFQEmodel::getDP()  {if (!calculatedDP)  calculateDP();  return *DP;  }
PressVec&   FTTFQEmodel::getDPT() {if (!calculatedDPT) calculateDPT(); return *DPT; }
PressVec&   FTTFQEmodel::getDPC() {if (!calculatedDPC) calculateDPC(); return *DPC; }
EnergyVec&  FTTFQEmodel::getDE()  {if (!calculatedDE)  calculateDE();  return *DE;  }
EnergyVec&  FTTFQEmodel::getDET() {if (!calculatedDET) calculateDET(); return *DET; }
EnergyVec&  FTTFQEmodel::getDEC() {if (!calculatedDEC) calculateDEC(); return *DEC; }
EntropyVec& FTTFQEmodel::getDS()  {if (!calculatedDS)  calculateDS();  return *DS;  }
EntropyVec& FTTFQEmodel::getDST() {if (!calculatedDST) calculateDST(); return *DST; }
EntropyVec& FTTFQEmodel::getDSC() {if (!calculatedDSC) calculateDSC(); return *DSC; }
ChemPotVec& FTTFQEmodel::getDM()  {if (!calculatedDM)  calculateDM();  return *DM;  }
ChemPotVec& FTTFQEmodel::getDMT() {if (!calculatedDMT) calculateDMT(); return *DMT; }
ChemPotVec& FTTFQEmodel::getDMC() {if (!calculatedDMC) calculateDMC(); return *DMC; }

PressVec&   FTTFQEmodel::getP0()  { return fttf.getP();  }
PressVec&   FTTFQEmodel::getPT0() { return fttf.getPT(); }
PressVec&   FTTFQEmodel::getPC0() { return fttf.getPC(); }
EnergyVec&  FTTFQEmodel::getE0()  { return fttf.getE();  }
EnergyVec&  FTTFQEmodel::getET0() { return fttf.getET(); }
EnergyVec&  FTTFQEmodel::getEC0() { return fttf.getEC(); }
EntropyVec& FTTFQEmodel::getS0()  { return fttf.getS();  }
EntropyVec& FTTFQEmodel::getST0() { return fttf.getST(); }
EntropyVec& FTTFQEmodel::getSC0() { return fttf.getSC(); }
ChemPotVec& FTTFQEmodel::getM0()  { return fttf.getM();  }
ChemPotVec& FTTFQEmodel::getMT0() { return fttf.getMT(); }
ChemPotVec& FTTFQEmodel::getMC0() { return fttf.getMC(); }

void FTTFQEmodel::calculateAll() {
	Int Vsize = V->size();
	Int Tsize = T->size();
	// prepare necessary data in FTTF model
	std::string FTTFinput = "";

	if ( needP0 ||  needP) FTTFinput += "P ";
	if (needPT0 || needPT) FTTFinput += "PT ";
	if (needPC0 || needPC) FTTFinput += "PC ";
	if ( needE0 ||  needE) FTTFinput += "E ";
	if (needET0 || needET) FTTFinput += "ET ";
	if (needEC0 || needEC) FTTFinput += "EC ";
	if ( needS0 ||  needS) FTTFinput += "S ";
	if (needST0 || needST) FTTFinput += "ST ";
	if (needSC0 || needSC) FTTFinput += "SC ";
	if ( needM0 ||  needM) FTTFinput += "M ";
	if (needMT0 || needMT) FTTFinput += "MT ";
	if (needMC0 || needMC) FTTFinput += "MC ";
	
	if (!FTTFinput.empty()) fttf.prepareOutput(FTTFinput);

	if ( needP0 ||  needP)   P0 = &fttf.getP();
	if (needPT0 || needPT) { PT0 = &fttf.getPT(); P0 = &fttf.getP(); PC0 = &fttf.getPC(); }
	if (needPC0 || needPC)   PC0 = &fttf.getPC();
	if ( needE0 ||  needE)   E0 = &fttf.getE();
	if (needET0 || needET) { ET0 = &fttf.getET(); E0 = &fttf.getE(); EC0 = &fttf.getEC(); } 
	if (needEC0 || needEC)   EC0 = &fttf.getEC();
	if ( needS0 ||  needS)   S0 = &fttf.getS();
	if (needST0 || needST) { ST0 = &fttf.getST(); S0 = &fttf.getS(); SC0 = &fttf.getSC(); }
	if (needSC0 || needSC)   SC0 = &fttf.getSC();
	if ( needM0 ||  needM)   M0 = &fttf.getM();
	if (needMT0 || needMT) { MT0 = &fttf.getMT(); M0 = &fttf.getM(); MC0 = &fttf.getMC(); }
	if (needMC0 || needMC)   MC0 = &fttf.getMC();

    // allocate memory for the necessary arrays
	if (  needP) { 
		if (DP == NULL) DP = new PressVec(Vsize*Tsize); 
		if (P0 == NULL) P0 = &fttf.getP(); 
		P = new PressVec(Vsize*Tsize); 
	}
	if ( needPT) { 
		if (DP  == NULL) DP  = new PressVec(Vsize*Tsize); 
		if (P0  == NULL) P0  = &fttf.getP();  
		if (DPC == NULL) DPC = new PressVec(Vsize*Tsize); 
		if (PC0 == NULL) PC0 = &fttf.getPC();
		PT = new PressVec(Vsize*Tsize); 
	}
	if ( needPC) {
		if (DPC == NULL) DPC = new PressVec(Vsize); 
		if (PC0 == NULL) PC0 = &fttf.getPC(); 
		PC = new PressVec(Vsize);
	}
	if (  needE) {
		if (DE == NULL) DE = new EnergyVec(Vsize*Tsize); 
		if (E0 == NULL) E0 = &fttf.getE(); 
		E = new EnergyVec(Vsize*Tsize);	
	} 
	if ( needET) { 
		if (DE  == NULL) DE  = new EnergyVec(Vsize*Tsize); 
		if (E0  == NULL) E0  = &fttf.getE();  
		if (DEC == NULL) DEC = new EnergyVec(Vsize*Tsize); 
		if (EC0 == NULL) EC0 = &fttf.getEC();
		ET = new EnergyVec(Vsize*Tsize); 
	}
	if ( needEC) { 
		if (DEC == NULL) DEC = new EnergyVec(Vsize); 
		if (EC0 == NULL) EC0 = &fttf.getEC(); 
		EC = new EnergyVec(Vsize); 
	}
	if (  needS) { 
		if (DS == NULL) DS = new EntropyVec(Vsize*Tsize); 
		if (S0 == NULL) S0 = &fttf.getS(); 
		S = new EntropyVec(Vsize*Tsize); 
	}
	if ( needST) { 
		if (DS  == NULL) DS  = new EntropyVec(Vsize*Tsize); 
		if (S0  == NULL) S0  = &fttf.getS();  
		if (DSC == NULL) DSC = new EntropyVec(Vsize*Tsize); 
		if (SC0 == NULL) SC0 = &fttf.getSC();
		ST = new EntropyVec(Vsize*Tsize); 
	}
	if ( needSC) { 
		if (DSC == NULL) DSC = new EntropyVec(Vsize); 
		if (SC0 == NULL) SC0 = &fttf.getSC(); 
		SC = new EntropyVec(Vsize); 
	}
	if (  needM) { 
		if (DM == NULL) DM = new ChemPotVec(Vsize*Tsize); 
		if (M0 == NULL) M0 = &fttf.getM(); 
		M = new ChemPotVec(Vsize*Tsize); 
	}
	if ( needMT) { 
		if (DM  == NULL) DM  = new ChemPotVec(Vsize*Tsize); 
		if (M0  == NULL) M0  = &fttf.getM();  
		if (DMC == NULL) DMC = new ChemPotVec(Vsize*Tsize); 
		if (MC0 == NULL) MC0 = &fttf.getMC();
		MT = new ChemPotVec(Vsize*Tsize); 
	}
	if ( needMC) { 
		if (DMC == NULL) DMC = new ChemPotVec(Vsize); 
		if (MC0 == NULL) MC0 = &fttf.getMC(); 
		MC = new ChemPotVec(Vsize); 
	}

	if ( needDP) DP = new PressVec(Vsize*Tsize);
	if (needDPT) { 
		if (DP  == NULL) DP  = new PressVec(Vsize*Tsize); 
		if (DPC == NULL) DPC = new PressVec(Vsize); 
		DPT = new PressVec(Vsize*Tsize); 
	}
	if (needDPC) DPC = new PressVec(Vsize);
	if ( needDE) DE = new EnergyVec(Vsize*Tsize);
	if (needDET) { 
		if (DE  == NULL) DE  = new EnergyVec(Vsize*Tsize); 
		if (DEC == NULL) DEC = new EnergyVec(Vsize); 
		DET = new EnergyVec(Vsize*Tsize); 
	}
	if (needDEC) DEC = new EnergyVec(Vsize);
	if ( needDS) DS = new EntropyVec(Vsize*Tsize);
	if (needDST) { 
		if (DS  == NULL) DS  = new EntropyVec(Vsize*Tsize); 
		if (DSC == NULL) DSC = new EntropyVec(Vsize); 
		DST = new EntropyVec(Vsize*Tsize); 
	}
	if (needDSC) DSC = new EntropyVec(Vsize);
	if ( needDM) DM = new ChemPotVec(Vsize*Tsize);
	if (needDMT) { 
		if (DM  == NULL) DM  = new ChemPotVec(Vsize*Tsize); 
		if (DMC == NULL) DMC = new ChemPotVec(Vsize); 
		DMT = new ChemPotVec(Vsize*Tsize); 
	}
	if (needDMC) DMC = new ChemPotVec(Vsize);

	if (needQEoS) { 
		if (P  == NULL) P  = new PressVec(Vsize*Tsize); 
		if (E  == NULL) E  = new EnergyVec(Vsize*Tsize); 
		if (DP == NULL) DP  = new PressVec(Vsize*Tsize); 
		if (DE == NULL) DE  = new EnergyVec(Vsize*Tsize); 
		QEoS = new EnergyVec(Vsize*Tsize); 
	}
	if (needTEoS) { 
		if (P  == NULL) P  = new PressVec(Vsize*Tsize); 
		if (DP == NULL) DP  = new PressVec(Vsize*Tsize); 
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

		calculatedDPC = false; 
		calculatedDEC = false; 
		calculatedDSC = false; 
		calculatedDMC = false; 
		calculatedPC  = false; 
		calculatedEC  = false; 
		calculatedSC  = false; 
		calculatedMC  = false;

		coldPhi.setParameters((*V)[v](), coldT(), Z);
		coldPsi.setParameters((*V)[v](), coldT(), Z);

		if (needDPC) calculateDPC(v);
		if (needDEC) calculateDEC(v);
		if (needDSC) calculateDSC(v);
		if (needDMC) calculateDMC(v);
		if ( needPC)  calculatePC(v);
		if ( needEC)  calculateEC(v);
		if ( needSC)  calculateSC(v);
		if ( needMC)  calculateMC(v);
		for (Int t = 0; t < Tsize; ++t) {
			
			calculatedDP  = false; 
			calculatedDPT = false; 
			calculatedDE  = false; 
			calculatedDET = false; 
			calculatedDS  = false; 
			calculatedDST = false; 
			calculatedDM  = false; 
			calculatedDMT = false; 

			calculatedP   = false; 
			calculatedPT  = false; 
			calculatedE   = false; 
			calculatedET  = false; 
			calculatedS   = false; 
			calculatedST  = false; 
			calculatedM   = false; 
			calculatedMT  = false; 

			calculatedQEoS = false; 
			calculatedTEoS = false; 

			phi.setParameters((*V)[v](), (*T)[t](), Z);
			psi.setParameters((*V)[v](), (*T)[t](), Z);

			if ( needDP)  calculateDP(v, t); 
			if (needDPT) calculateDPT(v, t); 
			if ( needDE)  calculateDE(v, t); 
			if (needDET) calculateDET(v, t); 
			if ( needDS)  calculateDS(v, t); 
			if (needDST) calculateDST(v, t); 
			if ( needDM)  calculateDM(v, t); 
			if (needDMT) calculateDMT(v, t); 

			if (  needP)   calculateP(v, t); 
			if ( needPT)  calculatePT(v, t); 
			if (  needE)   calculateE(v, t); 
			if ( needET)  calculateET(v, t); 
			if (  needS)   calculateS(v, t); 
			if ( needST)  calculateST(v, t); 
			if (  needM)   calculateM(v, t); 
			if ( needMT)  calculateMT(v, t); 

			if (needQEoS) calculateQEoS(v, t); 
			if (needTEoS) calculateTEoS(v, t); 
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// corrections
//////////////////////////////////////////////////////////////////////////////////////////////////
void FTTFQEmodel::calculateDP() {
	Int Vsize = V->size();
	Int Tsize = T->size();
	Int Psize = Vsize*Tsize;
	DP = new PressVec(Psize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			phi.setParameters((*V)[v](), (*T)[t](), Z);
			psi.setParameters((*V)[v](), (*T)[t](), Z);
			calculateDP(v, t);
		}
	}
	calculatedDP  = true;
}
void FTTFQEmodel::calculateDPT() {
	if (!calculatedDP) calculateDP();
	if (!calculatedDPC) calculateDPC();
	Int DPTsize = DP->size();
	Int Vsize = V->size();
	Int Tsize = T->size();
	DPT = new PressVec(DPTsize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			calculateDPT(v, t);
		}
	}
	calculatedDPT = true;
}
void FTTFQEmodel::calculateDPC() {
	Int Vsize = V->size();
	DPC = new PressVec(Vsize);
	for (Int v = 0; v < Vsize; ++v) {
		coldPhi.setParameters((*V)[v](), coldT(), Z);
		coldPsi.setParameters((*V)[v](), coldT(), Z);
		calculateDPC(v);
	}
	calculatedDPC = true;
}
void FTTFQEmodel::calculateDE() {
    Int Vsize = V->size();
	Int Tsize = T->size();
    DE = new EnergyVec(Vsize*Tsize);

	energySolver.SetOutput(energyData);
	energySolver.SetTolerance(0.0, eps/10);

    for (Int v = 0; v < Vsize; ++v) {
    	for (Int t = 0; t < Tsize; ++t) {
    		phi.setParameters((*V)[v](), (*T)[t](), Z);
    		psi.setParameters((*V)[v](), (*T)[t](), Z);
    		calculateDE(v, t);
    	}
    }
	calculatedDE  = true;
}
void FTTFQEmodel::calculateDET() {
	if (!calculatedDE)  calculateDE();
	if (!calculatedDEC) calculateDEC();
	Int DETsize = DE->size();
	Int Vsize = V->size();
	Int Tsize = T->size();
	DET = new EnergyVec(DETsize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			calculateDET(v, t);
		}
	}
	calculatedDET = true;
}
void FTTFQEmodel::calculateDEC() {
	Int Vsize = V->size();
	DEC = new EnergyVec(Vsize);

	coldEnergySolver.SetOutput(coldEnergyData);
	coldEnergySolver.SetTolerance(0.0, eps/10);

	for (Int v = 0; v < Vsize; ++v) {
		coldPhi.setParameters((*V)[v](), coldT(), Z);
		coldPsi.setParameters((*V)[v](), coldT(), Z);
		calculateDEC(v);
	}
	calculatedDEC = true;
}
void FTTFQEmodel::calculateDS() {
	Int Vsize = V->size();
	Int Tsize = T->size();
    DS = new EntropyVec(Vsize*Tsize);

	entropySolver.SetOutput(entropyData);
	entropySolver.SetTolerance(0.0, eps/10);

    for (Int v = 0; v < Vsize; ++v) {
    	for (Int t = 0; t < Tsize; ++t) {
    		phi.setParameters((*V)[v], (*T)[t], Z);
    		psi.setParameters((*V)[v], (*T)[t], Z);
			calculateDS(v, t);					
    	}
    }
	calculatedDS  = true;
}
void FTTFQEmodel::calculateDST() {
	if (!calculatedDS)  calculateDS();
	if (!calculatedDSC) calculateDSC();
	Int DSTsize = DS->size();
	Int Vsize = V->size();
	Int Tsize = T->size();
	DST = new EntropyVec(DSTsize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			calculateDST(v, t);
		}
	}
	calculatedDST = true;
}
void FTTFQEmodel::calculateDSC() {
	Int Vsize = V->size();
	DSC = new EntropyVec(Vsize);
	coldEntropySolver.SetOutput(coldEntropyData);
	coldEntropySolver.SetTolerance(0.0, eps/10);
	for (Int v = 0; v < Vsize; ++v) {
		coldPhi.setParameters((*V)[v], coldT, Z);
		coldPsi.setParameters((*V)[v], coldT, Z);
		calculateDSC(v);
	}
	calculatedDSC = true;
}
void FTTFQEmodel::calculateDM() {
	Int Vsize = V->size();
	Int Tsize = T->size();
	Int DMsize = Vsize*Tsize;
	DM = new ChemPotVec(DMsize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			phi.setParameters((*V)[v], (*T)[t], Z);
			psi.setParameters((*V)[v], (*T)[t], Z);
			calculateDM(v, t);
		}
	}
	calculatedDM  = true;
}
void FTTFQEmodel::calculateDMT() {
	if (!calculatedDM)  calculateDM();
	if (!calculatedDMC) calculateDMC();
	Int DMTsize = DM->size();
	Int Vsize = V->size();
	Int Tsize = T->size();
	DMT = new ChemPotVec(DMTsize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			calculateDMT(v, t);
		}
	}
	calculatedDMT = true;
}
void FTTFQEmodel::calculateDMC() {
	Int Vsize = V->size();
	DMC = new ChemPotVec(Vsize);
	for (Int v = 0; v < Vsize; ++v) {
		coldPhi.setParameters((*V)[v], coldT, Z);
		coldPsi.setParameters((*V)[v], coldT, Z);
		calculateDMC(v);
	}
	calculatedDMC = true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
// FTTFQE
/////////////////////////////////////////////////////////////////////////////////////////////////////
void FTTFQEmodel::calculateP() {
	if (!calculatedDP) calculateDP();
	Int Vsize = V->size();
	Int Tsize = T->size();
	Int Psize = Vsize*Tsize;
	P = new PressVec(Psize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			calculateP(v, t);
		}
	}
	calculatedP  = true;
}
void FTTFQEmodel::calculatePT() {
	if (!calculatedDP) calculateDP();
	if (!calculatedDPC) calculateDPC();
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
void FTTFQEmodel::calculatePC() {
	if (!calculatedDPC) calculateDPC();
	Int Vsize = V->size();
	PC = new PressVec(Vsize);
	for (Int v = 0; v < Vsize; ++v) {
		calculatePC(v);
	}
	calculatedPC = true;
}
void FTTFQEmodel::calculateE() {
	if (!calculatedDE) calculateDE();
	Int Vsize = V->size();
	Int Tsize = T->size();
	Int Esize = Vsize*Tsize;
	E = new EnergyVec(Esize);
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			calculateE(v, t);
		}
	}
	calculatedE = true;
}
void FTTFQEmodel::calculateET() {
	if (!calculatedDE)  calculateDE();
	if (!calculatedDEC) calculateDEC();
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
void FTTFQEmodel::calculateEC() {
	if (!calculatedDEC) calculateDEC();
	Int Vsize = V->size();
	EC = new EnergyVec(Vsize);
	Double currentEC;
	for (Int v = 0; v < Vsize; ++v) {
		calculateEC(v);
	}
	calculatedEC = true;
}
void FTTFQEmodel::calculateS() {
	if (!calculatedDS) calculateDS();
	Int Vsize = V->size();
	Int Tsize = T->size();
	Int Ssize = Vsize*Tsize;
	S = new EntropyVec(Ssize);
	Double currentS;
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			calculateS(v, t);
		}
	}
	calculatedS = true;
}
void FTTFQEmodel::calculateST() {
	if (!calculatedDS)  calculateDS();
	if (!calculatedDSC) calculateDSC();
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
void FTTFQEmodel::calculateSC() {
	if (!calculatedDSC) calculateDSC();
	Int Vsize = V->size();
	SC = new EntropyVec(Vsize);
	Double currentSC;
	for (Int v = 0; v < Vsize; ++v) {
		calculateSC(v);
	}
	calculatedSC = true;
}
void FTTFQEmodel::calculateM() {
	if (!calculatedDM) calculateDM();
	Int Vsize = V->size();
	Int Tsize = T->size();
	Int Msize = Vsize*Tsize;
	M = new ChemPotVec(Msize);
	Double currentM;
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			calculateM(v, t);
		}
	}
	calculatedM  = true;
}
void FTTFQEmodel::calculateMT() {
	if (!calculatedDM)  calculateDM();
	if (!calculatedDMC) calculateDMC();
	Int MTsize = M->size();
	Int Vsize = V->size();
	Int Tsize = T->size();
	MT = new ChemPotVec(MTsize);
	Double currentMT;
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			calculateMT(v, t);
		}
	}
	calculatedMT = true;
}
void FTTFQEmodel::calculateMC() {
	if (!calculatedDMC) calculateDMC();
	Int Vsize = V->size();
	MC = new ChemPotVec(Vsize);
	Double currentMC;
	for (Int v = 0; v < Vsize; ++v) {
		calculateMC(v);
	}
	calculatedMC = true;
}
void FTTFQEmodel::calculateQEoS() {
	if (!calculatedE) calculateE();
	if (!calculatedP) calculateP();
	Int Qsize = E->size();
	Int Vsize = V->size();
	Int Tsize = T->size();
	QEoS = new EnergyVec(Qsize);
	Double currentQ;
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			calculateQEoS(v, t);
		}
	}
	calculatedQEoS = true;
}
void FTTFQEmodel::calculateTEoS() {
	if (!calculatedP) calculateP();
	Int TEoSsize = P->size();
	Int Vsize = V->size();
	Int Tsize = T->size();
	TEoS = new DoubleVec(TEoSsize);
	Double currentTEoS;
	for (Int v = 0; v < Vsize; ++v) {
		for (Int t = 0; t < Tsize; ++t) {
			calculateTEoS(v, t);
		}
	}
	calculatedTEoS = true;
}

/////////////////////////////////////////////////////////////////////////////////////
// single point calculation functions
/////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
// corrections
//////////////////////////////////////////////////////////////////////////////////////////////////
void FTTFQEmodel::calculateDP(Int v, Int t) {
	Int Tsize = T->size();
	Double currentDP;
	currentDP = pow((*T)[t](), 2.0)/3.0/M_PI/M_PI/M_PI
				* (  FDhalf(phi.valueAt_1())*psi.valueAt_1() 
					  + Y(phi.valueAt_1()) );
	(*DP)[v*Tsize + t].setValue(currentDP);
	calculatedDP  = true;
}
void FTTFQEmodel::calculateDPT(Int v, Int t) {
	if (!calculatedDP) calculateDP(v, t);
	if (!calculatedDPC) calculateDPC(v);
	Int Tsize = T->size();
	Double currentDPT;
	currentDPT = (*DP)[v*Tsize + t]() - (*DPC)[v]();
	(*DPT)[v*Tsize + t].setValue(currentDPT);
	calculatedDPT = true;
}
void FTTFQEmodel::calculateDPC(Int v) {
	Double currentDPC;
	currentDPC = pow(2*coldT(), 2.0)/3.0/M_PI/M_PI/M_PI
				* (  FDhalf(coldPhi.valueAt_1())*coldPsi.valueAt_1() 
					  + Y(coldPhi.valueAt_1()) );
	(*DPC)[v].setValue(currentDPC);
	calculatedDPC = true;
}
void FTTFQEmodel::calculateDE(Int v, Int t) {
	Double currentDE;
    DoubleVec startDE(0.0, rhsFTTFQEenergy::dim);
    Double a, b;
    Double V1, T1;
    Double xFrom = 1.0;
    Double xTo = 0.0;
	Int Tsize = T->size();

	V1 = (*V)[v]()*Z;
	T1 = (*T)[t]()*pow(Z, -4.0/3.0);

    startDE[0] = phi.valueAt_1();
	startDE[1] = startDE[0];
	startDE[2] = psi.valueAt_1();
    startDE[3] = startDE[2];
    startDE[4] = sqrt(0.5*T1)/3.0/M_PI*psi.derivativeAt_0() + 0.269900170;

    a =   pow(2.0, 7.0/6.0)
        * pow(3.0, 2.0/3.0)
        * pow(M_PI, -5.0/3.0)
        * sqrt(T1)*pow(V1, 2.0/3.0);
    b = V1*pow(T1, 2.0)/M_PI/M_PI/M_PI;

	rhsEnergy.updateParameters(a, b);

	energySolver.Integrate(rhsEnergy, startDE, xFrom, xTo);
	currentDE = startDE[4];
	(*DE)[v*Tsize + t].setValue(currentDE*pow(Z, 5.0/3.0));
	calculatedDE  = true;
}
void FTTFQEmodel::calculateDET(Int v, Int t) {
	if (!calculatedDE)  calculateDE(v, t);
	if (!calculatedDEC) calculateDEC(v);
	Int DETsize = DE->size();
	Int Tsize = T->size();
	Double currentDET;
	currentDET = (*DE)[v*Tsize + t]() - (*DEC)[v]();
	(*DET)[v*Tsize + t].setValue(currentDET);
	calculatedDET = true;
}
void FTTFQEmodel::calculateDEC(Int v) {
	Double currentDEC;
	DoubleVec startDE(0.0, rhsFTTFQEenergy::dim);
	Double a, b;
	Double V1, T1;
	Double xFrom = 1.0;
	Double xTo = 0.0;

	T1 = coldT()*pow(Z, -4.0/3.0);
	V1 = (*V)[v]()*Z;

    startDE[0] = coldPhi.valueAt_1();
	startDE[1] = startDE[0];
	startDE[2] = coldPsi.valueAt_1();
    startDE[3] = startDE[2];
    startDE[4] = sqrt(0.5*T1)/3.0/M_PI*coldPsi.derivativeAt_0() + 0.269900170;

    a =   pow(2.0, 7.0/6.0)
        * pow(3.0, 2.0/3.0)
        * pow(M_PI, -5.0/3.0)
        * sqrt(T1)*pow(V1, 2.0/3.0);
    b = V1*pow(T1, 2.0)/M_PI/M_PI/M_PI;

	rhsColdEnergy.updateParameters(a, b);

	coldEnergySolver.Integrate(rhsColdEnergy, startDE, xFrom, xTo);
	currentDEC = startDE[4];
	(*DEC)[v].setValue(currentDEC*pow(Z, 5.0/3.0));

	calculatedDEC = true;
}
void FTTFQEmodel::calculateDS(Int v, Int t) {
	Double currentDS;
    DoubleVec startDS(0.0, rhsFTTFQEentropy::dim);
    Double a, b;
    Double V1, T1;
    Double xFrom = 1.0;
    Double xTo = 0.0;
	Int Tsize = T->size();

	V1 = (*V)[v]()*Z;
	T1 = (*T)[t]()*pow(Z, -4.0/3.0);

	startDS[0] = phi.valueAt_1();
	startDS[1] = startDS[0];
	startDS[2] = psi.valueAt_1();
	startDS[3] = startDS[2];
	startDS[4] = 1.0/(sqrt(2.0*T1)*3.0*M_PI)*psi.derivativeAt_0();
    			 
    a =   pow(2.0, 7.0/6.0)
        * pow(3.0, 2.0/3.0)
        * pow(M_PI, -5.0/3.0)
        * sqrt(T1)*pow(V1, 2.0/3.0);
    b = T1*V1/M_PI/M_PI/M_PI;
	rhsEntropy.updateParameters(a, b);

	entropySolver.Integrate(rhsEntropy, startDS, xFrom, xTo);
	currentDS = startDS[4];
	(*DS)[v*Tsize + t].setValue(currentDS*pow(Z, 1.0/3.0));
    
	calculatedDS  = true;
}
void FTTFQEmodel::calculateDST(Int v, Int t) {
	if (!calculatedDS)  calculateDS(v, t);
	if (!calculatedDSC) calculateDSC(v);
	Int DSTsize = DS->size();
	Int Tsize = T->size();
	Double currentDST;
	currentDST = (*DS)[v*Tsize + t]() - (*DSC)[v]();
	(*DST)[v*Tsize + t].setValue(currentDST);
	calculatedDST = true;
}
void FTTFQEmodel::calculateDSC(Int v) {
	Double currentDSC;
	DoubleVec startDS(0.0, rhsFTTFQEentropy::dim);
	Double a, b;
	Double V1, T1;
	Double xFrom = 1.0;
	Double xTo = 0.0;
	
	T1 = coldT()*pow(Z, -4.0/3.0);
	V1 = (*V)[v]()*Z;

	startDS[0] = coldPhi.valueAt_1();
	startDS[1] = startDS[0];
	startDS[2] = coldPsi.valueAt_1();
	startDS[3] = startDS[2];
	startDS[4] = 1.0/(sqrt(2.0*T1)*3.0*M_PI)*coldPsi.derivativeAt_0();
    			 
    a =   pow(2.0, 7.0/6.0)
        * pow(3.0, 2.0/3.0)
        * pow(M_PI, -5.0/3.0)
        * sqrt(T1)*pow(V1, 2.0/3.0);
    b = T1*V1/M_PI/M_PI/M_PI;
	rhsColdEntropy.updateParameters(a, b);

	coldEntropySolver.Integrate(rhsEntropy, startDS, xFrom, xTo);
	currentDSC = startDS[4];
	(*DSC)[v].setValue(currentDSC*pow(Z, 1.0/3.0));

	calculatedDSC = true;
}
void FTTFQEmodel::calculateDM(Int v, Int t) {
	Int Tsize = T->size();
	Double currentDM;
	currentDM = sqrt((*T)[t]()*0.5)/3.0/M_PI
        	  * (0.5*FDmhalf(phi.valueAt_1()) 
         	     + psi.valueAt_1());
	(*DM)[v*Tsize + t].setValue(currentDM);
	calculatedDM  = true;
}
void FTTFQEmodel::calculateDMT(Int v, Int t) {
	if (!calculatedDM)  calculateDM(v, t);
	if (!calculatedDMC) calculateDMC(v);
	Int Tsize = T->size();
	Double currentDMT;
	currentDMT = (*DM)[v*Tsize + t]() - (*DMC)[v]();
	(*DMT)[v*Tsize + t].setValue(currentDMT);
	calculatedDMT = true;
}
void FTTFQEmodel::calculateDMC(Int v) {
	Double currentDMC;
	currentDMC = sqrt(coldT()*0.5)/3.0/M_PI
				* (0.5*FDmhalf(coldPhi.valueAt_1()) 
					+ coldPsi.valueAt_1());
	(*DMC)[v].setValue(currentDMC);
	calculatedDMC = true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
// FTTFQE
/////////////////////////////////////////////////////////////////////////////////////////////////////
void FTTFQEmodel::calculateP(Int v, Int t) {
	if (!calculatedDP) calculateDP(v, t);
	Int Tsize = T->size();
	Double currentP;
	currentP = (*P0)[v*Tsize + t]() + (*DP)[v*Tsize + t]();
	(*P)[v*Tsize + t].setValue(currentP);
	calculatedP  = true;
}
void FTTFQEmodel::calculatePT(Int v, Int t) {
	if (!calculatedDP) calculateDP(v, t);
	if (!calculatedDPC) calculateDPC(v);
	Int Tsize = T->size();
	Double currentPT;
	currentPT = (*P0)[v*Tsize + t]() - (*PC0)[v]() + (*DP)[v*Tsize + t]() - (*DPC)[v]();
	(*PT)[v*Tsize + t].setValue(currentPT);
	calculatedPT = true;
}
void FTTFQEmodel::calculatePC(Int v) {
	if (!calculatedDPC) calculateDPC(v);
	Double currentPC;
	currentPC = (*PC0)[v]() + (*DPC)[v]();
	(*PC)[v].setValue(currentPC);
	calculatedPC = true;
}
void FTTFQEmodel::calculateE(Int v, Int t) {
	if (!calculatedDE) calculateDE(v, t);
	Int Tsize = T->size();
	Double currentE;
	currentE = (*E0)[v*Tsize + t]() + (*DE)[v*Tsize + t]();
	(*E)[v*Tsize + t].setValue(currentE);
	calculatedE = true;
}
void FTTFQEmodel::calculateET(Int v, Int t) {
	if (!calculatedDE)  calculateDE(v, t);
	if (!calculatedDEC) calculateDEC(v);
	Int Tsize = T->size();
	Double currentET;
	currentET = (*E0)[v*Tsize + t]() - (*EC0)[v]() + (*DE)[v*Tsize + t]() - (*DEC)[v]();
	(*ET)[v*Tsize + t].setValue(currentET);
	calculatedET = true;
}
void FTTFQEmodel::calculateEC(Int v) {
	if (!calculatedDEC) calculateDEC(v);
	Double currentEC;
	currentEC = (*EC0)[v]() + (*DEC)[v]();
	(*EC)[v].setValue(currentEC);
	calculatedEC = true;
}
void FTTFQEmodel::calculateS(Int v, Int t) {
	if (!calculatedDS) calculateDS(v, t);
	Int Tsize = T->size();
	Double currentS;
	currentS = (*S0)[v*Tsize + t]() + (*DS)[v*Tsize + t]();
	(*S)[v*Tsize + t].setValue(currentS);
	calculatedS = true;
}
void FTTFQEmodel::calculateST(Int v, Int t) {
	if (!calculatedDS)  calculateDS(v, t);
	if (!calculatedDSC) calculateDSC(v);
	Int Tsize = T->size();
	Double currentST;
	currentST = (*S0)[v*Tsize + t]() - (*SC0)[v]() + (*DS)[v*Tsize + t]() - (*DSC)[v]();
	(*ST)[v*Tsize + t].setValue(currentST);
	calculatedST = true;
}
void FTTFQEmodel::calculateSC(Int v) {
	if (!calculatedDSC) calculateDSC(v);
	Double currentSC;
	currentSC = (*SC0)[v]() + (*DSC)[v]();
	(*SC)[v].setValue(currentSC);
	calculatedSC = true;
}
void FTTFQEmodel::calculateM(Int v, Int t) {
	if (!calculatedDM) calculateDM(v, t);
	Int Tsize = T->size();
	Double currentM;
	currentM = (*M0)[v*Tsize + t]() + (*DM)[v*Tsize + t]();
	(*M)[v*Tsize + t].setValue(currentM);
	calculatedM  = true;
}
void FTTFQEmodel::calculateMT(Int v, Int t) {
	if (!calculatedDM)  calculateDM(v, t);
	if (!calculatedDMC) calculateDMC(v);
	Int Tsize = T->size();
	Double currentMT;
	currentMT = (*M0)[v*Tsize + t]() - (*MC0)[v]() + (*DM)[v*Tsize + t]() - (*DMC)[v]();
	(*MT)[v*Tsize + t].setValue(currentMT);
	calculatedMT = true;
}
void FTTFQEmodel::calculateMC(Int v) {
	if (!calculatedDMC) calculateDMC(v);
	Double currentMC;
	currentMC = (*MC0)[v]() + (*DMC)[v]();
	(*MC)[v].setValue(currentMC);
	calculatedMC = true;
}
void FTTFQEmodel::calculateQEoS(Int v, Int t) {
	if (!calculatedE) calculateE(v, t);
	if (!calculatedP) calculateP(v, t);
	Int Tsize = T->size();
	Double currentQ;
	currentQ = (*E)[v*Tsize + t]() - (*P)[v*Tsize + t]()*(*V)[v]();
	(*QEoS)[v*Tsize + t].setValue(currentQ);
	calculatedQEoS = true;
}
void FTTFQEmodel::calculateTEoS(Int v, Int t) {
	if (!calculatedP) calculateP(v, t);
	Int Tsize = T->size();
	Double currentTEoS;
	currentTEoS = (*P)[v*Tsize + t]()*(*V)[v]()/(*T)[t]();
	(*TEoS)[v*Tsize + t] = currentTEoS;
	calculatedTEoS = true;
}