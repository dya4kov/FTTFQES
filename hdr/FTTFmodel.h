#pragma once
#include "../lib/numeric/hdr/ODE/ODEsolver.h"
#include "../lib/numeric/hdr/ODE/ODEdata.h"
#include "../lib/numeric/hdr/ODE/steppers/ODEstepperPD853.h"
#include "../lib/numeric/hdr/numericTypes.h"
#include "../lib/numeric/hdr/mathUtils.h"
#include "../lib/numeric/hdr/SF/FermiDirac.h"
#include "../hdr/ThermodynamicFunction.h"
#include "../hdr/FTTFpotential.h"
/**
* @brief This class implements interface for Thomas-Fermi model.
* @details It implements virtual methods for calculating thermodynamic values.
*/
class FTTFmodel {
public:
	/**
	* @brief A constructor of Thomas-Fermi model.
	* @details It allocates memory for ODE solvers of Thomas-Fermi energy and entropy
	* and for Thomas-Fermi potential.
	*/
    FTTFmodel(Double Z = 1.0);

    setParameters()

    /**
	* @brief Problem for calculating energy.
	* @details The calculation goes through solving the ODE with parameters 
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
		static const Int dim = 2;
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
	* @details The calculation goes through solving the ODE with parameters
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
		static const Int dim = 2;
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

	

protected:
	/**
	* @brief Implemented method for calculating pressure in the Thomas-Fermi model.
	* @details The formula for calculation: 
	* @f[
	*	P = \frac{(2T)^{5/2}}{6\pi^2}I_{3/2}(\phi(1)).
	* @f]
	*/
    void calculatePressure();
	/**
	* @brief Implemented method for calculating energy in the Thomas-Fermi model.
	* @details The calculation goes through solving the ODE with parameters 
	* @f$ p_1 = 2^{7/6}3^{2/3}\pi^{-5/3}T^{1/2}V^{2/3} @f$,
	* @f$ p_2 = 3\sqrt{2}\pi^{-2}VT^{5/2} @f$,
	* @f$ E_0 = -0.76874512422 @f$. Thomas-Fermi energy is solution of the ODE:
	* @f{eqnarray*}{
	*  \frac{d^2\phi}{dx^2} &=& p_1 x I_{1/2}
	*	  \left(
	*		\frac{\phi(x)}{x}
	*	  \right), \\
	*  E'(x) &=& p_2 x^2 I_{3/2}
	*				\left(
	*				  \frac{\phi(x)}{x}
	*				\right), \\
	* \phi(1) &=& \phi'(1), \\
	* E(1) &=& \frac{2\sqrt{2}VT^{5/2}}{\pi^2}I_{3/2}(\phi(1)) - E_0,
	* @f}
	* Finally, @f$ E = E(0) @f$.
	*/
    void calculateEnergy();
	/**
	* @brief Implemented method for calculating entropy in the Thomas-Fermi model.
	* @details The calculation goes through solving the ODE with parameters
	* @f$ p_1 = 2^{7/6}3^{2/3}\pi^{-5/3}T^{1/2}V^{2/3} @f$,
	* @f$ p_2 = 7\sqrt{2}\pi^{-2}VT^{3/2} @f$,
	* Thomas-Fermi entropy is solution of the ODE:
	* @f{eqnarray*}{
	*  \frac{d^2\phi}{dx^2} &=& p_1 x I_{1/2}
	*	  \left(
	*		\frac{\phi(x)}{x}
	*	  \right), \\
	*  S'(x) &=& p_2 x^2 I_{3/2}
	*				\left(
	*				  \frac{\phi(x)}{x}
	*				\right), \\
	* \phi(1) &=& \phi'(1), \\
	* S(1) &=& \frac{4\sqrt{2}VT^{3/2}}{\pi^2}I_{3/2}(\phi(1)) - \phi'(0),
	* @f}
	* Finally, @f$ S = S(0) @f$.
	*/
    void calculateEntropy();
	/**
	* @brief Implemented method for calculating chemical potential 
	* in the Thomas-Fermi model.
	* @details The formula for calculation:
	* @f[
	*	\mu = T\phi(1).
	* @f]
	*/
    void calculateChemPotential();
	/**
	* @brief Implemented method for calculating thermal pressure in the Thomas-Fermi model.
	* @details The formula for calculation:
	* @f[
	*	P_T = P - \left.P\right|_{T = 0}.
	* @f]
	*/
    void calculateThermalP();
	/**
	* @brief Implemented method for calculating thermal energy in the Thomas-Fermi model.
	* @details The formula for calculation:
	* @f[
	*	E_T = E - \left.E\right|_{T = 0}.
	* @f]
	*/
    void calculateThermalE();
	/**
	* @brief Implemented method for calculating thermal entropy in the Thomas-Fermi model.
	* @details The formula for calculation:
	* @f[
	*	S_T = S - \left.S\right|_{T = 0}.
	* @f]
	*/
    void calculateThermalS();
	/**
	* @brief Implemented method for calculating thermal chemical potential in the Thomas-Fermi model.
	* @details The formula for calculation:
	* @f[
	*	\mu_T = \mu - \left.\mu\right|_{T = 0}.
	* @f]
	*/
    void calculateThermalCP();
	/**
	* @brief Implemented method which sets all necessary parameters for Thomas-Fermi potential.
	*/
    void getReady();

private:

    FTTFpotential potential;
    FTTFpotential coldPotential;
    ODESolver<TF::Energy>* energySolver;
    ODESolver<TF::Entropy>* entropySolver;

    Temperature& makeChargeless(Temperature& T);
    Volume& makeChargeless(Volume& V);
    Pressure& makeChargeless(Pressure& P);
    Energy& makeChargeless(Energy& E);
    Entropy& makeChargeless(Entropy& E);
    ChemicalPotential& makeChargeless(ChemicalPotential& M);

    Temperature& makeChargeful(Temperature& T);
    Volume& makeChargeful(Volume& V);
    Pressure& makeChargeful(Pressure& P);
    Energy& makeChargeful(Energy& E);
    Entropy& makeChargeful(Entropy& E);
    ChemicalPotential& makeChargeful(ChemicalPotential& M);
    
    double coldP();
    double coldE();
    double coldS();
    double coldM();
};