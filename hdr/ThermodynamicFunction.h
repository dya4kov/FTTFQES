#pragma once
#include "../lib/numeric/hdr/numericTypes.h"
#include "../lib/numeric/hdr/mathUtils.h"
#include "../hdr/Units.h"
#include <vector>
/**
 * @brief This class provides interface for thermodynamic functions.
 * @details Contains value, unit and scale (linear or logarithmic) of the thermodynamic function. 
 */
class PhysicalQuantity {
public:
    /**
    * @brief Constructor of the thermodynamic function.
    */
    PhysicalQuantity(Double value = 0.0, Unit unit = Atomic, Scaling scale = lin) {
        setValue(value, unit, scale);
    };
    /**
    * @brief Set the value of the thermodynamic function.
    */
    void setValue(Double value, Unit unit = Atomic, Scaling scale = lin) {
        Double newValue;
        if (scale == lg) newValue = pow(10, value);
        else newValue = value;
        switch (unit) {
            case Atomic : 
                defaultValue = newValue;
                break;
            case cmc : 
                defaultValue = newValue/(aVol*1e+6);
                break;
            case mc : 
                defaultValue = newValue/aVol;
                break;
            case gOverCmc : 
                defaultValue = newValue/aDens;
                break;
            case kgOverMc : 
                defaultValue = newValue/(1000*aDens);
                break;
            case Kelvin : 
                defaultValue = newValue*kBoltzmann;
                break;
            case oneOverCmc : 
                defaultValue = newValue*aVol*1e+6;
                break;
            case oneOverMc : 
                defaultValue = newValue*aVol;
                break;
            case MBar : 
                defaultValue = newValue/(aPress*1e-11);
                break;
            case GPa : 
                defaultValue = newValue/(aPress*1e-9);
                break;
            case Pa:
                defaultValue = newValue/aPress;
                break;
            case eV : 
                defaultValue = newValue/Ehartree;
                break;
        }
    }
    /**
	* @brief Select units and get the value of the thermodynamic function.
	*/
    Double operator() (Unit unit = Atomic, Scaling scale = lin) const {
        Double returnValue;
        switch (unit) {
            case Atomic : 
                returnValue = defaultValue;
                break;
            case cmc : 
                returnValue = defaultValue*aVol*1e+6;
                break;
            case mc :
                returnValue = defaultValue*aVol;
                break;
            case gOverCmc :
                returnValue = defaultValue*aDens;
                break;
            case kgOverMc :
                returnValue = defaultValue*aDens*1000;
                break;
            case Kelvin : 
                returnValue = defaultValue/kBoltzmann;
                break;
            case oneOverCmc : 
                returnValue = defaultValue/(aVol*1e+6);
                break;
            case oneOverMc : 
                returnValue = defaultValue/aVol;
                break;
            case GPa : 
                returnValue = defaultValue*aPress*1e-9;
                break;
            case MBar :
                returnValue = defaultValue*aPress*1e-11;
                break;
            case Pa:
                returnValue = defaultValue*aPress;
                    break;
            case eV: 
                returnValue = defaultValue*Ehartree;
                break;
        }
        if (scale == lg) returnValue = log10(abs(returnValue));
        return returnValue;
    }
private:
    Double defaultValue;
};

typedef PhysicalQuantity PhysQ;
typedef std::vector<PhysQ> PhysQvec;
