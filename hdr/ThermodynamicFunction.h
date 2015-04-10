#pragma once
#include "../lib/numeric/hdr/numericTypes.h"
#include "../lib/numeric/hdr/mathUtils.h"
#include "../hdr/Units.h"
#include <vector>
/**
 * @brief This class provides interface for thermodynamic functions.
 * @details Contains value, unit and scale (linear or logarithmic) of the thermodynamic function. 
 */
template<class F>
class ThermodynamicFunction {
public:
    /**
    * @brief Constructor of the thermodynamic function.
    */
    ThermodynamicFunction(Double value = 0.0, Unit unit = Atomic, Scaling scale = lin) : instance(value, unit, scale) {};
    /**
    * @brief Set the value of the thermodynamic function.
    */
    void setValue(Double value, Unit unit = Atomic, Scaling scale = lin) {
        instance.SetValue(value, unit, scale);
    }
    /**
	* @brief Select units and get the value of the thermodynamic function.
	*/
    Double operator() (Unit unit = Atomic, Scaling scale = lin) const {
        return instance.GetValue(unit, scale);
    }
private:
    F instance;
};

struct volume {
    volume(Double value, Unit unit = Atomic, Scaling scale = lin) {
        SetValue(value, unit, scale);
    }
    void SetValue(Double value, Unit unit, Scaling scale) {
        Double newValue;
        if (scale == lg) newValue = pow(10, value);
        else newValue = value;
        switch (unit) {
            case Atomic : 
                defaultValue = newValue;
                break;
            case cmc : 
                defaultValue = newValue/1.4818e-25;
                break;
            case mc : 
                defaultValue = newValue/1.4818e-31;
                break;
        }
    }
    Double GetValue(Unit unit, Scaling scale) const {
        Double returnValue;
        switch (unit) {
            case Atomic : 
                returnValue = defaultValue;
                break;
            case cmc : 
                returnValue = defaultValue*1.4818e-25;
                break;
            case mc :
                returnValue = defaultValue*1.4818e-31;
                break;
        }
        if (scale == lg) returnValue = log10(abs(returnValue));
        return returnValue;
    }
private:
    Double defaultValue;
};

struct density {
    density(Double value, Unit unit = Atomic, Scaling scale = lin) { 
        SetValue(value, unit, scale);
    }
    void SetValue(Double value, Unit unit, Scaling scale) {
        Double newValue;
        if (scale == lg) newValue = pow(10, value);
        else newValue = value;

        switch (unit) {
            case gOverCmc : 
                defaultValue = newValue;
                break;
            case kgOverMc : 
                defaultValue = newValue*1000;
                break;
        }
    }
    Double GetValue(Unit unit, Scaling scale) const {
        Double returnValue;
        switch (unit) {
            case gOverCmc :
                returnValue = defaultValue;
                break;
            case kgOverMc :
                returnValue = defaultValue*1000;
                break;
        }
        if (scale == lg) returnValue = log10(abs(returnValue));
        return returnValue;   
    }
private:
    Double defaultValue;
};

struct concentration {
    concentration(Double value, Unit unit = Atomic, Scaling scale = lin) { 
        SetValue(value, unit, scale);
    }
    void SetValue(Double value, Unit unit, Scaling scale) {
        Double newValue;
        if (scale == lg) newValue = pow(10, value);
        else newValue = value;

        switch (unit) {
            case Atomic : 
                defaultValue = newValue;
                break;
            case oneOverCmc : 
                defaultValue = newValue*1.4818e-25;
                break;
            case oneOverMc : 
                defaultValue = newValue*1.4818e-31;
                break;
        }
    }
    Double GetValue(Unit unit, Scaling scale) const {
        Double returnValue;
        switch (unit) {
            case Atomic : 
                returnValue = defaultValue;
                break;
            case oneOverCmc : 
                returnValue = defaultValue/1.4818e-25;
                break;
            case oneOverMc : 
                returnValue = defaultValue/1.4818e-31;
                break;
        }
        if (scale == lg) returnValue = log10(abs(returnValue));
        return returnValue;
    }
private:
    Double defaultValue;    
};

struct temperature {
    temperature(Double value, Unit unit = Atomic, Scaling scale = lin) { 
        SetValue(value, unit, scale);
    }
    void SetValue(Double value, Unit unit, Scaling scale) {
        Double newValue;
        if (scale == lg) newValue = pow(10, value);
        else newValue = value;
        switch (unit) {
            case Atomic : 
                defaultValue = newValue;
                break;
            case eV : 
                defaultValue = newValue/27.229;
                break;
            case Kelvin : 
                defaultValue = newValue*8.6173324e-5/27.229;
                break;
        }
    }
    Double GetValue(Unit unit, Scaling scale) const {
        Double returnValue;
        switch (unit) {
            case Atomic : 
                returnValue = defaultValue;
                break;
            case eV : 
                returnValue = defaultValue*27.229;
                break;
            case Kelvin : 
                returnValue = defaultValue*27.229/8.6173324e-5;
                break;
        }
        if (scale == lg) returnValue = log10(abs(returnValue));
        return returnValue;        
    }
private:
    Double defaultValue;
};

struct pressure {
    pressure(Double value, Unit unit = Atomic, Scaling scale = lin) { 
        SetValue(value, unit, scale);
    }
    void SetValue(Double value, Unit unit, Scaling scale) {
        Double newValue;
        if (scale == lg) newValue = pow(10, value);
        else newValue = value;

        switch (unit) {
            case Atomic : 
                defaultValue = newValue;
                break;
            case MBar : 
                defaultValue = newValue/294.18;
                break;
            case GPa : 
                defaultValue = newValue/29418.0;
                break;
            case Pa:
                defaultValue = newValue/29418e+9;
                break;
        }
    }
    Double GetValue(Unit unit, Scaling scale) const {
        Double returnValue;
        switch (unit) {
            case Atomic : 
                returnValue = defaultValue;
                break;
            case GPa : 
                returnValue = defaultValue*29418;
                break;
            case MBar :
                returnValue = defaultValue*294.18;
                break;
            case Pa:
                returnValue = defaultValue*29418e+9;
                    break;
        }
        if (scale == lg) returnValue = log10(abs(returnValue));
        return returnValue;
    }
private:
    Double defaultValue;    
};

struct energy {
    energy(Double value, Unit unit = Atomic, Scaling scale = lin) { 
        SetValue(value, unit, scale);
    }
    void SetValue(Double value, Unit unit, Scaling scale) {
        Double newValue;
        if (scale == lg) newValue = pow(10, value);
        else newValue = value;
        switch (unit) {
            case Atomic : 
                defaultValue = newValue;
                break;
            case eV : 
                defaultValue = newValue/27.229;
                break;
        }
    }
    Double GetValue(Unit unit, Scaling scale) const {
        Double returnValue;
        switch (unit) {
            case Atomic : 
                returnValue = defaultValue;
                break;
            case eV : 
                returnValue = defaultValue*27.229;
                break;
        }
        if (scale == lg) returnValue = log10(abs(returnValue));
        return returnValue;   
    }
private:
    Double defaultValue;
};

struct entropy {
    entropy(Double value, Unit unit = Atomic, Scaling scale = lin) { 
        SetValue(value, unit, scale);
    }
    void SetValue(Double value, Unit unit, Scaling scale) {
        Double newValue;
        if (scale == lg) newValue = pow(10, newValue);
        else newValue = value;
        switch (unit) {
            case Atomic : 
                defaultValue = newValue;
                break;
        }
    }
    Double GetValue(Unit unit, Scaling scale) const {
        Double returnValue;
        switch (unit) {
            case Atomic : 
                returnValue = defaultValue;
                break;
            }
        if (scale == lg) returnValue = log10(abs(returnValue));
        return returnValue;
    }
private:
    Double defaultValue;
};

struct chemicalPotential {
    chemicalPotential(Double value, Unit unit = Atomic, Scaling scale = lin) { 
        SetValue(value, unit, scale);
    }
    void SetValue(Double value, Unit unit, Scaling scale) {
        Double newValue;
        if (scale == lg) newValue = pow(10, value);
        else newValue = value;

        switch (unit) {
            case Atomic : 
                defaultValue = newValue;
                break;
            case eV : 
                defaultValue = newValue/27.229;
                break;
        }
    }
    Double GetValue(Unit unit, Scaling scale) const {
        Double returnValue;
        switch (unit) {
            case Atomic : 
                returnValue = defaultValue;
                break;
            case eV : 
                returnValue = defaultValue*27.229;
                break;
        }
        if (scale == lg) returnValue = log10(abs(returnValue));
        return returnValue;
    }
private:
    Double defaultValue;
};

typedef ThermodynamicFunction<volume> Volume;
typedef ThermodynamicFunction<density> Density;
typedef ThermodynamicFunction<temperature> Temperature;
typedef ThermodynamicFunction<concentration> Concentration;
typedef ThermodynamicFunction<pressure> Pressure;
typedef ThermodynamicFunction<energy> Energy;
typedef ThermodynamicFunction<entropy> Entropy;
typedef ThermodynamicFunction<chemicalPotential> ChemicalPotential;

typedef std::vector<Volume> VolVec;
typedef std::vector<Temperature> TempVec;
typedef std::vector<Density> DensVec;
typedef std::vector<Concentration> ConVec;
typedef std::vector<Pressure> PressVec;
typedef std::vector<Energy> EnergyVec;
typedef std::vector<Entropy> EntropyVec;
typedef std::vector<ChemicalPotential> ChemPotVec;