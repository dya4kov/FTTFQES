#pragma once
#include <string>

enum Unit {
    dimensionless,
    Atomic,
    eV,
    Kelvin,
    cmc,
    mc,
    Pa,
    GPa,
    MBar,
    gOverCmc,
    kgOverMc,
    oneOverCmc,
    oneOverMc,
    undefinedUnit,
    totalUnits
};

enum Scaling {
    lin = 0,
    lg = 1
};

std::string unitToString(Unit unit) {
    std::string s;
    switch (unit) {
        case dimensionless : 
            return s = "-";
        case Atomic :
            return s = "atm";
        case eV :
            return s = "eV";
        case Kelvin :
            return s = "K";
        case cmc :
            return s = "cm^{3}";
        case mc :
            return s = "m^{3}";
        case Pa :
            return s = "Pa";
        case GPa :
            return s = "GPa";
        case MBar :
            return s = "MBar";
        case gOverCmc :
            return s = "g/cm^3";
        case kgOverMc :
            return s = "kg/m^{3}";
        case oneOverCmc :
            return s = "cm^{-3}";
        case oneOverMc :
            return s = "m^{-3}";
        case undefinedUnit :
            return s = "-";
    }
}

Unit stringToUnit(std::string unit) {
    if (!unit.compare("-")) return dimensionless;
    if (!unit.compare("atm")) return Atomic;
    if (!unit.compare("eV")) return eV;
    if (!unit.compare("K")) return Kelvin;
    if (!unit.compare("cm^{3}")) return cmc;
    if (!unit.compare("m^{3}")) return mc;
    if (!unit.compare("Pa")) return Pa;
    if (!unit.compare("GPa")) return GPa;
    if (!unit.compare("MBar")) return MBar;
    if (!unit.compare("g/cm^3")) return gOverCmc;
    if (!unit.compare("kg/m^{3}")) return kgOverMc;
    if (!unit.compare("cm^{-3}")) return oneOverCmc;
    if (!unit.compare("m^{-3}")) return oneOverMc;
    return undefinedUnit;
}

Scaling stringToScale(std::string scale) {
    if (!scale.compare("lin")) return lin;
    if (!scale.compare("log")) return lg;
}

std::string scaleToString(Scaling scale) {
    std::string s;
    switch (scale) {
        case lin : 
            return s = "lin";
        case lg :
            return s = "log";
    }
}

const double Avogadro = 6.02214129e+23;          // [mol^{-1}]
const double eMass = 9.10938291e-28;             // [g]
const double rBohr = 5.2917720859e-11;           // [m]
const double aVol = rBohr*rBohr*rBohr;           // [m^3]
const double Ehartree = 27.2113845;              // [eV]
const double kBoltzmann = 8.6173324e-5/Ehartree; // [hartree*K^{-1}]
const double qElectron = 1.60217657e-19;         // [C] coulomb
const double aPress = Ehartree*qElectron/aVol;   // [Pa]
const double aDens = eMass/(aVol*1e+6);          // [g*cm^{-3}]