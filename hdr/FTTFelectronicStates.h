#pragma once
#include "../lib/numeric/hdr/ODE/ODEsolver.h"
#include "../lib/numeric/hdr/ODE/ODEdata.h"
#include "../lib/numeric/hdr/ODE/steppers/ODEstepperPD853.h"
#include "../lib/numeric/hdr/numericTypes.h"
#include "../lib/numeric/hdr/mathUtils.h"
#include "../lib/numeric/hdr/SF/FermiDirac.h"
#include "../hdr/FTTFpotential.h"
#include "../hdr/ThermodynamicFunction.h"
#include "../hdr/FDhalfInc.h"
#include <fstream>
#include <iostream>

class FTTFelectronicStates {
public:
	FTTFelectronicStates(const Int _nMax);
	void setParameters(const Volume &V, const Temperature &T, const Double Z);
    // energy range for finding levels
    void setEnergyRange(const Double Yleft, const Double Yright);
    // get energy level n,l
    Double operator() (const Int n, const Int l);
    // number of states with user BE
    Double NTF(Double BE);
    Double N(Double BE);
    Double DN(Double BE);
    // number of states with optimal BE
    Double NTF();
    Double N();
	Double DN();

	struct RHSAction {
		static const Int dim = 3;
		FermiDirac<Half> FDhalf;
		RHSAction(Double _a = 0, Double eArg = 0, Double lArg = 0) : 
			a(_a), energyArg(eArg), lambdaArg(lArg) {}
		void setA(const Double _a) { a = _a; }
		void setLArg(const Double lArg) { lambdaArg = lArg; }
		void setEArg(const Double eArg) { energyArg = eArg; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
		    dydx[0] = y[1];
		    if (x > 0) {
		        dydx[1] = a*x*FDhalf(y[0] / x);
		        p = energyArg + y[0] / x - lambdaArg / x / x;
		        dydx[2] = p > 0 ? -sqrt(p) : 0;
		    }
		    else {
		        dydx[1] = 1e+10;
		        dydx[2] = 0;
		    }
		}
	private:
		Double a;
		Double energyArg;
		Double lambdaArg;
		Double p;
	};

	struct RHSRP2 {
		static const Int dim = 2;
		FermiDirac<Half> FDhalf;
		Double xUp;
		Double xDown;
		Double phi_xUp;
		Double dphi_xUp;
		Double dphi_xDown;
		Double phi_xDown;
		RHSRP2(Double _a = 0, Double eArg = 0, Double lArg = 0) : 
			a(_a), energyArg(eArg), lambdaArg(lArg) {
                yOld = -1;
                dyOld = -1;
                p1 = -1;
                p2 = -1;
                xOld = 1.0;
            }
		void setA(const Double _a) { a = _a; }
		void setLArg(const Double lArg) { lambdaArg = lArg; }
		void setEArg(const Double eArg) { energyArg = eArg; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {

		    dydx[0] = y[1];
		    if (x > 0) {
                dydx[1] = a*x*FDhalf(y[0] / x);  
                p2 = energyArg + y[0] / x - lambdaArg / x / x;
            } 
            else {
                dydx[1] = 1e+10;
            }
            
            if (p1 < 0 && p2 >= 0 && x < xOld) {
                phi_xUp = yOld;
                dphi_xUp = dyOld;
                phi_xDown = y[0];
                dphi_xDown = y[1];
                xUp = xOld;
                xDown = x;
            }

		    p1 = p2;
		    yOld = y[0];
		    dyOld = y[1];
		    xOld = x;
		}
	private:
		Double a;
		Double energyArg;
		Double lambdaArg;
		Double p1;
		Double p2;
		Double yOld;
		Double dyOld;
		Double xOld;
	};

	struct RHSRP1 {
		Double xUp;
		Double xDown;
		Double phi_xUp;
		Double dphi_xUp;
		Double dphi_xDown;
		Double phi_xDown;
		static const Int dim = 2;
		FermiDirac<Half> FDhalf;
		RHSRP1(Double _a = 0, Double eArg = 0, Double lArg = 0) : 
			a(_a), energyArg(eArg), lambdaArg(lArg) {
                yOld = -1;
                dyOld = -1;
                p1 = -1;
                p2 = -1;
            }
		void setA(const Double _a) { a = _a; }
		void setLArg(const Double lArg) { lambdaArg = lArg; }
		void setEArg(const Double eArg) { energyArg = eArg; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {

		    dydx[0] = y[1];
		    if (x > 0) dydx[1] = a*x*FDhalf(y[0] / x);
		    else {
		        dydx[1] = 1e+10;
		    }

		    p2 = energyArg + y[0] / x - lambdaArg / x / x;

		    if (p1 >= 0 && p2 < 0 && x < xOld) {
		        phi_xUp = yOld;
		        dphi_xUp = dyOld;
		        phi_xDown = y[0];
		        dphi_xDown = y[1];
		        xUp = xOld;
		        xDown = x;
		    }

		    p1 = p2;
		    yOld = y[0];
		    dyOld = y[1];
		    xOld = x;
		}
	private:
		Double a;
		Double energyArg;
		Double lambdaArg;
		Double p1;
		Double p2;
		Double yOld;
		Double dyOld;
		Double xOld;
	};

	struct RHSBoundEnergy {
		static const Int dim = 3;
		FermiDirac<Half> FDhalf;
		RHSBoundEnergy(Double _a = 0, Double _boundEnergy = 0) : 
			a(_a), boundEnergy(_boundEnergy) {}
		void setA(const Double _a) { a = _a; }
		void setBE(const Double BE) { boundEnergy = BE; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
		    dydx[0] = y[1];
		    if (x > 0) {
		        dydx[1] = a*x*FDhalf(y[0] / x);
		        if (boundEnergy + y[0] / x <= 0) dydx[2] = 0.0;
		        else {
					dydx[2] = -2.0 / 3.0*x*x*pow(boundEnergy + y[0] / x, 3.0 / 2.0);
		        }
		    }
		    else {
		        dydx[1] = 1e+10;
		        dydx[2] = 0;
    		}
		}
	private:
		Double a;
		Double boundEnergy;
	};

	struct RHSFTTFstates {
		static const Int dim = 3;
		FermiDirac<Half> FDhalf;
		FDhalfInc FDhalfI;
		RHSFTTFstates(Double _a = 0, Double _boundEnergy = 0) : 
			a(_a), boundEnergy(_boundEnergy), FDhalfI(1e-7) {}
		void setA(const Double _a) { a = _a; }
		void setBE(const Double BE) { boundEnergy = BE; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
		    dydx[0] = y[1];
		    if (x > 0) {
		        dydx[1] = a*x*FDhalf(y[0] / x);
		        if (boundEnergy + y[0] / x <= 0) dydx[2] = 0.0;
		        else    
		            dydx[2] = -x*x*FDhalfI(y[0] / x, 0, boundEnergy + y[0] / x);
		    }
		    else {
		        dydx[1] = 1e+10;
		        dydx[2] = 0;
		    }
		}
	private:
		Double a;
		Double boundEnergy;	
	};

private:
	void calculateAction(Double eArg, Double lArg); Double action;
	void calculateRP2(Double eArg, Double lArg);    Double rp2; Double phi_rp2; Double dphi_rp2;
	void calculateRP1(Double eArg, Double lArg);    Double rp1; 

	void calculateNTF(Double boundEnergy); Double NTFval; bool calculated_NTF;
	void calculateN(Double boundEnergy);   Double Nval;   bool calculated_N;
	void calculateBoundEnergy(); Double boundEnergy; bool calculated_BE;
    void calculateEnergyLevel(Int n, Int l); DoubleMat energyLevel; bool** calculatedLevel;
    void calculateEnergyLevels(); 

    Double boundEnergyRoot(Double e1, Double e2); DoubleVec BERoot;
	Double mixBound();
	Double pseudoStates(Double energy);  
	Double pseudoStatesTF(Double energy);

    // ode solvers
    RHSAction rhsAction;
    ODEsolver<ODEstepperPD853<RHSAction> > actionSolver;
    ODEdata actionData;

    RHSRP1 rhsRP1;
    ODEsolver<ODEstepperPD853<RHSRP1> > RP1Solver;
    ODEdata RP1Data;

    RHSRP2 rhsRP2;
    ODEsolver<ODEstepperPD853<RHSRP2> > RP2Solver;
    ODEdata RP2Data;

    RHSBoundEnergy rhsBoundEnergy;
    ODEsolver<ODEstepperPD853<RHSBoundEnergy> > boundEnergySolver;
    ODEdata boundEnergyData;

    RHSFTTFstates rhsFTTFstates;
    ODEsolver<ODEstepperPD853<RHSFTTFstates> > FTTFstatesSolver;
    ODEdata FTTFstatesData;

    FTTFpotential phi;
        
    Double startYleft;
    Double startYright;

    Double r0;

    Volume V;
    Temperature T;

    Int nMax;
    Double eps; 
};

FTTFelectronicStates::FTTFelectronicStates(const Int _nMax) : 
                nMax(_nMax), 
                actionData(-1),
                actionSolver(1e-6),
                rhsAction(0.0),
                RP1Data(-1),
                RP1Solver(1e-6),
                rhsRP1(0.0),
                RP2Data(-1),
                RP2Solver(1e-6),
                rhsRP2(0.0),
                boundEnergyData(-1),
                boundEnergySolver(1e-6),
                rhsBoundEnergy(0.0),
                FTTFstatesData(-1),
                FTTFstatesSolver(1e-6),
                rhsFTTFstates(0.0),
                energyLevel(_nMax + 1, _nMax),
                BERoot(0.0, _nMax),
                eps(1e-7),
                startYleft(-1000.0),
                startYright(100.0) {
                    actionSolver.SetOutput(actionData);
                    RP1Solver.SetOutput(RP1Data);
                    RP2Solver.SetOutput(RP2Data);
                    boundEnergySolver.SetOutput(boundEnergyData);
                    FTTFstatesSolver.SetOutput(FTTFstatesData);
                    calculatedLevel = new bool*[_nMax + 1];
                    for (Int i = 0; i <= _nMax; ++i) {
                        calculatedLevel[i] = new bool[_nMax];
                        for (Int j = 0; j < _nMax; ++j) {
                            calculatedLevel[i][j] = false;
                        }
                    }
                }

void FTTFelectronicStates::setParameters(const Volume &_V, const Temperature &_T, const Double Z) {
    
    if (abs(log10(_V()) - log10(V())) > 1e-10
        || abs(log10(_T()) - log10(T())) > 1e-10) 
    {
    	V = _V;
    	T = _T; 
        phi.setParameters(V, T, Z);
        phi.setTolerance(1e-10);
        Double a =   pow(2.0, 7.0/6.0)
                   * pow(3.0, 2.0/3.0)
                   * pow(M_PI, -5.0/3.0)
                   * sqrt(T())*pow(V(), 2.0/3.0);
        rhsAction.setA(a);
        rhsBoundEnergy.setA(a);
        rhsRP2.setA(a);
        rhsRP1.setA(a);
        rhsFTTFstates.setA(a);
        r0 = pow(3.0*V() / 4.0 / M_PI, 1.0 / 3.0);
        // reset levels
        for (Int i = 0; i <= nMax; ++i) {
            calculatedLevel[i] = new bool[nMax];
            for (Int j = 0; j < nMax; ++j) {
                calculatedLevel[i][j] = false;
            }
        }
        calculated_BE = false;
        calculated_NTF = false;
        calculated_N = false;
    }
}

void FTTFelectronicStates::setEnergyRange(const Double Yleft, const Double Yright) {
    startYleft = Yleft;
    startYright = Yright;
}

Double FTTFelectronicStates::operator() (const Int n, const Int l) {
    if (l >= n) {
        std::cerr << "l = " << l << " is greater than n = " << n << 
        ". Program stops. " << std::endl;
        exit(0);
    }
    if (!calculatedLevel[n][l]) calculateEnergyLevel(n, l);
    return energyLevel[n][l];
}

void FTTFelectronicStates::calculateAction(Double eArg, Double lArg) {
    DoubleVec y0(0.0, RHSAction::dim);
    
    rhsAction.setEArg(eArg);
    rhsAction.setLArg(lArg);
    
    y0[0] = phi_rp2;
    y0[1] = dphi_rp2;
    y0[2] = 0.0;
    
    actionSolver.SetTolerance(1e-10, 0.0);
    if (rp2 > rp1) {
        actionSolver.Integrate(rhsAction, y0, rp2, rp1);
    }
    
    action = r0*sqrt(2 * T())*y0[2];
}

void FTTFelectronicStates::calculateRP2(Double eArg, Double lArg) {
    DoubleVec y0(0.0, RHSRP2::dim);
    Double h = 0.0;

    rhsRP2.setEArg(eArg);
    rhsRP2.setLArg(lArg);

    y0[0] = y0[1] = phi.valueAt_1();
    rhsRP2.xUp = 0;
    rhsRP2.xDown = 0;
    
    if (eArg + y0[0] - lArg > 0) {
        phi_rp2 = y0[0];
        dphi_rp2 = y0[1];
        rp2 = 1.0;
    }
    else {
        RP2Solver.SetTolerance(0.0, 1e-10);
        RP2Solver.Integrate(rhsRP2, y0, 1.0, 0.0);
        while (abs(rhsRP2.xUp - rhsRP2.xDown) > 1e-12) {

            y0[0] = rhsRP2.phi_xUp;
            y0[1] = rhsRP2.dphi_xUp;
            h =  (rhsRP2.xUp - rhsRP2.xDown) / 11.0;
            RP2Solver.SetStep(h);
            RP2Solver.SetTolerance(0.0, h / 2.0);
            RP2Solver.Integrate(rhsRP2, y0, rhsRP2.xUp, rhsRP2.xDown);
        }
        phi_rp2 = rhsRP2.phi_xDown;
        dphi_rp2 = rhsRP2.dphi_xDown;
        rp2 = rhsRP2.xDown;
    }
}

void FTTFelectronicStates::calculateRP1(Double eArg, Double lArg) {
    DoubleVec y0(0.0, RHSRP1::dim);
    Double h = 0.0;

    rhsRP1.setEArg(eArg);
    rhsRP1.setLArg(lArg);

    y0[0] = phi_rp2;
    y0[1] = dphi_rp2;

    rhsRP1.xUp = 0;
    rhsRP1.xDown = 0;
    
    RP1Solver.SetTolerance(0.0, 1e-10);
    if (rp2 > 1e-12) {
        RP1Solver.Integrate(rhsRP1, y0, rp2, 0.0);

        while (abs(rhsRP1.xUp - rhsRP1.xDown) > 1e-12) {

            y0[0] = rhsRP1.phi_xUp;
            y0[1] = rhsRP1.dphi_xUp;
            h = (rhsRP1.xUp - rhsRP1.xDown) / 11.0;
            RP1Solver.SetStep(h);
            RP1Solver.SetTolerance(0.0, h / 2.0);
            RP1Solver.Integrate(rhsRP1, y0, rhsRP1.xUp, rhsRP1.xDown);
        }
    }
    rp1 = rhsRP1.xUp;
}

Double FTTFelectronicStates::boundEnergyRoot(Double e1, Double e2) {
    Double TFStates;
    Double trueStates;
    Double boundEnergyLeft = e1;
    Double boundEnergyRight = e2;
    Double currentBoundEnergy;

    while (abs(boundEnergyLeft - boundEnergyRight) / abs(boundEnergyLeft + boundEnergyRight) > eps/10)
    {
        currentBoundEnergy = (boundEnergyLeft + boundEnergyRight) / 2;

        trueStates = pseudoStates(currentBoundEnergy);
        TFStates = pseudoStatesTF(currentBoundEnergy);

        if (trueStates > TFStates) boundEnergyLeft += (boundEnergyRight - boundEnergyLeft) / 2;
        else boundEnergyRight -= (boundEnergyRight - boundEnergyLeft) / 2;

    }
    currentBoundEnergy = (boundEnergyLeft + boundEnergyRight) / 2;
    return currentBoundEnergy;
}

void FTTFelectronicStates::calculateBoundEnergy() {
    for (Int i = 0; i < nMax - 1; ++i) {
        BERoot[i] = boundEnergyRoot(energyLevel[i + 1][0], energyLevel[i + 2][i + 1]);
        if (abs(BERoot[i] - energyLevel[i + 1][0]) < eps || 
            abs(BERoot[i] - energyLevel[i + 2][i + 1]) < eps)
            BERoot[i] = BERoot[i - 1];
        if (BERoot[i] > mixBound()) {
            boundEnergy = BERoot[i - 1];
            calculated_BE = true;
            return;
        }
    }
    boundEnergy = BERoot[nMax - 2];
    calculated_BE = true;
}

Double FTTFelectronicStates::mixBound() {
    for (Int n = 1; n <= nMax; ++n) {
        for (Int l = 1; l < n; ++l) {
            if (energyLevel[n][l] < energyLevel[n][l - 1]) {
                return energyLevel[n][n - 1];
            }
        }
    }
    return energyLevel[nMax][nMax - 1];
}

void FTTFelectronicStates::calculateEnergyLevel(Int n, Int l) {
    Double lambda;
    Double lArg;
    Double eArg;
    Double exactAction = M_PI*(n - l - 0.5);
    Double currentAction = 0;
    Double YnlLeft = startYleft;
    Double YnlRight = startYright;

    lambda = l + 0.5;
    lArg = 0.5*lambda*lambda / T() / r0 / r0;

    while (abs(YnlLeft - YnlRight) / abs(YnlLeft + YnlRight) > eps / 100) {
        eArg = (YnlRight + YnlLeft) / 2;
        calculateRP2(eArg, lArg);
        calculateRP1(eArg, lArg);
        calculateAction(eArg, lArg);
        if (action - exactAction > 0) {
            YnlRight -= (YnlRight - YnlLeft)*0.5;
        }
        else {
            YnlLeft += (YnlRight - YnlLeft)*0.5;
        }
    }

    energyLevel[n][l] = 0.5*(YnlRight + YnlLeft);
    calculatedLevel[n][l] = true;
    std::cout << "n = " << n <<
               "; l = " << l <<
             "; ynl = " << energyLevel[n][l] << std::endl;
}

void FTTFelectronicStates::calculateEnergyLevels() {
    for (int n = 1; n <= nMax; n++) {
        for (int l = 0; l < n; ++l) {
            if (!calculatedLevel[n][l]) calculateEnergyLevel(n, l);
        }
    }
}

void FTTFelectronicStates::calculateN(Double BE) {
    Double N = 0;
    Double Nhalf = 0;
    Double Nn;
    Double Nnl;

    for (Int n = 1; n <= nMax; ++n) {
        Nn = 0;
        for (Int l = 0; l < n; ++l) {
            if (energyLevel[n][l] < BE) {
                Nnl = (2 * l + 1) / (1 + exp(energyLevel[n][l]));
                Nn += Nnl;
            }
        }
        Nhalf += Nn;
    }
    N = Nhalf * 2;
    Nval = N;
    calculated_N = true;
}

void FTTFelectronicStates::calculateNTF(Double BE) {
    DoubleVec y0(0.0, RHSFTTFstates::dim);
    Double result;
    rhsFTTFstates.setBE(BE);

    y0[0] = phi.valueAt_1();
    y0[1] = y0[0];
    y0[2] = 0;

    FTTFstatesSolver.SetTolerance(1e-7, 0.0);
    FTTFstatesSolver.SetStep(1e-15);

    FTTFstatesSolver.Integrate(rhsFTTFstates, y0, 1.0, 0.0);
    result = y0[2];

    NTFval = result * 3.0 * T()*V()*sqrt(2.0 * T()) / M_PI / M_PI;
    calculated_NTF = true;
}

Double FTTFelectronicStates::pseudoStates(Double energy) {
    Double PS = 0;
    for (Int n = 1; n <= nMax; ++n) {
        for (Int l = 0; l < n; ++l) {
            if (energyLevel[n][l] < energy) PS += 2 * l + 1;
        }
    }
    return PS;
}

Double FTTFelectronicStates::pseudoStatesTF(Double energy) {
    DoubleVec y0(0.0, RHSBoundEnergy::dim);
    Double PS = 0;

    rhsBoundEnergy.setBE(energy);

    y0[0] = phi.valueAt_1();
    y0[1] = y0[0];
    y0[2] = 0;

    boundEnergySolver.SetTolerance(1e-8, 0.0);
    boundEnergySolver.SetStep(1e-15);
    boundEnergySolver.Integrate(rhsBoundEnergy, y0, 1.0, 0.0);
    PS = y0[2];

    PS = PS*1.5*T()*V()*sqrt(2.0 * T()) / M_PI / M_PI;

    return PS;
}

Double FTTFelectronicStates::N() {
    calculateEnergyLevels();
    if (!calculated_BE) calculateBoundEnergy();
    if (!calculated_N) calculateN(boundEnergy);
    return Nval;
}

Double FTTFelectronicStates::NTF() {
    calculateEnergyLevels();
    if (!calculated_BE) calculateBoundEnergy();
    if (!calculated_NTF) calculateNTF(boundEnergy);
    return NTFval;
}

Double FTTFelectronicStates::DN() {
    return N() - NTF();
}

Double FTTFelectronicStates::N(Double BE) {
    calculateEnergyLevels();
    if (!calculated_N) calculateN(BE);
    return Nval;
}

Double FTTFelectronicStates::NTF(Double BE) {
    if (!calculated_NTF) calculateNTF(boundEnergy);
    return NTFval;
}

Double FTTFelectronicStates::DN(Double BE) {
    return N(BE) - NTF(BE);
}