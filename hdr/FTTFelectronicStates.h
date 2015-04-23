#pragma once
#include "../lib/numeric/hdr/ODE/ODEsolver.h"
#include "../lib/numeric/hdr/ODE/ODEdata.h"
#include "../lib/numeric/hdr/ODE/steppers/ODEstepperPD853.h"
#include "../lib/numeric/hdr/numericTypes.h"
#include "../lib/numeric/hdr/mathUtils.h"
#include "../lib/numeric/hdr/SF/FermiDirac.h"
#include "../hdr/FTTFpotential.h"
#include "../hdr/ThermodynamicFunction.h"
#include <fstream>
#include <iostream>

class FTTFelectronicStates {
public:
	FTTFelectronicStates(const Int _nMax);

	void setParameters(const Volume &V, const Temperature &T, const Double Z);
	Double getDN();

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
		        dydx[1] = 1e+20;
		        dydx[2] = -1e+20;
		    }
		}
	private:
		Double a;
		Double energyArg;
		Double lambdaArg;
		Double p;
	};

	struct RHSRP2 {
		static const Int dim = 3;
		FermiDirac<Half> FDhalf;
		Double xUp;
		Double xDown;
		Double phi_xUp;
		Double dphi_xUp;
		Double dphi_xDown;
		Double phi_xDown;
		RHSRP2(Double _a = 0, Double eArg = 0, Double lArg = 0) : 
			a(_a), energyArg(eArg), lambdaArg(lArg) {}
		void setA(const Double _a) { a = _a; }
		void setLArg(const Double lArg) { lambdaArg = lArg; }
		void setEArg(const Double eArg) { energyArg = eArg; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {

		    dydx[0] = y[1];
		    if (x > 0) dydx[1] = a*x*FDhalf(y[0] / x);
		    else {
		        dydx[1] = 1e+20;
		    }

		    p2 = energyArg + y[0] / x - lambdaArg / x / x;

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
		static const Int dim = 3;
		FermiDirac<Half> FDhalf;
		RHSRP1(Double _a = 0, Double eArg = 0, Double lArg = 0) : 
			a(_a), energyArg(eArg), lambdaArg(lArg) {}
		void setA(const Double _a) { a = _a; }
		void setLArg(const Double lArg) { lambdaArg = lArg; }
		void setEArg(const Double eArg) { energyArg = eArg; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {

		    dydx[0] = y[1];
		    if (x > 0) dydx[1] = a*x*FDhalf(y[0] / x);
		    else {
		        dydx[1] = 1e+20;
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
		        dydx[1] = 1e+20;
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
		FermiDirac<HalfI> FDhalfI;
		RHSFTTFstates(Double _a = 0, Double _boundEnergy = 0) : 
			a(_a), boundEnergy(_boundEnergy) {}
		void setA(const Double _a) { a = _a; }
		void setBE(const Double BE) { boundEnergy = BE; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
		    dydx[0] = y[1];
		    if (x > 0) {
		        dydx[1] = a*x*FDhalf(y[0] / x);
		        if (boundEnergy + y[0] / x <= 0) dydx[2] = 0.0;
		        else    
		            dydx[2] = -x*x*FDhalfI(y[0] / x, 0 /**e10 + y[0] / t*/, boundEnergy + y[0] / x);
		    }
		    else {
		        dydx[1] = 1e+20;
		        dydx[2] = 0;
		    }
		}
	private:
		Double a;
		Double boundEnergy;	
	};

private:
	Double calculateAction(Double eArg, Double lArg);
	Double calculateRP2(Double eArg, Double lArg);
	Double calculateRP1(Double eArg, Double lArg);
	Double calculateDN();

	void calculateEnergyLevels();
	Double calculateEnergyLevel(Int n, Int l);
	Double calculateFTTFLowStates(Double boundEnergy);
	Double calculateExactLowStates(Double boundEnergy);
	Double calculateBoundEnergyRoot(Double e1, Double e2);
	Double calculateBoundEnergy();
	Double calculateMixBound();
	Double calculatePseudoStates(Double energy);
	Double calculatePseudoStatesTF(Double energy);

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

    Double action;

	Double rp1;
	Double rp2;
	Double phi_rp2;
	Double dphi_rp2;
        
    Double startYleft;
    Double startYright;

    Double r0;

    Volume V;
    Temperature T;
    Double dN;

    DoubleVec BERoot;
    DoubleMat energyLevel;
    Int nMax;
    Double eps;
};

FTTFelectronicStates::FTTFelectronicStates(const Int _nMax) : 
                nMax(_nMax), 
                actionData(-1),
                actionSolver(1e-6, 0.0),
                rhsAction(0.0),
                RP1Data(-1),
                RP1Solver(1e-6, 0.0),
                rhsRP1(0.0),
                RP2Data(-1),
                RP2Solver(1e-6, 0.0),
                rhsRP2(0.0),
                boundEnergyData(-1),
                boundEnergySolver(1e-6, 0.0),
                rhsBoundEnergy(0.0),
                FTTFstatesData(-1),
                FTTFstatesSolver(1e-6, 0.0),
                rhsFTTFstates(0.0),
                energyLevel(_nMax + 1, _nMax),
                BERoot(0.0, _nMax),
                eps(1e-7) {}

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
    }
}

Double FTTFelectronicStates::getDN() {return calculateDN(); }

Double FTTFelectronicStates::calculateAction(Double eArg, Double lArg) {
	DoubleVec y0(0.0, RHSAction::dim);
    Double result;

    rhsAction.setEArg(eArg);
    rhsAction.setLArg(lArg);

    y0[0] = phi_rp2;
    y0[1] = dphi_rp2;
    y0[2] = 0;

    actionSolver.SetOutput(actionData);
    actionSolver.SetTolerance(0.0, eps/10);
    actionSolver.Integrate(rhsAction, y0, rp2, rp1);

    result = y0[2];

    return action = r0*sqrt(2 * T())*result;
}

Double FTTFelectronicStates::calculateRP2(Double eArg, Double lArg) {
    DoubleVec y0(0.0, RHSRP2::dim);

    rhsRP2.setEArg(eArg);
    rhsRP2.setLArg(lArg);

    y0[0] = y0[1] = phi.valueAt_1();
    rhsRP2.xUp = 0;
    rhsRP2.xDown = 0;
    
    if (eArg + y0[0] - lArg > 0) {
        phi_rp2 = y0[0];
        dphi_rp2 = y0[1];
        return rp2 = 1.0;
    }
    else {
        RP2Solver.SetTolerance(0.0, 1e-10);
        RP2Solver.SetOutput(RP2Data);
        RP2Solver.Integrate(rhsRP2, y0, 1.0, 0.0);

        while (abs(rhsRP2.xUp - rhsRP2.xDown) > 1e-12) {

            y0[0] = rhsRP2.phi_xUp;
            y0[1] = rhsRP2.dphi_xUp;

            RP2Solver.SetStep((rhsRP2.xUp - rhsRP2.xDown) / 11.0);
            RP2Solver.SetTolerance(0.0, fabs(rhsRP2.xUp - rhsRP2.xDown) / 20.0);
            RP2Solver.Integrate(rhsRP2, y0, rhsRP2.xUp, rhsRP2.xUp - abs(rhsRP2.xUp - rhsRP2.xDown) / 11.0);
        }
        phi_rp2 = rhsRP2.phi_xDown;
        dphi_rp2 = rhsRP2.dphi_xDown;
        return rp2 = rhsRP2.xDown;
    }
}

Double FTTFelectronicStates::calculateRP1(Double eArg, Double lArg) {
    DoubleVec y0(0.0, RHSRP1::dim);

    rhsRP1.setEArg(eArg);
    rhsRP1.setLArg(lArg);

    y0[0] = y0[1] = phi.valueAt_1();
    rhsRP1.xUp = 0;
    rhsRP1.xDown = 0;
    
    RP1Solver.SetTolerance(0.0, 1e-10);
    RP1Solver.SetOutput(RP1Data);
    RP1Solver.Integrate(rhsRP1, y0, rp2, 0.0);

    while (abs(rhsRP1.xUp - rhsRP1.xDown) > 1e-12) {

        y0[0] = rhsRP1.phi_xUp;
        y0[1] = rhsRP1.dphi_xUp;

        RP1Solver.SetStep((rhsRP1.xUp - rhsRP1.xDown) / 11.0);
        RP1Solver.SetTolerance(0.0, fabs(rhsRP1.xUp - rhsRP1.xDown) / 20.0);
        RP1Solver.Integrate(rhsRP1, y0, rhsRP1.xUp, rhsRP1.xUp - abs(rhsRP1.xUp - rhsRP1.xDown) / 11.0);
    }

    return rp1 = rhsRP1.xUp;
}

Double FTTFelectronicStates::calculateFTTFLowStates(Double boundEnergy) {
    DoubleVec y0(0.0, RHSFTTFstates::dim);
    Double result;
    rhsFTTFstates.setBE(boundEnergy);

    y0[0] = phi.valueAt_1();
    y0[1] = y0[0];
    y0[2] = 0;

    FTTFstatesSolver.SetOutput(FTTFstatesData);
    FTTFstatesSolver.SetTolerance(1e-7, 0.0);
    FTTFstatesSolver.SetStep(1e-15);

    FTTFstatesSolver.Integrate(rhsFTTFstates, y0, 1.0, 0.0);
    result = y0[2];

    return result * 3.0 * T()*V()*sqrt(2.0 * T()) / M_PI / M_PI;
}

Double FTTFelectronicStates::calculateBoundEnergyRoot(Double e1, Double e2) {
    Double TFStates;
    Double trueStates;
    Double boundEnergyLeft = e1;
    Double boundEnergyRight = e2;
    Double currentBoundEnergy;

    while (abs(boundEnergyLeft - boundEnergyRight) / abs(boundEnergyLeft + boundEnergyRight) > eps/10)
    {
        currentBoundEnergy = (boundEnergyLeft + boundEnergyRight) / 2;

        trueStates = calculatePseudoStates(currentBoundEnergy);
        TFStates = calculatePseudoStatesTF(currentBoundEnergy);

        if (trueStates > TFStates) boundEnergyLeft += (boundEnergyRight - boundEnergyLeft) / 2;
        else boundEnergyRight -= (boundEnergyRight - boundEnergyLeft) / 2;

    }
    currentBoundEnergy = (boundEnergyLeft + boundEnergyRight) / 2;
    return currentBoundEnergy;
}

Double FTTFelectronicStates::calculateBoundEnergy() {
    Double mixBound = calculateMixBound();
    for (Int i = 0; i < nMax - 1; ++i) {
        BERoot[i] = calculateBoundEnergyRoot(energyLevel[i + 1][0], energyLevel[i + 2][i + 1]);
        if (abs(BERoot[i] - energyLevel[i + 1][0]) < eps || 
            abs(BERoot[i] - energyLevel[i + 2][i + 1]) < eps)
            BERoot[i] = BERoot[i - 1];
        if (BERoot[i] > mixBound) {
            return BERoot[i - 1];
        }
    }
    return BERoot[nMax - 2];
}

Double FTTFelectronicStates::calculateMixBound() {
    for (Int n = 1; n <= nMax; ++n) {
        for (Int l = 1; l < n; ++l) {
            if (energyLevel[n][l] < energyLevel[n][l - 1]) {
                return energyLevel[n][n - 1];
            }
        }
    }
    return energyLevel[nMax][nMax - 1];
}

Double FTTFelectronicStates::calculateEnergyLevel(Int n, Int l) {
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
        currentAction = calculateAction(eArg, lArg);
        if (currentAction - exactAction > 0) {
            YnlRight -= (YnlRight - YnlLeft)*0.5;
        }
        else {
            YnlLeft += (YnlRight - YnlLeft)*0.5;
        }
    }

    return energyLevel[n][l] = 0.5*(YnlRight + YnlLeft);
}

void FTTFelectronicStates::calculateEnergyLevels() {
    for (int n = 1; n <= nMax; n++) {
        for (int l = 0; l < n; ++l) {
            calculateEnergyLevel(n, l);
            std::cout << "n = " << n <<
                    "; l = " << l <<
                    "; ynl = " << energyLevel[n][l] << std::endl;
        }
    }
}

Double FTTFelectronicStates::calculateExactLowStates(Double boundEnergy) {
    Double N = 0;
    Double Nhalf = 0;
    Double Nn;
    Double Nnl;

    for (Int n = 1; n <= nMax; ++n) {
        Nn = 0;
        for (Int l = 0; l < n; ++l) {
            if (energyLevel[n][l] < boundEnergy) {
                Nnl = (2 * l + 1) / (1 + exp(energyLevel[n][l]));
                Nn += Nnl;
            }
        }
        Nhalf += Nn;
    }
    N = Nhalf * 2;

    return N;
}

Double FTTFelectronicStates::calculatePseudoStates(Double energy) {
    Double pseudoStates = 0;
    for (Int n = 1; n <= nMax; ++n) {
        for (Int l = 0; l < n; ++l) {
            if (energyLevel[n][l] < energy) pseudoStates += 2 * l + 1;
        }
    }
    return pseudoStates;
}

Double FTTFelectronicStates::calculatePseudoStatesTF(Double energy) {
    DoubleVec y0(0.0, RHSBoundEnergy::dim);
    Double pseudoStates = 0;

    rhsBoundEnergy.setBE(energy);

    y0[0] = phi.valueAt_1();
    y0[1] = y0[0];
    y0[2] = 0;

    boundEnergySolver.SetTolerance(1e-8, 0.0);
    boundEnergySolver.SetStep(1e-15);
    boundEnergySolver.Integrate(rhsBoundEnergy, y0, 1.0, 0.0);
    pseudoStates = y0[2];

    pseudoStates = pseudoStates*1.5*T()*V()*sqrt(2.0 * T()) / M_PI / M_PI;

    return pseudoStates;
}

Double FTTFelectronicStates::calculateDN() {
    Double N = 0;
    Double NTF = 0;
    Double boundEnergy;
    Double CP = phi.valueAt_1()*T();

    calculateEnergyLevels();
    boundEnergy = calculateBoundEnergy();

    std::cout << "EB = " << boundEnergy << std::endl;
    std::cout << "p = " << -phi.valueAt_1() << std::endl;

    N = calculateExactLowStates(boundEnergy);
    NTF = calculateFTTFLowStates(boundEnergy);
    
    std::cout << "N = " << N << std::endl;
    std::cout << "NTF = " << NTF << std::endl;

    return dN = N - NTF;
}
