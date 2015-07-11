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
#include "../hdr/Printer.h"
#include "../hdr/Timer.h"
#include <fstream>
#include <iostream>

class FTTFelectronicStates {
public:
	FTTFelectronicStates(const Int _nMax);
	void setParameters(const PhysQ &V, const PhysQ &T, const Double Z);
    void setPhiShift(const Double _phiShift = 0.0);
    // energy range for finding levels
    void setEnergyRange(const Double Yleft, const Double Yright);
    // get energy level n,l
    Double operator() (const Int n, const Int l);
    // number of states with user BE
    Double NTFlow(Double BE);
    Double Nlow(Double BE);
    Double DNlow(Double BE);
    Double N(Double BE); // Nlow + NTFhigh
    // number of states with optimal BE
    Double NTFlow();
    Double Nlow();
	Double DNlow();
    Double N(); // Nlow + NTFhigh
    // explicit calculation of energy levels
    void calculateEnergyLevels(); 
    /**
    * @brief Set tolerance eps for the further calculations. Default is @f$ 10^{-6} @f$.
    */
    void setTolerance(const Double _eps);
    // log and output
    /**
    * @brief Write log to specified stream.
    */
    void setMainLogStream(std::ofstream* _LOG);
    /**
    * @brief Disable writing log to specified stream.
    */
    void clearMainLogStream();
    /**
    * @brief Enable self-printing log output into file.
    */
    void setPrintMainLogOn();
    /**
    * @brief Disable self-printing log output into file.
    */
    void setPrintMainLogOff();
    /**
    * @brief Enable self-printing log with single level calculation into file.
    */
    void setPrintLevelLogOn();
    /**
    * @brief Disable self-printing log with single level calculation into file.
    */
    void setPrintLevelLogOff();

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
			a(_a), BE(_boundEnergy) {}
		void setA(const Double _a) { a = _a; }
		void setBE(const Double _BE) { BE = _BE; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
		    dydx[0] = y[1];
		    if (x > 0) {
		        dydx[1] = a*x*FDhalf(y[0] / x);
		        if (BE + y[0] / x <= 0) dydx[2] = 0.0;
		        else {
					dydx[2] = -2.0 / 3.0*x*x*pow(BE + y[0] / x, 3.0 / 2.0);
		        }
		    }
		    else {
		        dydx[1] = 1e+10;
		        dydx[2] = 0;
    		}
		}
	private:
		Double a;
		Double BE;
	};

	struct RHSFTTFstatesLow {
		static const Int dim = 3;
		FermiDirac<Half> FDhalf;
		FDhalfInc FDhalfI;
		RHSFTTFstatesLow(Double _a = 0, Double _BE = 0, Double _phiShift = 0) : 
			a(_a), BE(_BE), dphi(_phiShift), FDhalfI(1e-7) {}
		void setA(const Double _a) { a = _a; }
		void setBE(const Double _BE) { BE = _BE; }
        void setShift(const Double _phiShift) { dphi = _phiShift; }
		void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
		    dydx[0] = y[1];
		    if (x > 0) {
		        dydx[1] = a*x*FDhalf(y[0] / x + dphi);
		        if (BE + y[0] / x <= 0) dydx[2] = 0.0;
		        else    
		            dydx[2] = -x*x*FDhalfI(y[0] / x + dphi, 0, BE + y[0] / x);
		    }
		    else {
		        dydx[1] = 1e+10;
		        dydx[2] = 0;
		    }
		}
	private:
		Double a;
		Double BE;	
        Double dphi;
	};

    struct RHSFTTFstatesFull {
        static const Int dim = 3;
        FermiDirac<Half> FDhalf;
        RHSFTTFstatesFull(Double _a = 0, Double _phiShift = 0) : 
            a(_a), dphi(_phiShift) {}
        void setA(const Double _a) { a = _a; }
        void setShift(const Double _phiShift) { dphi = _phiShift; }
        void operator() (const Double x, DoubleVec &y, DoubleVec &dydx) {
            dydx[0] = y[1];
            if (x > 0) {
                dydx[1] = a*x*FDhalf(y[0] / x + dphi);
                dydx[2] = -x*x*FDhalf(y[0] / x + dphi);
            }
            else {
                dydx[1] = 1e+10;
                dydx[2] = 0;
            }
        }
    private:
        Double a;
        Double dphi;
    };

private:
	void calculateAction(Double eArg, Double lArg); Double action;
	void calculateRP2(Double eArg, Double lArg);    Double rp2; Double phi_rp2; Double dphi_rp2;
	void calculateRP1(Double eArg, Double lArg);    Double rp1; 

	void calculateNTFlow(Double boundEnergy);  Double NTFlowVal;  bool calculated_NTFlow;
    void calculateNTFfull();                   Double NTFfullVal; bool calculated_NTFfull;
    void calculateNlow(Double boundEnergy);    Double NlowVal;    bool calculated_Nlow;

	void calculateBoundEnergy(); Double boundEnergy; bool calculated_BE;
    void calculateEnergyLevel(Int n, Int l); DoubleMat energyLevel; bool** calculatedLevel;


    Double boundEnergyRoot(Double e1, Double e2); DoubleVec BERoot;
	Double MixBound();
	Double pseudoStates(Double energy);  
	Double pseudoStatesTF(Double energy);

    void reset();

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

    RHSFTTFstatesLow rhsFTTFstatesLow;
    ODEsolver<ODEstepperPD853<RHSFTTFstatesLow> > FTTFstatesLowSolver;
    ODEdata FTTFstatesLowData;

    RHSFTTFstatesFull rhsFTTFstatesFull;
    ODEsolver<ODEstepperPD853<RHSFTTFstatesFull> > FTTFstatesFullSolver;
    ODEdata FTTFstatesFullData;

    FTTFpotential phi;
        
    Double startYleft;
    Double startYright;

    Double r0;

    PhysQ V;
    PhysQ T;

    Int nMax;
    Double eps; 
    Double phiShift;

    // log and output
    Printer printer;
    Timer mainTimer;
    Timer levelTimer;
    bool printMainLogOn;
    bool printLevelLogOn;
    bool mainLogStreamIsSet;
    bool levelLogStreamIsSet;
    std::ofstream* mainLOG;
    std::ofstream* levelLOG;
    std::ofstream* OUT;
    Int precision;
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
    FTTFstatesLowData(-1),
    FTTFstatesLowSolver(1e-6),
    rhsFTTFstatesLow(0.0),
    FTTFstatesFullData(-1),
    FTTFstatesFullSolver(1e-6),
    rhsFTTFstatesFull(0.0),
    energyLevel(_nMax + 1, _nMax),
    BERoot(0.0, _nMax),
    eps(1e-7),
    phiShift(0.0),
    startYleft(-1000.0),
    startYright(100.0) {
        actionSolver.SetOutput(actionData);
        RP1Solver.SetOutput(RP1Data);
        RP2Solver.SetOutput(RP2Data);
        boundEnergySolver.SetOutput(boundEnergyData);
        FTTFstatesLowSolver.SetOutput(FTTFstatesLowData);
        FTTFstatesFullSolver.SetOutput(FTTFstatesFullData);
        actionSolver.SetTolerance(1e-10, 0.0);
        calculatedLevel = new bool*[_nMax + 1];
        for (Int i = 0; i <= _nMax; ++i) {
            calculatedLevel[i] = new bool[_nMax];
            for (Int j = 0; j < _nMax; ++j) {
                calculatedLevel[i][j] = false;
            }
        }
        mainLogStreamIsSet = false;
        levelLogStreamIsSet = false;
        printMainLogOn = false;
        printLevelLogOn = false;
        OUT = NULL;
        mainLOG = NULL;
        levelLOG = NULL;
        precision = 7;
    }

void FTTFelectronicStates::setTolerance(const Double _eps) {
    eps = _eps;
    precision = static_cast<int>(-log10(eps));
    if (printMainLogOn) {
        *mainLOG << "electronic states accepted new tolerance, eps = ";
        printer.printSciDouble(*mainLOG, eps, precision);
        *mainLOG << std::endl;
    }
    phi.setTolerance(eps);
}

void FTTFelectronicStates::setMainLogStream(std::ofstream* _LOG) {
    if (printMainLogOn) setPrintMainLogOff();
    mainLOG = _LOG;
    mainLogStreamIsSet = true;
    setPrintMainLogOn();
}

void FTTFelectronicStates::clearMainLogStream() {
    setPrintMainLogOff();
    if (mainLogStreamIsSet) {
        mainLOG = NULL;
        mainLogStreamIsSet = false;
    }
}

void FTTFelectronicStates::setPrintMainLogOn() {
    if (!mainLogStreamIsSet) {
        mainLOG = new std::ofstream;
    }
    mainTimer.start();
    printMainLogOn = true;
}

void FTTFelectronicStates::setPrintMainLogOff() {
    if (!mainLogStreamIsSet) {
        if (printMainLogOn) {
            if (mainLOG->is_open()) mainLOG->close();
            delete mainLOG;
            mainLOG = NULL;
        }
    }
    mainTimer.stop();
    printMainLogOn = false;
}

void FTTFelectronicStates::setPrintLevelLogOn() {
    levelLOG = new std::ofstream;
    levelLogStreamIsSet = true;
    printLevelLogOn = true;
    levelTimer.start();
}

void FTTFelectronicStates::setPrintLevelLogOff() {
    if (printLevelLogOn) {
        if (levelLOG->is_open()) levelLOG->close();
        delete levelLOG;
        levelLOG = NULL;
        levelLogStreamIsSet = false;
    }
    printLevelLogOn = false;
    levelTimer.stop();
}

void FTTFelectronicStates::setParameters(const PhysQ &_V, const PhysQ &_T, const Double Z) {
    
    if (abs(log10(_V()) - log10(V())) > 1e-10
        || abs(log10(_T()) - log10(T())) > 1e-10) 
    {
    	V = _V;
    	T = _T;

        if (printMainLogOn) {
            if (!mainLogStreamIsSet) {
                std::stringstream filename;
                filename << "log/log_eStates;";
                filename << "_Z(" << Z << ")"; 
                filename << "_V(" << V() << ")_";
                filename << "_T(" << T() << ")_";
                filename << mainTimer.getCurrentDatetime();
                filename << ".txt";
                mainLOG->open(filename.str().c_str(), std::ios::out);
            }
            mainTimer.start();
        }

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
        rhsFTTFstatesLow.setA(a);
        rhsFTTFstatesFull.setA(a);
        rhsFTTFstatesLow.setShift(0.0);
        rhsFTTFstatesFull.setShift(0.0);
        r0 = pow(3.0*V() / 4.0 / M_PI, 1.0 / 3.0);
        reset();
        /********************main LOG section******************************************/
        if (printMainLogOn) {                                                         //
            (*mainLOG) <<                                                             //
                "FTTF electronic states accepted new parameters:"                     //
                << std::endl;                                                         //
            printer.printString((*mainLOG), "V_{Z}[Atomic] = ", 22, right);           //
            printer.printSciDouble((*mainLOG), V(), precision, 15, left);             //
            (*mainLOG) << std::endl;                                                  //
            printer.printString((*mainLOG), "T_{Z}[Hartree] = ", 23, right);          //
            printer.printSciDouble((*mainLOG), T(), precision, 15, left);             //
            (*mainLOG) << std::endl;                                                  //
            printer.printString((*mainLOG), "Z = ", 10, right);                       //
            printer.printDouble((*mainLOG), Z, 3, 15, left);                          //
            (*mainLOG) << std::endl;                                                  //
            (*mainLOG) << "Recalculate parameters:" << std::endl;                     //
            printer.printString((*mainLOG),                                           //
                                "V_{1}[Atomic] = V_{Z}*Z = ",                         //
                                32, right);                                           //
            printer.printSciDouble((*mainLOG), V()*Z, precision, 15, left);           //
            (*mainLOG) << std::endl;                                                  //
            printer.printString((*mainLOG),                                           //
                                "T_{1}[Hartree] = T_{Z}*Z^{-4/3} = ",                 //
                                40, right);                                           //
            printer.printSciDouble((*mainLOG), T()*pow(Z, -4.0/3.0),                  //
                                    precision, 15, left);                             //
            (*mainLOG) << std::endl;                                                  //
            printer.printString((*mainLOG),                                           //
                                "a = 2^{7/6}*3^{2/3}*\\pi^{-5/3}*T^{1/2}*V^{2/3} = ", //
                                 55, right);                                          //
            printer.printSciDouble((*mainLOG), a, precision, 15, left);               //
            (*mainLOG) << std::endl;                                                  //
        }                                                                             //
        /******************************************************************************/
    }
}

void FTTFelectronicStates::reset() {
    for (Int i = 0; i <= nMax; ++i) {
        calculatedLevel[i] = new bool[nMax];
        for (Int j = 0; j < nMax; ++j) {
            calculatedLevel[i][j] = false;
        }
    }
    calculated_BE = false;
    calculated_NTFlow = false;
    calculated_NTFfull = false;
    calculated_Nlow = false;
    /********************main LOG section*********************************/
    if (printMainLogOn) {                                                //
        (*mainLOG) << "FTTF electron states data cleared." << std::endl; //
    }                                                                    //
    /*********************************************************************/
}

void FTTFelectronicStates::setPhiShift(const Double _phiShift) {
    rhsFTTFstatesLow.setShift(_phiShift);
    rhsFTTFstatesFull.setShift(_phiShift);
    phiShift = _phiShift;
    calculated_NTFlow = false;
    calculated_NTFfull = false;
    calculated_Nlow = false;

    /********************main LOG section*****************************************/
    if (printMainLogOn) {                                                        //
        (*mainLOG) << "Potential boundary condition was shifted: " << std::endl; //
        (*mainLOG) << "\\delta\\phi(1) = ";                                      //
        printer.printSciDouble(*mainLOG, _phiShift, precision, 15, left);        //
        (*mainLOG) << std::endl;                                                 //
    }                                                                            //
    /*****************************************************************************/
}

void FTTFelectronicStates::setEnergyRange(const Double Yleft, const Double Yright) {
    startYleft = Yleft;
    startYright = Yright;
    /********************main LOG section*****************************************/
    if (printMainLogOn) {                                                        //
        (*mainLOG) << "Energy range for calculating energy levels updated: "     //
                   << std::endl;                                                 //
        (*mainLOG) << "y = (\\varepsilon - \\mu)/T" << std::endl;                //
        (*mainLOG) << "y^{left} = ";                                             //
        printer.printSciDouble(*mainLOG, Yleft, precision, 15, left);            //
        (*mainLOG) << std::endl;                                                 //
        (*mainLOG) << "y^{left} = ";                                             //
        printer.printSciDouble(*mainLOG, Yright, precision, 15, left);           //
        (*mainLOG) << std::endl;                                                 //
    }                                                                            //
    /*****************************************************************************/
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
    Double actTimeStart;
    Double actTimeEnd;
    /************************level LOG section***************************/
    if (printLevelLogOn) {                                              //
        (*levelLOG) << "Calculating action:" << std::endl;              //
        (*levelLOG) << "energy parameter: ";                            //
        (*levelLOG) << " y^{test} = (\\varepsilon^{test} - \\mu)/T = "; //
        printer.printSciDouble(*levelLOG, eArg, precision, 15);         //
        (*levelLOG) << "lambda parameter:";                             //
        (*levelLOG) << "\\lambda = l + 1/2" << std::endl;               //
        (*levelLOG) << "parameter = \\lambda^2/(2T r_0^2) = ";          //
        printer.printSciDouble(*levelLOG, lArg, precision, 15);         //
        (*levelLOG) << std::endl;                                       //
        (*levelLOG) << "rotate points:" << std::endl;                   //
        (*levelLOG) << "x_1 = ";                                        //
        printer.printSciDouble(*levelLOG, rp1, precision, 15);          //
        (*levelLOG) << std::endl;                                       //
        (*levelLOG) << "x_2 = ";                                        //
        printer.printSciDouble(*levelLOG, rp2, precision, 15);          //
        (*levelLOG) << std::endl;                                       //
        actTimeStart = levelTimer.getElapsedTimeInMilliSec();           //
    }                                                                   //
    /********************************************************************/
    if (rp2 > rp1) {
        rhsAction.setEArg(eArg);
        rhsAction.setLArg(lArg);
        y0[0] = phi_rp2;
        y0[1] = dphi_rp2;
        y0[2] = 0.0;
        actionSolver.Integrate(rhsAction, y0, rp2, rp1);
        action = r0*sqrt(2 * T())*y0[2];
        actTimeEnd = levelTimer.getElapsedTimeInMilliSec();
        /***************************level LOG section************************************/
        if (printLevelLogOn) {                                                          //
            actTimeEnd = levelTimer.getElapsedTimeInMilliSec();                         //
            (*levelLOG) << "output steps during integration:" << std::endl;             //
            Int size = actionData.Count();                                              //
            printer.printString((*levelLOG), "x");                                      //
            printer.printString((*levelLOG), "phi");                                    //
            printer.printString((*levelLOG), "dphi");                                   //
            printer.printString((*levelLOG), "S(x)/(r_0*(2T)^{1/2})");                  //
            (*levelLOG) << std::endl;                                                   //
            for (Int i = 0; i < size; ++i) {                                            //
                printer.printSciDouble((*levelLOG), actionData.xSave[i], precision);    //
                printer.printSciDouble((*levelLOG), actionData.ySave[0][i], precision); //
                printer.printSciDouble((*levelLOG), actionData.ySave[1][i], precision); //
                printer.printSciDouble((*levelLOG), actionData.ySave[2][i], precision); //
                (*levelLOG) << std::endl;                                               //
            }                                                                           //
            (*levelLOG) << std::endl;                                                   //
            (*levelLOG) << "S = S(0) = ";                                               //
            printer.printSciDouble((*levelLOG), action, precision, 15, left);           //
            (*levelLOG) << std::endl;                                                   //
            (*levelLOG) << "Elapsed time: ";                                            //
            printer.printSciDouble((*levelLOG), actTimeEnd - actTimeStart,              //
                                    precision, 15, left);                               //
            (*levelLOG) << " ms" << std::endl;                                          //
        }                                                                               //
        /********************************************************************************/
    }
    else {
        action = 0;
        /***************************level LOG section****************************/
        if (printLevelLogOn) {                                                  //
            (*levelLOG) << "x_2 <= x_1: do not perform calculation, set S = 0"; //
            (*levelLOG) << std::endl;                                           //
        }                                                                       //
        /************************************************************************/   
    }
}

void FTTFelectronicStates::calculateRP2(Double eArg, Double lArg) {
    DoubleVec y0(0.0, RHSRP2::dim);
    Double h = 0.0;
    Double rp2TimeStart;
    Double rp2TimeEnd;
    Double p;
    Double error;
    Double epsNew;
    Int nSteps = 0;

    /***************************level LOG section*************************/
    if (printLevelLogOn) {                                               //
        (*levelLOG) << "Calculating right rotate point:" << std::endl;   //
        (*levelLOG) << "energy parameter: ";                             //
        (*levelLOG) << " y^{test} = (\\varepsilon^{test} - \\mu)/T = ";  //
        printer.printSciDouble(*levelLOG, eArg, precision, 15);          //
        (*levelLOG) << "lambda parameter:";                              //
        (*levelLOG) << "\\lambda = l + 1/2" << std::endl;                //
        (*levelLOG) << "parameter = \\lambda^2/(2T r_0^2) = ";           //
        printer.printSciDouble(*levelLOG, lArg, precision, 15);          //
        (*levelLOG) << std::endl;                                        //
        rp2TimeStart = levelTimer.getElapsedTimeInMilliSec();                 //
    }                                                                    //
    /*********************************************************************/
    rhsRP2.setEArg(eArg);
    rhsRP2.setLArg(lArg);
    
    y0[0] = y0[1] = phi(1);
    rhsRP2.xUp = 0;
    rhsRP2.xDown = 0;
    p = eArg + y0[0] - lArg;
    
    if (p > 0) {
        phi_rp2 = y0[0];
        dphi_rp2 = y0[1];
        rp2 = 1.0;
        /***************************level LOG section********************/
        if (printLevelLogOn) {                                          //
            (*levelLOG) << "Radial momentum p_r^2(1) =";                //
            printer.printSciDouble(*levelLOG, p, precision, 15);        //
            (*levelLOG) << " > 0, then set right rotate point x_2 = 1"; //
            (*levelLOG) << std::endl;                                   //
        }                                                               //
        /****************************************************************/
    }
    else {
        /***************************level LOG section**********************/
        if (printLevelLogOn) {                                            //
            (*levelLOG) << "Radial momentum p_r^2(1) =";                  //
            printer.printSciDouble(*levelLOG, p, precision, 15);          //
            (*levelLOG) << " < 0, then start searching rotate point x_2"; //
            (*levelLOG) << std::endl;                                     //
        }                                                                 //
        /******************************************************************/
        h = 1e-6;
        epsNew = 1e-10;
        RP2Solver.SetTolerance(0.0, epsNew);
        RP2Solver.Integrate(rhsRP2, y0, 1.0, 0.0);
        /***************************level LOG section************************/
        if (printLevelLogOn) {                                              //
            (*levelLOG) << "x_2 location detected:";                        //
            (*levelLOG) << std::endl;                                       //
            (*levelLOG) << "x_2^{down} = ";                                 //
            printer.printSciDouble(*levelLOG, rhsRP2.xDown, precision, 15); //
            (*levelLOG) << std::endl;                                       //
            (*levelLOG) << "x_2^{up} = ";                                   //
            printer.printSciDouble(*levelLOG, rhsRP2.xUp, precision, 15);   //
            (*levelLOG) << std::endl;                                       //
            (*levelLOG) << "Begin iteration procedure:";                    //
            printer.printString(*levelLOG, "nSteps");                       //
            printer.printString(*levelLOG, "xDown");                        //
            printer.printString(*levelLOG, "xUp");                          //
            printer.printString(*levelLOG, "\\phi(xUp)");                   //
            printer.printString(*levelLOG, "\\phi'(xUp)");                  //
            printer.printString(*levelLOG, "hStart");                       //
            printer.printString(*levelLOG, "epsRel");                       //
            printer.printString(*levelLOG, "error");                        //
            (*levelLOG) << std::endl;                                       //
        }                                                                   //
        /********************************************************************/
        error = abs(rhsRP2.xUp - rhsRP2.xDown);
        while (error > 1e-12) {
            y0[0] = rhsRP2.phi_xUp;
            y0[1] = rhsRP2.dphi_xUp;
            h =  (rhsRP2.xUp - rhsRP2.xDown) / 11.0;
            epsNew = h/2.0;
            nSteps++;
            /***********************level LOG section********************/
            if (printLevelLogOn) {                                      //
                printer.printInt(*levelLOG, nSteps);                    //
                printer.printSciDouble(*levelLOG, rhsRP2.xDown, 12);    //
                printer.printSciDouble(*levelLOG, rhsRP2.xUp, 12);      //
                printer.printSciDouble(*levelLOG, rhsRP2.phi_xUp, 12);  //
                printer.printSciDouble(*levelLOG, rhsRP2.dphi_xUp, 12); //
                printer.printSciDouble(*levelLOG, h, 12);               //
                printer.printSciDouble(*levelLOG, epsNew, 12);          //
                printer.printSciDouble(*levelLOG, error, 12);           //
                (*levelLOG) << std::endl;                               //
            }                                                           //
            /************************************************************/
            RP2Solver.SetStep(h);
            RP2Solver.SetTolerance(0.0, epsNew);
            RP2Solver.Integrate(rhsRP2, y0, rhsRP2.xUp, rhsRP2.xDown);
            error = abs(rhsRP2.xUp - rhsRP2.xDown);
        }
        phi_rp2 = rhsRP2.phi_xDown;
        dphi_rp2 = rhsRP2.dphi_xDown;
        rp2 = rhsRP2.xDown;
        /***************************level LOG section****************/
        if (printLevelLogOn) {                         //
            rp2TimeEnd = levelTimer.getElapsedTimeInMilliSec();     //
            (*levelLOG) << "Finally calculated value x_2 = ";       //
            printer.printSciDouble(*levelLOG, rp2, 12);             //
            (*levelLOG) << std::endl;                               //
            (*levelLOG) << "Elapsed time:";                         //
            printer.printSciDouble(*levelLOG,                       //
                                    rp2TimeEnd - rp2TimeStart, 12); //
            (*levelLOG) << " ms" << std::endl;                      //
        }                                                           //
        /************************************************************/
    }
}

void FTTFelectronicStates::calculateRP1(Double eArg, Double lArg) {
    DoubleVec y0(0.0, RHSRP1::dim);
    Double h = 0.0;
    Double rp1TimeStart;
    Double rp1TimeEnd;
    Double error;
    Double epsNew;
    Int nSteps = 0;
    /***************************level LOG section**********************/
    if (printLevelLogOn) {                                            //
        (*levelLOG) << "Calculating left rotate point:" << std::endl; //
        rp1TimeStart = levelTimer.getElapsedTimeInMilliSec();         //
    }                                                                 //
    /******************************************************************/

    if (rp2 > 1e-12) {
        rhsRP1.setEArg(eArg);
        rhsRP1.setLArg(lArg);

        y0[0] = phi_rp2;
        y0[1] = dphi_rp2;

        rhsRP1.xUp = 0;
        rhsRP1.xDown = 0;
        
        RP1Solver.SetTolerance(0.0, 1e-10);
        RP1Solver.Integrate(rhsRP1, y0, rp2, 0.0);
        /***************************level LOG section************************/
        if (printLevelLogOn) {                                              //
            (*levelLOG) << "x_1 location detected:";                        //
            (*levelLOG) << std::endl;                                       //
            (*levelLOG) << "x_1^{down} = ";                                 //
            printer.printSciDouble(*levelLOG, rhsRP1.xDown, precision, 15); //
            (*levelLOG) << std::endl;                                       //
            (*levelLOG) << "x_1^{up} = ";                                   //
            printer.printSciDouble(*levelLOG, rhsRP1.xUp, precision, 15);   //
            (*levelLOG) << std::endl;                                       //
            (*levelLOG) << "Begin iteration procedure:";                    //
            printer.printString(*levelLOG, "nSteps");                       //
            printer.printString(*levelLOG, "xDown");                        //
            printer.printString(*levelLOG, "xUp");                          //
            printer.printString(*levelLOG, "\\phi(xUp)");                   //
            printer.printString(*levelLOG, "\\phi'(xUp)");                  //
            printer.printString(*levelLOG, "hStart");                       //
            printer.printString(*levelLOG, "epsRel");                       //
            printer.printString(*levelLOG, "error");                        //
            (*levelLOG) << std::endl;                                       //
        }                                                                   //
        /********************************************************************/
        error = abs(rhsRP1.xUp - rhsRP1.xDown);
        while (error > 1e-12) {
            y0[0] = rhsRP1.phi_xUp;
            y0[1] = rhsRP1.dphi_xUp;
            h = (rhsRP1.xUp - rhsRP1.xDown) / 11.0;
            epsNew = h/2.0;
            nSteps++;
            /***********************level LOG section********************/
            if (printLevelLogOn) {                                      //
                printer.printInt(*levelLOG, nSteps);                    //
                printer.printSciDouble(*levelLOG, rhsRP1.xDown, 12);    //
                printer.printSciDouble(*levelLOG, rhsRP1.xUp, 12);      //
                printer.printSciDouble(*levelLOG, rhsRP1.phi_xUp, 12);  //
                printer.printSciDouble(*levelLOG, rhsRP1.dphi_xUp, 12); //
                printer.printSciDouble(*levelLOG, h, 12);               //
                printer.printSciDouble(*levelLOG, epsNew, 12);          //
                printer.printSciDouble(*levelLOG, error, 12);           //
                (*levelLOG) << std::endl;                               //
            }                                                           //
            /************************************************************/
            RP1Solver.SetStep(h);
            RP1Solver.SetTolerance(0.0, epsNew);
            RP1Solver.Integrate(rhsRP1, y0, rhsRP1.xUp, rhsRP1.xDown);
            error = abs(rhsRP1.xUp - rhsRP1.xDown);
        }
        rp1 = rhsRP1.xUp;
        /***********************level LOG section********************/
        if (printLevelLogOn) {                                      //
            rp1TimeEnd = levelTimer.getElapsedTimeInMilliSec();     //
            (*levelLOG) << "Finally calculated value x_1 = ";       //
            printer.printSciDouble(*levelLOG, rp1, 12);             //
            (*levelLOG) << std::endl;                               //
            (*levelLOG) << "Elapsed time:";                         //
            printer.printSciDouble(*levelLOG,                       //
                                    rp1TimeEnd - rp1TimeStart, 12); //
            (*levelLOG) << " ms" << std::endl;                      //
        }                                                           //
        /************************************************************/
    }
    else {
        rp1 = 0;
        /***********************level LOG section******************/
        if (printLevelLogOn) {                                    //
            (*levelLOG) << "Right rotate point x_2 = ";           //
            printer.printSciDouble(*levelLOG, rp2, 12, 20, left); //
            (*levelLOG) << std::endl;                             //
            (*levelLOG) << "Set left rotate point x_1 = 0";       //
            (*levelLOG) << std::endl;                             //
        }                                                         //
        /**********************************************************/   
    }

}

Double FTTFelectronicStates::boundEnergyRoot(Double e1, Double e2) {
    Double TFStates;
    Double trueStates;
    Double boundEnergyLeft = e1;
    Double boundEnergyRight = e2;
    Double currentBoundEnergy;
    Int nStep = 0;
    DoubleVec timeSteps(0.0, 100);
    Double timePerStep;
    Double overallTime;
    Double averageTime;
    Double err;
    /************************main LOG section*********************/
    if (printMainLogOn) {                                        //
        printer.printString((*mainLOG), "nStep",      10, left); //
        printer.printString((*mainLOG), "BEleft",     20, left); //
        printer.printString((*mainLOG), "BEright",    20, left); //
        printer.printString((*mainLOG), "BEcenter",   20, left); //
        printer.printString((*mainLOG), "Nstates",    20, left); //
        printer.printString((*mainLOG), "NstatesTF",  20, left); //
        printer.printString((*mainLOG), "error",      20, left); //
        printer.printString((*mainLOG), "time[ms]",   20, left); //
        (*mainLOG) << std::endl;                                 //
    }                                                            //
    /*************************************************************/
    err = abs(boundEnergyLeft - boundEnergyRight) / abs(boundEnergyLeft + boundEnergyRight);
    while (err > eps/10) {
        /************************main LOG section*****************************************/
        if (printMainLogOn) {                                                            //
            timePerStep = mainTimer.getElapsedTimeInMilliSec();                          //
        }                                                                                //
        /*********************************************************************************/
        currentBoundEnergy = (boundEnergyLeft + boundEnergyRight) / 2;

        trueStates = pseudoStates(currentBoundEnergy);
        TFStates = pseudoStatesTF(currentBoundEnergy);
        /************************main LOG section*****************************************/
        if (printMainLogOn) {                                                            //
            printer.printInt((*mainLOG), nStep, 10, left);                               //
            printer.printSciDouble((*mainLOG), boundEnergyLeft, precision, 20, left);    //
            printer.printSciDouble((*mainLOG), boundEnergyRight, precision, 20, left);   //
            printer.printSciDouble((*mainLOG), currentBoundEnergy, precision, 20, left); //
            printer.printSciDouble((*mainLOG), trueStates, precision, 20, left);         //
            printer.printSciDouble((*mainLOG), TFStates, precision, 20, left);           //
        }                                                                                //
        /*********************************************************************************/
        if (trueStates > TFStates) boundEnergyLeft += (boundEnergyRight - boundEnergyLeft) / 2;
        else boundEnergyRight -= (boundEnergyRight - boundEnergyLeft) / 2;

        err = abs(boundEnergyLeft - boundEnergyRight) / abs(boundEnergyLeft + boundEnergyRight);
        /************************main LOG section*********************************/
        if (printMainLogOn) {                                                    //
            printer.printSciDouble((*mainLOG), err, precision, 20, left);        //
            timePerStep = mainTimer.getElapsedTimeInMilliSec() - timePerStep;    //
            timeSteps[nStep] = timePerStep;                                      //
            printer.printSciDouble((*mainLOG), timePerStep, precision, 20, left);//
            *mainLOG << std::endl;                                               //
            ++nStep;                                                             //
        }                                                                        //
        /*************************************************************************/
    }
    currentBoundEnergy = (boundEnergyLeft + boundEnergyRight) / 2;
    /************************main LOG section*****************************/
    if (printMainLogOn) {                                                //
        overallTime = 0;                                                 //
        for (Int i = 0; i <= nStep; ++i) {                               //
            overallTime += timeSteps[i];                                 //
        }                                                                //
        averageTime = overallTime/(nStep + 1);                           //
        *mainLOG << "calculated root: "                                  //
                 << currentBoundEnergy << std::endl;                     //
        *mainLOG << "overall time: "                                     //
                 << overallTime << " [ms]" << std::endl;                 //
        *mainLOG << "average time: "                                     //
                 << averageTime << " [ms]" << std::endl;                 //
    }                                                                    //
    /*********************************************************************/
    return currentBoundEnergy;
}

void FTTFelectronicStates::calculateBoundEnergy() {
    Double mixBound = MixBound();
    for (Int i = 0; i < nMax - 1; ++i) {
        /************************main LOG section*****************************/
        if (printMainLogOn) {                                                //
            *mainLOG <<                                                      //
                    "searching root for bound energy between energy levels:" //
                     << std::endl;                                           //
            *mainLOG << "e[n=" << i + 1 << "][l=" << 0 << "] = "             //
                     << energyLevel[i + 1][0] << std::endl;                  //
            *mainLOG << "e[n=" << i + 2 << "][l=" << i + 1 << "] = "         //
                     << energyLevel[i + 2][i + 1] << std::endl;              //
        }                                                                    //
        /*********************************************************************/
        // calculating next root for bound energy
        BERoot[i] = boundEnergyRoot(energyLevel[i + 1][0], 
                                    energyLevel[i + 2][i + 1]);
        // check whether the root is between initial values,
        // otherwise there is no root and we use previous root
        // at this step

        if (abs(BERoot[i] - energyLevel[i + 1][0]) < eps || 
            abs(BERoot[i] - energyLevel[i + 2][i + 1]) < eps)
            BERoot[i] = BERoot[i - 1];
        // check whether current root is less than boundary
        // of mixing energy levels
        if (BERoot[i] > MixBound()) {
            boundEnergy = BERoot[i - 1];
            /*********************************main LOG section*******************************************/
            if (printMainLogOn) {                                                                       //
                (*mainLOG) << "Mixing levels boundary: "                                                //
                           << mixBound << std::endl;                                                    //
                (*mainLOG) << "All found roots:" << std::endl;                                          //
                printer.printString((*mainLOG), "n",      10, left);                                    //
                printer.printString((*mainLOG), "e1",     20, left);                                    //
                printer.printString((*mainLOG), "e2",     20, left);                                    //
                printer.printString((*mainLOG), "BEroot", 20, left);                                    //
                (*mainLOG) << std::endl;                                                                //
                for (Int j = 0; j <= i; ++j) {                                                          //
                    printer.printInt((*mainLOG), j, 10, left);                                          //
                    printer.printSciDouble((*mainLOG), energyLevel[j + 1][0], precision, 20, left);     //
                    printer.printSciDouble((*mainLOG), energyLevel[j + 2][j + 1], precision, 20, left); //
                    printer.printSciDouble((*mainLOG), BERoot[j], precision, 20, left);                 //
                    (*mainLOG) << std::endl;                                                            //
                }                                                                                       //
                (*mainLOG) << "Selected root before mixing: "                                           //
                           << boundEnergy << std::endl;                                                 //
            }                                                                                           //
            /********************************************************************************************/
            calculated_BE = true;
            return;
        }
    }
    boundEnergy = BERoot[nMax - 2];
    /*********************************main LOG section*******************************************/
    if (printMainLogOn) {                                                                       //
        (*mainLOG) << "Mixing levels boundary: "                                                //
                   << mixBound << std::endl;                                                    //
        (*mainLOG) << "All found roots:" << std::endl;                                          //
        printer.printString((*mainLOG), "n",      10, left);                                    //
        printer.printString((*mainLOG), "e1",     20, left);                                    //
        printer.printString((*mainLOG), "e2",     20, left);                                    //
        printer.printString((*mainLOG), "BEroot", 20, left);                                    //
        for (Int j = 0; j < nMax; ++j) {                                                        //
            printer.printInt((*mainLOG), j, 10, left);                                          //
            printer.printSciDouble((*mainLOG), energyLevel[j + 1][0], precision, 20, left);     //
            printer.printSciDouble((*mainLOG), energyLevel[j + 2][j + 1], precision, 20, left); //
            printer.printSciDouble((*mainLOG), BERoot[j], precision, 20, left);                 //
            (*mainLOG) << std::endl;                                                            //
        }                                                                                       //
        (*mainLOG) << "Selected root before mixing: BE = "                                      //
                   << boundEnergy << std::endl;                                                 //
    }                                                                                           //
    /********************************************************************************************/
    calculated_BE = true;
}

Double FTTFelectronicStates::MixBound() {
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

    DoubleVec timeSteps(0.0, 100);
    Double timePerStep;
    Double overallTime;
    Double averageTime;
    Double err;
    Int nStep = 0;
    /**********main LOG section*********************/
    if (printMainLogOn) {                          //
        printer.printInt((*mainLOG), n, 10, left); //
        printer.printInt((*mainLOG), l, 10, left); //
    }                                              //
    /***********************************************/
    err = abs(YnlLeft - YnlRight) / abs(YnlLeft + YnlRight);
    while (err > eps / 100) {
        if (printMainLogOn) {
            timePerStep = mainTimer.getElapsedTimeInMilliSec();
        }
        eArg = (YnlRight + YnlLeft) / 2;
        calculateRP2(eArg, lArg);
        calculateRP1(eArg, lArg);
        calculateAction(eArg, lArg);
        /***************************main LOG section****************************/
        if (printMainLogOn && nStep > 0) {                                     //
            printer.printString((*mainLOG), "", 20, left, ' ');                //
        }                                                                      //
        if (printMainLogOn) {                                                  //
            printer.printInt((*mainLOG), nStep, 10, left);                     //
            printer.printSciDouble((*mainLOG), YnlLeft, precision, 20, left);  //
            printer.printSciDouble((*mainLOG), YnlRight, precision, 20, left); //
            printer.printSciDouble((*mainLOG), eArg, precision, 20, left);     //
            printer.printSciDouble((*mainLOG), action, precision, 20, left);   //
            printer.printSciDouble((*mainLOG), rp1, precision, 20, left);      //
            printer.printSciDouble((*mainLOG), rp2, precision, 20, left);      //
        }                                                                      //
        /***********************************************************************/
        if (action - exactAction > 0) {
            YnlRight -= (YnlRight - YnlLeft)*0.5;
        }
        else {
            YnlLeft += (YnlRight - YnlLeft)*0.5;
        }
        err = abs(YnlLeft - YnlRight) / abs(YnlLeft + YnlRight);
        /***************************main LOG section*******************************/
        if (printMainLogOn) {                                                     //
            timePerStep = mainTimer.getElapsedTimeInMilliSec() - timePerStep;     //
            timeSteps[nStep] = timePerStep;                                       //
            printer.printSciDouble((*mainLOG), err, precision, 20, left);         //
            printer.printSciDouble((*mainLOG), timePerStep, precision, 20, left); //
            *mainLOG << std::endl;                                                //
            ++nStep;                                                              //
        }                                                                         //
        /**************************************************************************/
    }

    energyLevel[n][l] = 0.5*(YnlRight + YnlLeft);
    calculatedLevel[n][l] = true;
    /********************************main LOG section********************************/
    if (printMainLogOn) {                                                           //
        printer.printString((*mainLOG), "", 20, left, ' ');                         //
        printer.printString((*mainLOG), "", 170, left, '-');                        //
        *mainLOG << std::endl;                                                      //
        printer.printString((*mainLOG), "Final results: ",  20, left, ' ');         //
        printer.printString((*mainLOG), "Ynl",              20, left);              //
        printer.printString((*mainLOG), "Sexact",           20, left);              //
        printer.printString((*mainLOG), "Scalc",            20, left);              //
        printer.printString((*mainLOG), "rp1",              20, left);              //
        printer.printString((*mainLOG), "rp2",              20, left);              //
        printer.printString((*mainLOG), "time[ms]",         20, left);              //
        printer.printString((*mainLOG), "timePerIter[ms]",  20, left);              //
        *mainLOG << std::endl;                                                      //
        printer.printString((*mainLOG), "", 20, left, ' ');                         //
        printer.printString((*mainLOG), "", 170, left, '-');                        //
        *mainLOG << std::endl;                                                      //
        printer.printString((*mainLOG), "",  20, left, ' ');                        //
        printer.printSciDouble((*mainLOG), energyLevel[n][l], precision, 20, left); //
        printer.printSciDouble((*mainLOG), exactAction,       precision, 20, left); //
        printer.printSciDouble((*mainLOG), action,            precision, 20, left); //
        printer.printSciDouble((*mainLOG), rp1,               precision, 20, left); // 
        printer.printSciDouble((*mainLOG), rp2,               precision, 20, left); //
        overallTime = 0;                                                            //
        for (Int i = 0; i <= nStep; ++i) {                                          //
            overallTime += timeSteps[i];                                            //
        }                                                                           //
        averageTime = overallTime/(nStep + 1);                                      //
        printer.printSciDouble((*mainLOG), overallTime,  precision, 20, left);      //
        printer.printSciDouble((*mainLOG), averageTime,  precision, 20, left);      //
        *mainLOG << std::endl;                                                      //
        printer.printString((*mainLOG), "", 190, left, '-');                        //
        *mainLOG << std::endl;                                                      //
    }                                                                               //
    /********************************************************************************/
}

void FTTFelectronicStates::calculateEnergyLevels() {
    Double overallTime;
    /**********************main LOG section**********************/
    if (printMainLogOn) {                                       //
        *mainLOG << "calculating energy levels: " << std::endl; //
        overallTime = mainTimer.getElapsedTimeInMilliSec();     //
    }                                                           //
    /************************************************************/
    for (Int n = 1; n <= nMax; n++) {
        for (Int l = 0; l < n; ++l) {
            if (!calculatedLevel[n][l]) {
                /**********************main LOG section**********************/
                if (printMainLogOn) {                                       //
                    printer.printString((*mainLOG), "n",         10, left); //
                    printer.printString((*mainLOG), "l",         10, left); //
                    printer.printString((*mainLOG), "nStep",     10, left); //
                    printer.printString((*mainLOG), "Ynl_left",  20, left); //
                    printer.printString((*mainLOG), "Ynl_right", 20, left); //
                    printer.printString((*mainLOG), "Ynl",       20, left); //
                    printer.printString((*mainLOG), "Scalc",     20, left); //
                    printer.printString((*mainLOG), "rp1",       20, left); //
                    printer.printString((*mainLOG), "rp2",       20, left); //
                    printer.printString((*mainLOG), "error",     20, left); //
                    printer.printString((*mainLOG), "time[ms]",  20, left); // 
                    (*mainLOG) << std::endl;                                //
                    printer.printString((*mainLOG), "",  190, left, '-');   //
                    (*mainLOG) << std::endl;                                //
                }                                                           //
                /************************************************************/
                calculateEnergyLevel(n, l);
            }
        }
    }
    /**********************main LOG section********************************/
    if (printMainLogOn) {                                                 //
        *mainLOG << "finished calculating energy levels, overall time: "; //
        overallTime = mainTimer.getElapsedTimeInMilliSec() - overallTime; //
        *mainLOG << overallTime << std::endl;                             //
    }                                                                     //
    /**********************************************************************/
}

void FTTFelectronicStates::calculateNlow(Double BE) {
    Double Nhalf = 0;
    Double Nn;
    Double Nnl;

    for (Int n = 1; n <= nMax; ++n) {
        Nn = 0;
        for (Int l = 0; l < n; ++l) {
            if (energyLevel[n][l] < BE) {
                Nnl = (2 * l + 1) / (1 + exp(energyLevel[n][l] - phiShift));
                Nn += Nnl;
            }
        }
        Nhalf += Nn;
    }
    NlowVal = Nhalf * 2;
    calculated_Nlow = true;
    /*************************************main LOG section********************************/
    if (printMainLogOn) {                                                                //
        *mainLOG << "Calculating exact number of electronic states below bound energy: " //
                 << std::endl;                                                           //
        *mainLOG << "Nlow = " << NlowVal << std::endl;                                   //
    }                                                                                    //
    /*************************************************************************************/
}

void FTTFelectronicStates::calculateNTFlow(Double BE) {
    DoubleVec y0(0.0, RHSFTTFstatesLow::dim);
    Double result;
    rhsFTTFstatesLow.setBE(BE);

    y0[0] = phi(1);
    y0[1] = y0[0];
    y0[2] = 0;

    FTTFstatesLowSolver.SetTolerance(1e-7, 0.0);
    FTTFstatesLowSolver.SetStep(1e-15);

    FTTFstatesLowSolver.Integrate(rhsFTTFstatesLow, y0, 1.0, 0.0);
    result = y0[2];

    NTFlowVal = result * 3.0 * T()*V()*sqrt(2.0 * T()) / M_PI / M_PI;
    calculated_NTFlow = true;
    /*************************main LOG section*******************************/
    if (printMainLogOn) {                                                   //
        *mainLOG << "Thomas-Fermi states below bound energy:" << std::endl; //
        *mainLOG << "NTFlow = " << NTFlowVal << std::endl;                  //
    }                                                                       //
    /************************************************************************/
}

void FTTFelectronicStates::calculateNTFfull() {
    DoubleVec y0(0.0, RHSFTTFstatesFull::dim);
    Double result;

    y0[0] = phi(1);
    y0[1] = y0[0];
    y0[2] = 0;

    FTTFstatesFullSolver.SetTolerance(1e-7, 0.0);
    FTTFstatesFullSolver.SetStep(1e-15);

    FTTFstatesFullSolver.Integrate(rhsFTTFstatesFull, y0, 1.0, 0.0);
    result = y0[2];

    NTFfullVal = result * 3.0 * T()*V()*sqrt(2.0 * T()) / M_PI / M_PI;
    calculated_NTFfull = true;
    /*************************main LOG section*****************/
    if (printMainLogOn) {                                     //
        *mainLOG << "Full Thomas-Fermi states:" << std::endl; //
        *mainLOG << "NTFfull = " << NTFfullVal << std::endl;  //
    }                                                         //
    /**********************************************************/
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

    y0[0] = phi(1);
    y0[1] = y0[0];
    y0[2] = 0;

    boundEnergySolver.SetTolerance(1e-8, 0.0);
    boundEnergySolver.SetStep(1e-15);
    boundEnergySolver.Integrate(rhsBoundEnergy, y0, 1.0, 0.0);
    PS = y0[2];

    PS = PS*1.5*T()*V()*sqrt(2.0 * T()) / M_PI / M_PI;

    return PS;
}

Double FTTFelectronicStates::Nlow() {
    calculateEnergyLevels();
    if (!calculated_BE) calculateBoundEnergy();
    if (!calculated_Nlow) calculateNlow(boundEnergy);
    return NlowVal;
}

Double FTTFelectronicStates::NTFlow() {
    calculateEnergyLevels();
    if (!calculated_BE) calculateBoundEnergy();
    if (!calculated_NTFlow) calculateNTFlow(boundEnergy);
    return NTFlowVal;
}

Double FTTFelectronicStates::DNlow() {
    return Nlow() - NTFlow();
}

Double FTTFelectronicStates::N() {
    calculateEnergyLevels();
    if (!calculated_Nlow) calculateNlow(boundEnergy);
    if (!calculated_NTFlow) calculateNTFlow(boundEnergy);
    if (!calculated_NTFfull) calculateNTFfull();
    return NlowVal + NTFfullVal - NTFlowVal;
}

Double FTTFelectronicStates::Nlow(Double BE) {
    calculateEnergyLevels();
    if (!calculated_Nlow) calculateNlow(BE);
    return NlowVal;
}

Double FTTFelectronicStates::NTFlow(Double BE) {
    if (!calculated_NTFlow) calculateNTFlow(BE);
    return NTFlowVal;
}

Double FTTFelectronicStates::DNlow(Double BE) {
    return Nlow(BE) - NTFlow(BE);
}

Double FTTFelectronicStates::N(Double BE) {
    calculateEnergyLevels();
    if (!calculated_Nlow) calculateNlow(BE);
    if (!calculated_NTFlow) calculateNTFlow(BE);
    if (!calculated_NTFfull) calculateNTFfull();
    return NlowVal + NTFfullVal - NTFlowVal;
}