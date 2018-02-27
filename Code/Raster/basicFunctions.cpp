/*******************************************************************************
* basicFunctions.cpp                                                           *
*                                                                              *
* Created by: Jeff Tracey                                                      *
* Date Created: 8/22/2000                                                      *
*                                                                              *
* Implementation file containing function definitions for functions            *
* used to simulate animal movement.                                            *
*                                                                              *
*******************************************************************************/

/**** PREPROCESSOR DIRECTIVES ****/
#include "/home/jeff/projects/myCPPlib/basicFunctions.h"

using namespace std;


/*******************************************************************************
* Angles                                                                       *
*******************************************************************************/



/*******************************************************************************
* Switch functions                                                             *
*******************************************************************************/

/*------------------------------------------------------------------------------
* Logistic equation:
* x is the independent variable
* c1 is value when x = -inf
* c2 is value when x = +inf
* c3 is slope at inflection point
* c4 is location of inflection point
* VERIFIED *********************************************************************
------------------------------------------------------------------------------*/

bool thresholdGate(double x) {
    if (x >= 0.0) {
        return true;
    }
    else {
        return false;
    }
}

double logisticGateInv(double y) {
    if((y < 1.00000000000)&&(y > 0.00000000000)) {
        double a, out_val;
        a = (1.0 - y)/y;
        out_val = -1.0*log(a);
        return out_val;
    }
    else {
        cerr << "Error in logisticGateInv: y is not < 1.0 and > 0.0" << endl;
        exit(1);
    }
}

double logisticGateDeriv(double x) {
    double numer, denom_rt, out_val;
    numer = exp(-1*x);
    denom_rt = (1 + numer);
    out_val = numer/(denom_rt*denom_rt);
    return out_val;
}

/*******************************************************************************
* Length-direction cross correlation                                           *
*******************************************************************************/

/*------------------------------------------------------------------------------
* Calculate mean length as a function of an angle (cross-correlation)          *
* Throws an error if amplitude is greater than mean level so check first       *
* Does not catch parameters that are out of range                              *
* VERIFIED *********************************************************************
------------------------------------------------------------------------------*/
double meanMoveLength(double amplitude, double acrophase, double mean_lvl, 
	double angle) {
	double mv_length;
	double diff;
	
	if (amplitude > mean_lvl) {
		// amplitude cannot be greater than mean level as this would allow negative mean
		// move lengths
		cerr << "Error in meanMoveLength(): Amplitude is greater than mean level." << endl;
		exit (1);
	}
	diff = angleSum(angle, -acrophase);
	mv_length = mean_lvl + amplitude*cos(diff);
	return mv_length;
}


