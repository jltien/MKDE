/************************************************************************
* basicFunctions.h                                                      *
*                                                                       *
* Created by: Jeff Tracey                                               *
* Date Created: 8/22/2000                                               *
*                                                                       *
* Header file containing constants, data types, and function prototypes *
* used for simulating animal movement.                                  *
*                                                                       *
************************************************************************/

/**** PREPROCESSOR DIRECTIVES ****/

#ifndef BASIC_FUNC_H
#define BASIC_FUNC_H

#include "/home/jeff/projects/myCPPlib/basicDefinitions.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <ctime>

// number of days per month for normal and leap years
const int dayPerMonthNormYr[] = {31,28,31,30,31,30,31,31,30,31,30,31};
const int dayPerMonthLeapYr[] = {31,29,31,30,31,30,31,31,30,31,30,31};

/*-----------------------------------------------------------------------------
* functions that handle floating point problems                               *
-----------------------------------------------------------------------------*/
inline double protectedLog(double x) {
    if (!isnan(x)) {
        if (x >= 0.0) {
            double y = log(x);
            int finCheck = isinf(y);
            if (finCheck == 1) return DBL_MAX; // return the largest possible double
            else if (finCheck == -1) return DBL_MIN; // return the smallest possible double
            else return y;
        }
        else {
            std::cerr << "Error in protectedLog(): x is negative." << std::endl;
            exit(1);
        }
    }
    else {
        std::cerr << "Error in protectedLog(): x is NaN." << std::endl;
        exit(1);
    }
}

inline double protectedExp(double x) {
    if (!isnan(x)) {
        double y = exp(x);
        int finCheck = isinf(y);
        if (finCheck == 1) return DBL_MAX; // return the largest possible double
        else if (finCheck == -1) return DBL_MIN; // return the smallest possible double
        else return y;
    }
    else {
        std::cerr << "Error in protectedExp(): x is NaN." << std::endl;
        exit(1);
    }
}

inline double protectedDivision(double n, double d) {
    if ((!isnan(n))&&(!isnan(d))) {
        if (d == 0.0) {
            if (n == 0.0) {
                std::cerr << "Error in protectedDivision(): n and d are 0." << std::endl;
                std::cerr << "This will result in NaN." << std::endl;
                exit(1);
            }
            else {
                return DBL_MAX; // return max double instead of +inf
            }
        }
        else {
            return n/d;
        }
    }
    else {
        std::cerr << "Error in protectedDivision(): n and/or d is NaN." << std::endl;
        exit(1);
    }
}

// protected sqrt

/*-----------------------------------------------------------------------------
* ADDITIONAL FUNCTIONS THAT NEED TO BE CREATED:                               *
*                                                                             *
*                                                                             *
*                                                                             *
-----------------------------------------------------------------------------*/

inline bool floatValuesWithinTolerance(float u, float v) {
    if (fabs(double(v - u)) <= DBL_TOLERANCE) {
        return true;
    }
    else {
        return false;
    }
}

inline bool doubleValuesWithinTolerance(double u, double v) {
    if (fabs(v - u) <= DBL_TOLERANCE) {
        return true;
    }
    else {
        return false;
    }
}

/*-----------------------------------------------------------------------------
* general functions                                                           *
-----------------------------------------------------------------------------*/
// sum a vector of double
// mean of a vector of double
// mean of a vector of double

/*-----------------------------------------------------------------------------
* angle functions                                                             *
-----------------------------------------------------------------------------*/

// take any angle and return it on (-pi, pi]
inline double angleModulus(double ang) {
    double d = fmod(ang, 2.0*M_PI);
	if(d > M_PI) d -= 2.0*M_PI;
    if(d < -1.0*M_PI) d += 2.0*M_PI;
	return d;
}

inline double vectorToAngle(double delta_x, double delta_y) {
	if ((delta_x != 0.0)||(delta_y != 0.0)) return atan2(delta_y, delta_x);
	else return NAN;
}

inline double vectorLength(double delta_x, double delta_y) {
    double res = hypot(delta_x, delta_y);
    return res;
}

inline double angleSum(double angle1, double angle2) {
	double the_sum = angleModulus(angle1 + angle2);
	return the_sum;
}

// calculate mean vector from angles (use vector class?)
// weighted mean vector
// convert a vector to an angle and length
// convert an angle and length to a vector

/*-----------------------------------------------------------------------------
* move length functions                                                       *
-----------------------------------------------------------------------------*/
double meanMoveLength(double amplitude, double acrophase, double mean_lvl, double angle);


/*-----------------------------------------------------------------------------
* switch functions                                                            *
-----------------------------------------------------------------------------*/
inline double logisticFunc(double x, double c3, double c4) {
	double out_val;
	double lin_eq;
	double lin_to_e;
	
	lin_eq = -c3*(c4 - x);
	lin_to_e = protectedExp(lin_eq);
	out_val = (1/(1 + lin_to_e));
	
	return out_val;
}

inline double logisticFunc(double x, double c1, double c2, double c3, double c4) {
	double out_val;
	double lin_eq;
	double lin_to_e;
	
	lin_eq = -c3*(c4 - x);
	lin_to_e = protectedExp(lin_eq);
	out_val = (c1 - c2)*(1/(1 + lin_to_e)) + c2;
	
	return out_val;
}

inline double exponentialFunc(double x, double c1, double c2, double c3) {
	double exp_arg;
	double out_val;
	
	exp_arg = -x*c2;
	out_val = (c1 - c3)*protectedExp(exp_arg) + c3;
	
	return out_val;
}

inline double saturationFunction(double x, double halfSat, double scale) {
	// what is halfSat and x are 0?
    double res = scale*x/(halfSat + x);
    return res;
}

// contraint between minV and maxV
inline double constrainBetweenValues(double x, double minV, double maxV) {
    if (minV > maxV) {
        double tmpV = minV;
        minV = maxV;
        maxV = tmpV;
    }
    if (x > maxV) return maxV;
    else if (x < minV) return minV;
    else return x;
}

inline double constrainLowerBound(double x, double minV) {
    if (x < minV) return minV;
    else return x;
}

inline double constrainUpperBound(double x, double maxV) {
    if (x > maxV) return maxV;
    else return x;
}

inline int constrainBetweenValues(int x, int minV, int maxV) {
    if (minV > maxV) {
        int tmpV = minV;
        minV = maxV;
        maxV = tmpV;
    }
    if (x > maxV) return maxV;
    else if (x < minV) return minV;
    else return x;
}

inline int constrainLowerBound(int x, int minV) {
    if (x < minV) return minV;
    else return x;
}

inline int constrainUpperBound(int x, int maxV) {
    if (x > maxV) return maxV;
    else return x;
}

/*-----------------------------------------------------------------------------
* ANN node functions                                                          *
-----------------------------------------------------------------------------*/
bool thresholdGate(double x);                              //

inline double binaryGate(double x) {
	double out_val;
	if(x <= 0.0) out_val = 0.0;
	else out_val = 1.0;
	return out_val;
}

inline double peicewiseGate(double x) {
	double out_val;
	if(x < 0.0) out_val = 0.0;
	else if (x > 1.0) out_val = 1.0;
	else out_val = x;
	return out_val;
}

inline double logisticGate(double x) {
    double out_val;
	out_val = 1.0/(1.0 + protectedExp(-1*x));
	return out_val;
}

double logisticGateInv(double y);                          //
double logisticGateDeriv(double x);                        //

/*-----------------------------------------------------------------------------
* time related functions                                                      *
-----------------------------------------------------------------------------*/
// you know, it might be splendid to have a time class

// check to see if a year is a leap year, then return number of days in the year
inline int daysInYear(int yr) {
    if((yr%4 == 0)&&((yr%100 != 0)||(yr%400 == 0))) return 366;
    else return 365;
}

/******************************************************************************
* Calculate the difference betweenfrom time to to time in decimal minutes
******************************************************************************/
inline double timeDiffInMinutes(int fmYr, int fmDay, int fmHr, int fmMin, int fmSec, 
                                int toYr, int toDay, int toHr, int toMin, int toSec) {
    //
    double fmMinutes = double(fmDay)*1440.0 + double(fmHr)*60.0 + double(fmMin) + double(fmSec)/60.0;
    double toMinutes = double(toDay)*1440.0 + double(toHr)*60.0 + double(toMin) + double(toSec)/60.0;
    double res = 0.0;
    int i = 0, daySum = 0;
    if (fmYr == toYr) { // can ignore year
        // subtract fmMinutes from toMinutes; sign depends on args
        res = toMinutes - fmMinutes;
    }
    else if (fmYr < toYr) {
        // add up number of days from begining of fmYr to end of toYr-1
        for (i = fmYr; i < toYr; i++) daySum += daysInYear(i);
        // calculate minutes; will be a positive number (if args correct)
        res = toMinutes - fmMinutes + double(daySum)*1440.0;
    }
    else { // fmYr > toYr
        // add up number of days from begining of toYr to end of fmYr-1
        for (i = toYr; i < fmYr; i++) daySum += daysInYear(i);
        // calculate minutes; will be a negative number (if args correct)
        res = toMinutes - fmMinutes - double(daySum)*1440.0;
    }
    return res;
}

// FUNCTION TO CALCULATE JULIAN DAY GIVEN YEAR, MONTH, DAY_OF_MONTH

// PUT FUNCTIONS TO CALCULATE DAY LENGTH FROM OLD AGENT CLASS(AND LUNAR CYCLE?) HERE...

#endif
