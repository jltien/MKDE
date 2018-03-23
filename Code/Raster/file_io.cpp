/******************************************************************************
* file_io.cpp                                                                 *
*                                                                             *
* Functions to check machine properties, swap byte orders, and other file     *
* I/O relation operations.                                                    *
*                                                                             *
* Author: Jeff A. Tracey                                                      *
* License: GPL 2 or newer                                                     *
******************************************************************************/

using namespace std;
#include "file_io.h"

/******************************************************************************
*-----------------------------------------------------------------------------*
*| FUNCTIONS FOR MACHINE-SPECIFIC DATA HANDLING                              |*
*| Author: Jeff A. Tracey                                                    |*
*-----------------------------------------------------------------------------*
******************************************************************************/

/******************************************************************************
* Check characteristics of machine                                            *
******************************************************************************/

/*-----------------------------------------------------------------------------
* Test byte order of machine                                                  *
-----------------------------------------------------------------------------*/
bool isMachineBigEndian(void) {
	union {
		long num;
		unsigned char uc[sizeof(long)];
	} u;
	u.num = 1;
	if (u.uc[0] == 1) {
		return false;
	}
	else {
		return true;
	}
}

/*-----------------------------------------------------------------------------
* Check to see if long, double, and float data types are of the expected size *
-----------------------------------------------------------------------------*/
bool numDataBytesOK (void) {
	if ((sizeof(long)==4)&&(sizeof(double)==8)&&(sizeof(float)==4)&&(sizeof(char)==1)) {
		return true;
	}
	else {
		return false;
	}
}






