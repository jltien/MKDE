/******************************************************************************
* file_io.h                                                                   *
*                                                                             *
* Functions to check machine properties, swap byte orders, and other file     *
* I/O relation operations.                                                    *
*                                                                             *
* Author: Jeff A. Tracey                                                      *
* License: GPL 2 or newer                                                     *
******************************************************************************/

#ifndef FILE_IO_H
#define FILE_IO_H
#if __cplusplus > 199711L
#define register      // Deprecated in C++11.
#endif  // #if __cplusplus > 199711L
// header file includes

#include "basicDefinitions.h"

#ifdef _ANSI_STD_H       // using ANSI/ISO standard header file names
#include <string>
#include <iostream>      // not required by most systems
#include <fstream>
#include <iomanip>
#include <algorithm>
#else
#include <string.h>
#include <iostream.h>      // not required by most systems
#include <fstream.h>
#include <iomanip.h>
#endif

#include <cstdio>
#include <cstdlib>
//#include <malloc.h>      // for malloc_stats()
#include <sys/stat.h>    // to get file sizes using fstat

// EXPECTED SIZES OF DATA TYPES
const int LONG_BYTES = 4;          // the number of bytes in a long interger
const int FLOAT_BYTES = 4;         // the number of bytes in a floating point number
const int DOUBLE_BYTES = 8;        // the number of bytes in a double floating point number

/*******************************************************************************
* FUNCTIONS TO TEST AND SWAP BYTE ORDER IN BINARY DATA                         *
* Author: Jeff A. Tracey                                                       *
*******************************************************************************/

// machine check functions
bool numDataBytesOK (void);

#define SwapBytes(x) ByteSwap((unsigned char *) &x, sizeof(x))

inline void ByteSwap(unsigned char * b, int n) {
    register int i = 0;
    register int j = n - 1;
    while (i < j) {
        std::swap(b[i], b[j]);
        i++, j--;
    }
}



/******************************************************************************
* FUNCTIONS FOR TEXT FILES                                                    *
* Author: Jeff A. Tracey                                                      *
******************************************************************************/

// eats the rest of a line in an ASCII text file so input can continue on next line

inline void eatline(std::ifstream & file_name) {
	while (file_name.get() != '\n') continue;
}


// returns true if tmpString == T, false if tmpString == F, error otherwise
/*
inline bool boolStringIsTrue(char * tmpString) {
    if (strcmp(tmpString, "T") == 0) {
        return true;
    }
    else if (strcmp(tmpString, "F") == 0) {
        return false;
    }
    else {
        std::cerr << "Error in checkBoolString():" << std::endl;
        std::cerr << "\tArgument is not T or F." << std::endl;
        exit(1);
    }
}
 */

// checks to see if string is T or F; otherwise, returns false
/*
inline bool isProperBoolString(char * tmpString) {
    if ((strcmp(tmpString, "T") == 0)||(strcmp(tmpString, "F") == 0)) {
        return true;
    }
    else {
        return false;
    }
}
 */

/******************************************************************************
* FUNCTION PROTOTYPES                                                         *
* Author: Jeff A. Tracey                                                      *
******************************************************************************/




// load file data into memory
long getFileSize(char * file_name);
char * loadBinaryFile(char * f_name);
void testLBF(char * test_file);
bool isMachineBigEndian(void);

#endif

