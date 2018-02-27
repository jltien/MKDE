/*****************************************************
 * basicDefinitions.h                                *
 * Global C preprecessor definitions and global      *
 * constants.                                        *
 *                                                   *
 * Author: Jeff A. Tracey                            *
 * Created: 14 December 2002                         *
 * Copyright Jeff A. Tracey 2002.  All rights        *
 * reserved                                          *
 ****************************************************/

#ifndef BASIC_DEFS_H
#define BASIC_DEFS_H

// indicates whether or not I am using ANSI/ISO standard
// header file names
#define _ANSI_STD_H

// set to true to display debugging messages
const bool DO_IBMM_DEBUG = false;

/*******************************************************************************
* Mathematical Constants                                                       *
*******************************************************************************/

// Mathematical constants (note theses are also defined in cmath
const double PI = 3.14159265358979;
const double EE = 2.718281828;

/*******************************************************************************
* Tolerance Values (how far apart values can be before considered not equal)   *
*******************************************************************************/
const float FLT_TOLERANCE = 0.00001;
const double DBL_TOLERANCE = 0.00001;

// General
const char null_file[] = "no_file";

/*******************************************************************************
* CONSTANTS FOR ESRI BINARY RASTER FILES                                       *
* Author: Jeff A. Tracey                                                       *
*******************************************************************************/

const float FLT_NO_DATA_VALUE = -9999.0;
const int INT_NO_DATA_VALUE = -9999;

// STRINGS FOR READING AND WRITING HEADER FILES
const char row_s1[] = "nrows";
const char row_s2[] = "NROWS";
const char col_s1[] = "ncols";
const char col_s2[] = "NCOLS";
const char xll_s1A[] = "xllcorner";
const char xll_s1B[] = "XLLCORNER";
const char xll_s2A[] = "xllcenter";
const char xll_s2B[] = "XLLCENTER";
const char yll_s1A[] = "yllcorner";
const char yll_s1B[] = "YLLCORNER";
const char yll_s2A[] = "yllcenter";
const char yll_s2B[] = "YLLCENTER";
const char cel_s1[] = "cellsize";
const char cel_s2[] = "CELLSIZE";
const char nan_s1[] = "NODATA_value";
const char nan_s2[] = "nodata_value";
const char nan_s3[] = "NODATA_VALUE";
const char byt_s1[] = "byteorder";
const char byt_s2[] = "BYTEORDER";
const char byt_s1A[] = "MSBFIRST";
const char byt_s1B[] = "msbfirst";
const char byt_s2A[] = "LSBFIRST";
const char byt_s2B[] = "lsbfirst";

/*******************************************************************************
* CONSTANTS FOR ESRI BINARY RASTER FILES                                       *
* Author: Jeff A. Tracey                                                       *
*******************************************************************************/
const char ESRI_ASCII[] = "esri_ascii";
const char ESRI_BINARY[] = "esri_binary";



#endif
