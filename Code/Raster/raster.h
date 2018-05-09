/*******************************************************************************
*  raster.h                                                                    *
*                                                                              *
*  Data types and functions to store and query data from binary raster files   *
*  exported by ArcView spatial analyst.                                        *
*                                                                              *
* NOTE:                                                                        *
* Header file: *.hdr                                                           *
* Data file: *.flt                                                             *
*                                                                              *
*  Copyright Jeff Tracey 2001, All Rights Reserved.                            *
*  Created Nov 4, 2001                                                         *
*******************************************************************************/

#ifndef RASTER_H
#define RASTER_H

/*------------------------------------------------------------------------------
- PREPROCESSOR DIRECTIVES
------------------------------------------------------------------------------*/

#include "basicDefinitions.h"
#include "basicFunctions.h"
#include "file_io.h"    // byte swapping functions

#include <iostream>   // not required by most systems
#include <fstream>
#include <iomanip>
//#include <sys/stat.h>     // to get file sizes using fstat
//#include <fcntl.h>        // for open()
// c headers
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

const int NULL_RASTER_ID = -1;

/*-----------------------------------------------------------------------------
* TYPES OF RASTERS TO IMPLEMENT:
*     1.  FLOATING POINT (FOR CONTINUOUS)
* TYPES OF FILES TO READ FROM/WRITE TO:
*     1.  ESRI
*         A.  BINARY RASTER
*         B.  ASCII RASTER
* CREATE EMPTY RASTERS OF A GIVEN SIZE TO STORE DATA
-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
* NOTE ABOUT COORDINATES AND RASTER INDICES:                                  *
* The column indexes and the x-coordinates both increase from left to right   *
* The row indices increase from top to bottom, but the y-coordinates increase *
* as you go up.                                                               *
* In ESRI binary rasters, the header file specifies the lower-left corner of  *
* the raster, but the data in the flt files is as described above.            *
-----------------------------------------------------------------------------*/

/*============================================================================*
* A CLASS FOR BOUNDING BOX GEOMETRY                                           *
*============================================================================*/


class box {
	private:
        double x_min;
        double y_min;
        double x_max;
        double y_max;
		void correct();                                       //
		// 16 bytes of data
	public:
		box();                                                // default constructor
		box(double xMin, double yMin, double xMax, double yMax); // constructor
		box(const box & b);                                   // copy constructor
		void reset(box & b);                                  //
		void reset(double xMin, double yMin, double xMax, double yMax);  //
		void reset_xMin(double xx);                           //
		void reset_xMax(double xx);                           //
		void reset_yMin(double yy);                           //
		void reset_yMax(double yy);                           //
        bool pointIsInBox(double x, double y);                //
		double get_xMin();                                    //
		double get_yMin();                                    //
		double get_xMax();                                    //
		double get_yMax();                                    //
		double area();                                        //
		bool equals(box & b);                                 //
		void print();                                         //
};


/*-----------------------------------------------------------------------------
* Member Functions for box class
-----------------------------------------------------------------------------*/
inline box::box() {
    x_min = 0.0;
    y_min = 0.0;
    x_max = 0.0;
    y_max = 0.0;
}

inline box::box(double xMin, double yMin, double xMax, double yMax) {
    x_min = xMin;
    y_min = yMin;
    x_max = xMax;
    y_max = yMax;
	correct();
}

inline box::box(const box & b) {
    x_min = b.x_min;
    y_min = b.y_min;
    x_max = b.x_max;
    y_max = b.y_max;
	correct();
}

inline void box::reset(box & b) {
    x_min = b.x_min;
    y_min = b.y_min;
    x_max = b.x_max;
    y_max = b.y_max;
	correct();
}

inline void box::reset(double xMin, double yMin, double xMax, double yMax) {
	x_min = xMin;
    y_min = yMin;
    x_max = xMax;
    y_max = yMax;
	correct();
}

inline void box::reset_xMin(double xx) {
	x_min = xx;
	correct();
}

inline void box::reset_xMax(double xx) {
	x_max = xx;
	correct();
}

inline void box::reset_yMin(double yy) {
	y_min = yy;
	correct();
}

inline void box::reset_yMax(double yy) {
	y_max = yy;
	correct();
}

// fix the bounding box if it is incorrect.  Many of the functions
// for this class rely on xMin < xMax and yMin < yMax.  Calling this
// function after any changes will ensure this condition is true
inline void box::correct() {
	if (x_min > x_max) {
		double tmpx = x_max;
		x_max = x_min;
		x_min = tmpx;
	}
	if (y_min > y_max) {
		double tmpy = y_max;
		y_max = y_min;
		y_min = tmpy;
	}
}


inline bool box::pointIsInBox(double x, double y) {
	if((x < x_min)||(x > x_max)||(y < y_min)||(y > y_max)) {
		return false;
	}
	else {
		return true;
	}
}

inline double box::get_xMin() {
	return x_min;
}

inline double box::get_xMax() {
	return x_max;
}

inline double box::get_yMin() {
	return y_min;
}

inline double box::get_yMax() {
	return y_max;
}

inline double box::area() {
	double a = (x_max - x_min)*(y_max - y_min);
	return a;
}

inline bool box::equals(box & b) {
    if ((x_min == b.x_min) && (x_max == b.x_max) && (y_min == b.y_min) && (y_max = b.y_max)) {
        return true; // really should be a comparison +/- some tolerance
    }
    else {
        return false;
    }
}

inline void box::print() {
	std::cout << "BOX: LL: (" << x_min << ", " << y_min << ")";
	std::cout << " UR: (" << x_max << ", " << y_max << ")" << std::endl;
}


/*============================================================================*
* A CLASS FOR FLOATING-POINT RASTERS                                          *
*============================================================================*/
class gridFloat {
	private:
        // properties of the raster
	    long num_rows;                                               // number of rows
	    long num_cols;                                               // number of cols
	    double cell_size;                                            // cell width
	    float no_dat;                                                // no data value
        //      there may be others (e.g., row_major or col_major)
        // grid data
	    box grid_bbox;                                               // bounding box for grid
	    float ** data;                                               // grid data 2D dynamic Array
        int rasterID;                                                // an ID number fo rthis raster layer
	public:
	    gridFloat();                                                 // default constructor
        gridFloat(int rID, long nr, long nc, double xll, double yll, double cellsz);  // contruct and array of 0.0s
	    gridFloat(int rID, char * filenm, char * filetype);           // filetype = (esri_ascii, esri_binary)
        gridFloat(const gridFloat & inGrid);                         // copy constructor
        ~gridFloat();                                                // destructor
	
        bool isValidRowIndex(long i);                                //
        bool isValidColumnIndex(long j);                             //

        // access to grid information
        int getRasterID() { return rasterID; }                       // 
        box get_minBoundingBox() {return grid_bbox;}                 // 
        double getCellSize() { return cell_size; }                   // 
        long getNumberOfRows() {return num_rows; }                   // 
        long getNumberOfColumns() { return num_cols; }               // 
        float getNoDataValue() {return no_dat; }                     //

        // get indices into raster given location coordinates
        long getGridRow(double ycoord);                              // get row index
        long getGridCol(double xcoord);                              // get column index
        double getCellCenterXcoord(long colIndex);                   // get cell center xcoord given index
        double getCellCenterYcoord(long rowIndex);                   // get cell center xcoord given index

        // set raster values
        void setGridProperties(int rID, long nr, long nc, double xll, double yll, double cellsz); //
        void readDataFromFile(int rID, char * filenm, char * filetype);       //
        void readDataFromGrid(gridFloat & fgrd);                     // copy another grid into this one

        // set or change values
        void setGridValue(long rw, long cl, float val);              // set value
        void setGridValue(double xcoord, double ycoord, float val);  // set value

        // modify individual cells
        void addValueToGridCell(long rw, long cl, float val);        // add value to grid cell
        void addValueToGridCell(double xcoord, double ycoord, float val);  // add value to grid cell
        
        // modify entire grid
        void scaleGridCellValue(long rw, long cl, float val);        // multiplify grid cell valule
        void scaleGridCellValue(double xcoord, double ycoord, float val);  // multiplify grid cell valule
        void setAllCellsToZero(bool changeMissing=false);          // set all cells to 0.0
        void setAllCellsToNoData();                                  // set all cells to no data
        void setAllCellsTo(float val);                               // set all cells to given value    
        void scaleFromZeroToOne();                                   //

        // calculations and queries
        bool cellIsMissingData(long rw, long cl);                    //
        bool cellIsMissingData(double x, double y);                  //
        float getGridValue(long rw, long cl);                        // get value
        float getGridValue(double xcoord, double ycoord);            // get value

        // printing functions
        void summary();                                              // displays grid summary to standard output
        void display();                                              // print to standard output
        void printESRIascii(char * filen);                           // print to ESRI ASCII raster
        void printESRIbinary(char * filen);                          // print to ESRI BINARY raster
};

/*------------------------------------------------------------------------------
- INLINE FUNCTIONS
------------------------------------------------------------------------------*/

inline bool gridFloat::isValidRowIndex(long i) {
    if ((i >= 0)&&(i < num_rows)) return true;
    else return false;
}

inline bool gridFloat::isValidColumnIndex(long j) {
    if ((j >= 0)&&(j < num_cols)) return true;
    else return false;
}

// ----------------------------------------------------------------------------

inline long gridFloat::getGridRow(double ycoord) {
    // check to ensure point is in bounds
	if((ycoord >= grid_bbox.get_yMin())&&(ycoord <= grid_bbox.get_yMax())) {
		return long(floor(double(num_rows) - (ycoord - grid_bbox.get_yMin())/cell_size)); // if cell 0,0 is in lower left corner
	}
	else {
		std::cerr << "Error in gridFloat::getGridRow(): query point out of bounds." << std::endl;
        std::cerr << "Bounding box y min = " << grid_bbox.get_yMin() << "; bounding box y max = " << grid_bbox.get_yMax() << std::endl;
        std::cerr << "Arg y coordinate = " << ycoord << std::endl;
        exit(1);
	}
}

inline long gridFloat::getGridCol(double xcoord) {
    // check to ensure point is in bounds
	if ((xcoord >= grid_bbox.get_xMin())&&(xcoord <= grid_bbox.get_xMax())) {
		return (long)(floor((xcoord - grid_bbox.get_xMin())/cell_size));
	}
	else {
		std::cerr << "Error in gridFloat::getGridRow(): query point out of bounds." << std::endl;
        std::cerr << "Bounding box x min = " << grid_bbox.get_xMin() << "; bounding box x max = " << grid_bbox.get_xMax() << std::endl;
        std::cerr << "Arg x coordinate = " << xcoord << std::endl;
        exit(1);
	}
}


inline double gridFloat::getCellCenterXcoord(long colIndex) {
    double res;
    if ((colIndex >= 0)&&(colIndex< num_cols)) {
        res = grid_bbox.get_xMin() + cell_size*double(colIndex) + 0.5*cell_size;
    }
    else {
        // return an out-of-bounds location
        res = grid_bbox.get_xMin() - cell_size;
    }
    return res;
}

inline double gridFloat::getCellCenterYcoord(long rowIndex) {
    double res;
    if ((rowIndex >= 0)&&(rowIndex < num_rows)) {
        res = grid_bbox.get_yMax() - cell_size*double(rowIndex) - 0.5*cell_size;
    }
    else {
        // return an out-of-bounds location
        res = grid_bbox.get_yMin() - cell_size;
    }
    return res;
}

// ----------------------------------------------------------------------------

inline void gridFloat::setGridProperties(int rID, long nr, long nc, double xll, double yll, double cellsz) {
    if ((data == NULL)&&(num_rows == -1)&&(num_cols == -1)) {
        if (rID >= 0) {
            rasterID = rID;
        }
        else {
            rasterID = NULL_RASTER_ID;
        }
        if ((nr > 0)&&(nc > 0)&&(cellsz > 0.0)) {
            // 1.  set member variables
            num_rows = nr;
            num_cols = nc;
            cell_size = cellsz;
            // 2.  create dynamic array
            long i, j;
            data = new float*[num_rows];
            for (i = 0; i < num_rows; i++) {
                data[i] = new float[num_cols];
            }
            for (i = 0; i < num_rows; i++) {
                for (j = 0; j < num_cols; j++) {
                    data[i][j] = FLT_NO_DATA_VALUE;
                }
            }
            // set up bounding box
            double xur, yur;
            xur = xll + cellsz*long(num_cols);
            yur = yll + cellsz*long(num_rows);
            grid_bbox.reset(xll, yll, xur, yur);
        }
        else {
            std::cerr << "Error in setGridProperties(): Numbers of rows and columns, and cell size must be positive." << std::endl;
            exit(1);
        }
        
    }
    // else, it has already been set, so you can't reset it
}

inline void gridFloat::setGridValue(long rw, long cl, float val) {
    if((rw >= 0)&&(rw < num_rows)&&(cl >= 0)&&(cl < num_cols)) {
		data[rw][cl] = val;
	}
    // else, do nothing
}

inline void gridFloat::setGridValue(double xcoord, double ycoord, float val) {
    if ((xcoord >= grid_bbox.get_xMin())&&(xcoord <= grid_bbox.get_xMax())&&(ycoord >= grid_bbox.get_yMin())&&(ycoord <= grid_bbox.get_yMax())) {
        data[getGridRow(ycoord)][getGridCol(xcoord)] = val;
    }
    // else, do nothing
}

inline void gridFloat::addValueToGridCell(long rw, long cl, float val) {
    if((rw >= 0)&&(rw < num_rows)&&(cl >= 0)&&(cl < num_cols)) {
		if (data[rw][cl] != no_dat) {
            data[rw][cl] = data[rw][cl] + val;
        }
        // no data is still no data
	}
    // else, do nothing
}


inline void gridFloat::addValueToGridCell(double xcoord, double ycoord, float val) {
    if ((xcoord >= grid_bbox.get_xMin())&&(xcoord <= grid_bbox.get_xMax())&&(ycoord >= grid_bbox.get_yMin())&&(ycoord <= grid_bbox.get_yMax())) {
        long rw, cl;
        rw = getGridRow(ycoord);
        cl = getGridCol(xcoord);
        if (data[rw][cl] != no_dat) {
            data[rw][cl] = data[rw][cl] + val;
        }
        // no data is still no data
    }
    // else, do nothing
}

inline void gridFloat::scaleGridCellValue(long rw, long cl, float val) {
    if((rw >= 0)&&(rw < num_rows)&&(cl >= 0)&&(cl < num_cols)) {
		if (data[rw][cl] != no_dat) {
            data[rw][cl] = val*data[rw][cl];
        }
        // no data is still no data
	}
    // else, do nothing
}

inline void gridFloat::scaleGridCellValue(double xcoord, double ycoord, float val) {
    if ((xcoord >= grid_bbox.get_xMin())&&(xcoord <= grid_bbox.get_xMax())&&(ycoord >= grid_bbox.get_yMin())&&(ycoord <= grid_bbox.get_yMax())) {
        long rw, cl;
        rw = getGridRow(ycoord);
        cl = getGridCol(xcoord);
        if (data[rw][cl] != no_dat) {
            data[rw][cl] = val*data[rw][cl];
        }
        // no data is still no data
    }
    // else, do nothing
}

inline void gridFloat::setAllCellsToZero(bool changeMissing) {
    if ((data != NULL)&&(num_rows > 0)&&(num_cols > 0)) {
        long i, j;
        for (i = 0; i < num_rows; i++) {
            for (j = 0; j < num_cols; j++) {
                if ((data[i][j] != no_dat)||(changeMissing)) {
                    data[i][j] = 0.0;
                }
            }
        }
    }
    // else, do nothing
}

inline void gridFloat::setAllCellsToNoData() {
    if ((data != NULL)&&(num_rows > 0)&&(num_cols > 0)) {
        long i, j;
        for (i = 0; i < num_rows; i++) {
            for (j = 0; j < num_cols; j++) {
                data[i][j] = no_dat;
            }
        }
    }
    // else, do nothing
}

inline void gridFloat::setAllCellsTo(float val) {
    if ((data != NULL)&&(num_rows > 0)&&(num_cols > 0)) {
        long i, j;
        for (i = 0; i < num_rows; i++) {
            for (j = 0; j < num_cols; j++) {
                data[i][j] = val;
            }
        }
    }
    // else, do nothing
}


// ----------------------------------------------------------------------------

inline bool gridFloat::cellIsMissingData(long rw, long cl) {
    if((rw < 0)||(rw >= num_rows)||(cl < 0)||(cl >= num_cols)) {
		return true;
	}
	else {
		if (floatValuesWithinTolerance(data[rw][cl], no_dat)) {
            return true;
        }
        else {
            return false;
        }
	}
}

inline bool gridFloat::cellIsMissingData(double x, double y) {
    if((x <= grid_bbox.get_xMin())||(x >= grid_bbox.get_xMax())||(y <= grid_bbox.get_yMin())||(y >= grid_bbox.get_yMax())) {
		return true;
	}
	else {
        long rw = getGridRow(y);
        long cl = getGridCol(x);
		if (floatValuesWithinTolerance(data[rw][cl], no_dat)) {
            return true;
        }
        else {
            return false;
        }
	}
}

inline float gridFloat::getGridValue(long rw, long cl) {
	if((rw < 0)||(rw >= num_rows)||(cl < 0)||(cl >= num_cols)) {
        return no_dat;
	}
	else {
		return data[rw][cl];
	}
}


inline float gridFloat::getGridValue(double xcoord, double ycoord) {
    if((xcoord >= grid_bbox.get_xMin())&&(xcoord <= grid_bbox.get_xMax())&&(ycoord >= grid_bbox.get_yMin())&&(ycoord <= grid_bbox.get_yMax())) {
        return getGridValue(getGridRow(ycoord), getGridCol(xcoord));
	}
	else {
		return no_dat;
	}
}


#endif

