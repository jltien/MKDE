/*******************************************************************************
* raster_io.cpp                                                                *
*                                                                              *
* Binray ratser functions for individual-based movement models on              *
* heterogeneous landscapes.                                                    *
* Revised version.                                                             *
*                                                                              *
* Author: Jeff A. Tracey                                                       *
* Copyright (c)Jeff A. Tracey 2001 All Rights Reserved except where otherwise  *
* noted.                                                                       *
* Authorization to copy, use, distribute or modify this code in any form is    *
* by written consent from Jeff A. Tracey only.                                 *
*******************************************************************************/

// ONE PROBLEM IS THAT WE WILL NEED TO INCLUDE THE BYTE SWAPPING FUNCTION FOR FLOAT AND
// A BYTE ORDER, AND SIZE OF FLOAT CHECK

#include "raster.h"


/*******************************************************************************
* BOUNDING BOX                                                                 *
*******************************************************************************/



/*******************************************************************************
* FLOATING-POINT RASTER                                                        *
*******************************************************************************/

/*------------------------------------------------------------------------------
* Default Constructor                                                          *
------------------------------------------------------------------------------*/
gridFloat::gridFloat() {
    rasterID = NULL_RASTER_ID;
    num_rows = -1;
	num_cols = -1;
	cell_size = 0.0;
	no_dat = FLT_NO_DATA_VALUE;
	data = NULL;
}



/*------------------------------------------------------------------------------
* Constructor to make a raster of no data values                               *
------------------------------------------------------------------------------*/
gridFloat::gridFloat(int rID, long nr, long nc, double xll, double yll, double cellsz) {
    xmin = xll;
    ymin = yll;
    xsize = cellsz;
    ysize = cellsz;

    if (rID >= 0) {
            rasterID = rID;
        }
        else {
            rasterID = NULL_RASTER_ID;
    }

    // Note: bounding box set in initializer list
    int i, j;
    // 1. set number of rows
    if (nr > 0) {
        num_rows = nr;
    }
    else {
        std::cerr << "Error in gridFloat() constructor: Number of rows must be positive." << std::endl;
        exit(1);
    }

    // 2. set number of columns
    if (nc > 0) {
        num_cols = nc;
    }
    else {
        std::cerr << "Error in gridFloat() constructor: Number of columns must be positive." << std::endl;
        exit(1);
    }

    // 3. set no data value
    no_dat = FLT_NO_DATA_VALUE;

    // 4. set up 2D dynamic array
    data = new float*[num_rows];
    for (i = 0; i < num_rows; i++) {
        data[i] = new float[num_cols];
    }

    // 5. set array elements to no data value
    for (i = 0; i < num_rows; i++) {
        for (j = 0; j < num_cols; j++) {
            data[i][j] = FLT_NO_DATA_VALUE;
        }
    }

    // 6. set the cell size
    cell_size = cellsz;


    // 7. set the bounding box
    double xur_val, yur_val;
    xur_val = xll + cell_size*long(num_cols);
    yur_val = yll + cell_size*long(num_rows);
    grid_bbox.reset(xll, yll, xur_val, yur_val);

    xmax = xur_val;
    ymax = yur_val;
}

/*------------------------------------------------------------------------------
* Constructor that reads the data from a file                                  *
* ---------------------------------------------------------------------------- *
* NOTES FOR ESRI BINARY AND ASCII RASTER FILES:                                *
* 1.  Row corresponds to y-coordinate, and column corresponds to x-coordinate. *
* 2.  The first row and column position represents the lower-left corner of    *
*     the raster, so calculation of the row index from y-coordinate must       *
*     account for this.                                                        *
------------------------------------------------------------------------------*/
gridFloat::gridFloat(int rID, char * filenm, char * filetype) {
	char row_field[13];
	char col_field[13];
	char xll_field[13];
	char yll_field[13];
	char cell_field[13];
	char nan_field[14];
	char byte_field[15];
	char byte_order[15];
	bool center_flag = false;
    double xll_val, yll_val, xur_val, yur_val; // coordinate for corners of raster
    double half_w;
    bool binFileIsBigEndian;                   // true if data MSBF
    bool machIsBigEndian;                      // true if machine is big endian
	bool byteOrderNotSame;                     // true if the machine and the data are different byte order

    long i, j;                                 // loop indices
	long numf_bytes;        // DON'T NEED?

    if (rID >= 0) {
        rasterID = rID;
    }
    else {
        rasterID = NULL_RASTER_ID;
    }


    std::ifstream fin;


    // ESRI ASCII OR BINARY FLOATING-POINT RASTER =============================
    if ((strcmp(filetype, ESRI_ASCII) == 0)||(strcmp(filetype, ESRI_BINARY) == 0)) {
        // read header data
        // open file ----------------------------------------------------------
    	fin.open(filenm);
    	if (!fin.is_open()) {
    		std::cerr << "Error in gridFloat() constructor: File "<< filenm << " could not be opened." << std::endl;
    		exit(1);
    	}
    	
    	// read in header values ----------------------------------------------
    
        // 1. number of columns
    	fin >> col_field;
    	if ((strcmp(col_field, col_s1) == 0)||(strcmp(col_field, col_s2) == 0)) {
            fin >> num_cols;
       	}
    	else {
    		std::cerr << "Error in gridFloat() constructor [1]: Header file not formated correctly." << std::endl;
    		exit(1);
    	}
    	eatline(fin);
    
        // 2. number of rows
    	fin >> row_field;
    	if ((strcmp(row_field, row_s1) == 0)||(strcmp(row_field, row_s2) == 0)) {
            fin >> num_rows;
       	}
    	else {
            std::cerr << row_field << " != " << row_s1 << " or " << row_s2 << std::endl;
    		std::cerr << "Error in gridFloat() constructor [2]: Header file not formated correctly." << std::endl;
    		exit(1);
    	}
    	eatline(fin);
    
        // 3.  read lower left corner values
    	fin >> xll_field;
        fin >> xll_val;
        eatline(fin);
    	if ((strcmp(xll_field, xll_s1A) == 0)||(strcmp(xll_field, xll_s1B) == 0)) {
    		fin >> yll_field;
    		if ((strcmp(yll_field, yll_s1A) == 0)||(strcmp(yll_field, yll_s1B) == 0)) { 
    			center_flag = false;
                fin >> yll_val;
    		}
    		else {
    			std::cerr << "Error in gridFloat() constructor [3]: Header file not formated correctly." << std::endl;
                exit(1);
    		}
        }
        else if ((strcmp(xll_field, xll_s2A) == 0)||(strcmp(xll_field, xll_s2B) == 0)) {  // center value, not corner
            // NOTE: IN THE ESRI DOCUMENTATION IT IS NOT CLEAR TO ME THAT THIS IS THE CENTER
            // OF THE ENTIRE RASTER, OR THE CENTER OF THE LOWER LEFT CELL.  I ASSUME IT IS
            // THE CENTER OF THE LOWER LEFT CELL UNTIL I FIND OUT OTHERWISE.
            fin >> yll_field;
    		if ((strcmp(yll_field, yll_s2A) == 0)||(strcmp(yll_field, yll_s2B) == 0)) { 
    			center_flag = true;
                fin >> yll_val;
    		}
    		else {
    			std::cerr << "Error in gridFloat() constructor [4]: Header file not formated correctly." << std::endl;
                exit(1);
    		}
        }
        else { // incorrect flag
            std::cerr << "Error in gridFloat() constructor [5]: Header file not formated correctly." << std::endl;
            exit(1);
        }
    	eatline(fin);

        // 4, read in width (= size, = height) of cell
    	fin >> cell_field;
    	if ((strcmp(cell_field, cel_s1) == 0)||(strcmp(cell_field, cel_s2) == 0)) {
            fin >> cell_size;
    	}
    	else {
    		std::cerr << "Error in gridFloat() constructor [6]: Header file not formated correctly." << std::endl;
    		exit(1);
    	}
    	eatline(fin);
    
    	// 5. Read in no data value
        fin >> nan_field;
    	if ((strcmp(nan_field, nan_s1) == 0)||(strcmp(nan_field, nan_s2) == 0)||(strcmp(nan_field, nan_s3) == 0)) {
    		fin >> no_dat;
    	}
    	else {
    		std::cerr << "Error in gridFloat() constructor [7]: Header file not formated correctly." << std::endl;
    		exit(1);
    	}
    	eatline(fin);

        // 6. set up bounding box
        // if raster cell center values used, calculate lower left corner coordinates of ll cell
        if (center_flag) {
    		half_w = cell_size/(2.0);
    		xll_val = xll_val - half_w;
    		yll_val = yll_val - half_w;
    	}
        // set upper right corner values
        xur_val = xll_val + (double(num_cols))*cell_size;
    	yur_val = yll_val + (double(num_rows))*cell_size;
        // make bounding box for grid extent
    	grid_bbox.reset(xll_val, yll_val, xur_val, yur_val);

        // set up 2D dynamic array
        data = new float*[num_rows];
        for (i = 0; i < num_rows; i++) {
            data[i] = new float[num_cols];
        }

        if (strcmp(filetype, ESRI_ASCII) == 0) { // ASCII DATA
            // 7. Read in data
            // read raster data
            for(i = 0; i < num_rows; i++) {
                for(j = 0; j < num_cols; j++) {
                    fin >> data[i][j];
                }
                eatline(fin);
            }
    
    
            // 8. close file(s)
            fin.close();
        }

        else { // BINARY DATA
            // 7. read byte order field
            fin >> byte_field;
            fin >> byte_order;
            if ((strcmp(byte_order, byt_s1A) == 0)||(strcmp(byte_order, byt_s1B) == 0)) {
                binFileIsBigEndian = true;
            }
            else if ((strcmp(byte_order, byt_s2A) == 0)||(strcmp(byte_order, byt_s2B) == 0)) {
                binFileIsBigEndian = false;
            }
        	else {
        		std::cerr << "Error in gridFloat() constructor [8]: ";
        		std::cerr << "Header file not formated correctly." << std::endl;
        		exit(1);
        	}

            // 8. close header file
        	fin.close();

            // 9. set name of binary data file
            int name_len = strlen(filenm);
        	char flt_file[50];                   // char array for index file name
        	strcpy(flt_file, filenm);
        	flt_file[name_len - 3] = 'f';
        	flt_file[name_len - 2] = 'l';
        	flt_file[name_len - 1] = 't';

            // 10. Check byte order of machine and compare to file
            numf_bytes = long(FLOAT_BYTES);           // 4 bytes
            machIsBigEndian = isMachineBigEndian();   // check byte order of machine
            // compare byte order of machine to byte order of file
            if ((machIsBigEndian)&&(binFileIsBigEndian)) {
                byteOrderNotSame = false; // both machine and file data big endian
            }
            else if ((!machIsBigEndian)&&(!binFileIsBigEndian)) {
                byteOrderNotSame = false; // both machine and file data little endian
            }
            else {
                byteOrderNotSame = true; // machine and file have different byte orders, must swap
           }

            // 11. open binary data file
            std::ifstream dfin(flt_file, std::ios_base::binary);
            if (!dfin.is_open()) {
        		std::cerr << "Error in gridFloat() constructor [9]: Data file "<< flt_file << " could not be opened." << std::endl;
        		exit(1);
        	}
            
            // read data into array
            float tmp_float;
            for(i = 0; i < num_rows; i++) {
                for(j = 0; j < num_cols; j++) {
                    dfin.read((char *)&tmp_float, sizeof tmp_float);
                    if (byteOrderNotSame) SwapBytes(tmp_float);
                    data[i][j] = tmp_float;
                }
            }

            // 12. close data file
            dfin.close();
        }
    }
    // ERROR -- UNKNOWN FILE TYPE =============================================
    else {
        std::cerr << "Error in gridFloat() constructor: Unknown file type." << std::endl;
        exit(1);
    }
}

/*------------------------------------------------------------------------------
* Copy constructor                                                             *
------------------------------------------------------------------------------*/
gridFloat::gridFloat(const gridFloat & inGrid) : grid_bbox(inGrid.grid_bbox) {
    rasterID = inGrid.rasterID;
    // Note: bounding box set in initializer list
    int i, j;
    // 1. set number of rows
    num_rows = inGrid.num_rows;

    // 2. set number of columns
    num_cols = inGrid.num_cols;

    // 3. set no data value
    no_dat = inGrid.no_dat;

    if ((data != NULL)&&(num_rows > 0)&&(num_cols > 0)) {
        // 4. set up 2D dynamic array
        data = new float*[num_rows];
        for (i = 0; i < num_rows; i++) {
            data[i] = new float[num_cols];
        }
    
        // 5. set array elements to those in inGrid
        for (i = 0; i < num_rows; i++) {
            for (j = 0; j < num_cols; j++) {
                data[i][j] = inGrid.data[i][j];
            }
        }
    }
    else {
        data = NULL;
    }
    
    // 6. set the cell size
    cell_size = inGrid.cell_size;
}

/*------------------------------------------------------------------------------
* Destructor                                                                   *
------------------------------------------------------------------------------*/
gridFloat::~gridFloat() {
    long i;
    for (i = 0; i < num_rows; i++) {
        delete [] data[i];
        data[i]  = NULL;
    }
    delete [] data;
    data = NULL;
    num_rows = 0;
    num_cols = 0;
    // may want to zero out other stuff
}

void gridFloat::readDataFromFile(int rID, char * filenm, char * filetype) {
    char row_field[13];
	char col_field[13];
	char xll_field[13];
	char yll_field[13];
	char cell_field[13];
	char nan_field[14];
	char byte_field[15];
	char byte_order[15];
    long tmp_num_cols, tmp_num_rows;
    double tmp_cell_size;
    box tmp_bbox;
	bool center_flag = false;
    double xll_val, yll_val, xur_val, yur_val; // coordinate for corners of raster
    double half_w;
    bool binFileIsBigEndian;                   // true if data MSBF
    bool machIsBigEndian;                      // true if machine is big endian
	bool byteOrderNotSame;                     // true if the machine and the data are different byte order

    long i, j;                                 // loop indices
	long numf_bytes;        // DON'T NEED?

    if (rID >= 0) {
        rasterID = rID;
    }
    else {
        rasterID = NULL_RASTER_ID;
    }

    std::ifstream fin;

    // ESRI ASCII OR BINARY FLOATING-POINT RASTER =============================
    if ((strcmp(filetype, ESRI_ASCII) == 0)||(strcmp(filetype, ESRI_BINARY) == 0)) {
        // read header data
        // open file ----------------------------------------------------------
    	fin.open(filenm);
    	if (!fin.is_open()) {
    		std::cerr << "Error in readDataFromFile(): File "<< filenm << " could not be opened." << std::endl;
    		exit(1);
    	}
    	
    	// read in header values ----------------------------------------------
    
        // 1. number of columns
    	fin >> col_field;
    	if ((strcmp(col_field, col_s1) == 0)||(strcmp(col_field, col_s2) == 0)) {
            fin >> tmp_num_cols;
       	}
    	else {
    		std::cerr << "Error in readDataFromFile() [1]: Header file not formated correctly." << std::endl;
    		exit(1);
    	}
    	eatline(fin);
    
        // 2. number of rows
    	fin >> row_field;
    	if ((strcmp(row_field, row_s1) == 0)||(strcmp(row_field, row_s2) == 0)) {
            fin >> tmp_num_rows;
       	}
    	else {
            std::cerr << row_field << " != " << row_s1 << " or " << row_s2 << std::endl;
    		std::cerr << "Error in readDataFromFile() [2]: Header file not formated correctly." << std::endl;
    		exit(1);
    	}
    	eatline(fin);
    
        // 3.  read lower left corner values
    	fin >> xll_field;
        fin >> xll_val;
        eatline(fin);
    	if ((strcmp(xll_field, xll_s1A) == 0)||(strcmp(xll_field, xll_s1B) == 0)) {
    		fin >> yll_field;
    		if ((strcmp(yll_field, yll_s1A) == 0)||(strcmp(yll_field, yll_s1B) == 0)) { 
    			center_flag = false;
                fin >> yll_val;
    		}
    		else {
    			std::cerr << "Error in readDataFromFile() [3]: Header file not formated correctly." << std::endl;
                exit(1);
    		}
        }
        else if ((strcmp(xll_field, xll_s2A) == 0)||(strcmp(xll_field, xll_s2B) == 0)) {  // center value, not corner
            // NOTE: IN THE ESRI DOCUMENTATION IT IS NOT CLEAR TO ME THAT THIS IS THE CENTER
            // OF THE ENTIRE RASTER, OR THE CENTER OF THE LOWER LEFT CELL.  I ASSUME IT IS
            // THE CENTER OF THE LOWER LEFT CELL UNTIL I FIND OUT OTHERWISE.
            fin >> yll_field;
    		if ((strcmp(yll_field, yll_s2A) == 0)||(strcmp(yll_field, yll_s2B) == 0)) { 
    			center_flag = true;
                fin >> yll_val;
    		}
    		else {
    			std::cerr << "Error in readDataFromFile() [4]: Header file not formated correctly." << std::endl;
                exit(1);
    		}
        }
        else { // incorrect flag
            std::cerr << "Error in readDataFromFile() [5]: Header file not formated correctly." << std::endl;
            exit(1);
        }
    	eatline(fin);

        // 4, read in width (= size, = height) of cell
    	fin >> cell_field;
    	if ((strcmp(cell_field, cel_s1) == 0)||(strcmp(cell_field, cel_s2) == 0)) {
            fin >> tmp_cell_size;
    	}
    	else {
    		std::cerr << "Error in readDataFromFile() [6]: Header file not formated correctly." << std::endl;
    		exit(1);
    	}
    	eatline(fin);

        // 5. Calculate bounding box
        // if raster cell center values used, calculate lower left corner coordinates of ll cell
        if (center_flag) {
    		half_w = tmp_cell_size/(2.0);
    		xll_val = xll_val - half_w;
    		yll_val = yll_val - half_w;
    	}
        // set upper right corner values
        xur_val = xll_val + (double(tmp_num_cols))*tmp_cell_size;
    	yur_val = yll_val + (double(tmp_num_rows))*tmp_cell_size;
        // make bounding box for grid extent
    	tmp_bbox.reset(xll_val, yll_val, xur_val, yur_val);

        // 6.  Check to see if grid size has been set and if so whether if matches the size in the file
        bool gridIsSet = false;
        bool gridIsSameSize = false;
        if ((data == NULL)&&(num_rows == -1)&&(num_cols == -1)) {
            gridIsSet = false;
        }
        else {
            gridIsSet = true;
            if ((doubleValuesWithinTolerance(tmp_cell_size, cell_size))&&(tmp_num_rows == num_rows)&&(tmp_num_cols == num_cols)) {
                if (grid_bbox.equals(tmp_bbox)) {
                    gridIsSameSize = true;
                }
                else {
                    gridIsSameSize = false;
                }
            }
            else {
                gridIsSameSize = false;
            }
        }

        // 7. Make any necessary changes to the size of the grid
        if ((gridIsSet)&&(!gridIsSameSize)) {
            // delete old data so can resize grid
            for (i = 0; i < num_rows; i++) {
                delete [] data[i];
                data[i]  = NULL;
            }
            delete [] data;
            gridIsSet = false;
        }

        if (!gridIsSet) {
            // set grid properties
            num_rows = tmp_num_rows;
            num_cols = tmp_num_cols;
            cell_size = tmp_cell_size;
            grid_bbox.reset(tmp_bbox);
            // create dynamic array
            data = new float*[num_rows];
            for (i = 0; i < num_rows; i++) {
                data[i] = new float[num_cols];
            }
            gridIsSet = true;
        }

        // Continue to set data value and read in data

    	// 8. Read in no data value
        fin >> nan_field;
    	if ((strcmp(nan_field, nan_s1) == 0)||(strcmp(nan_field, nan_s2) == 0)||(strcmp(nan_field, nan_s3) == 0)) {
    		fin >> no_dat;
    	}
    	else {
    		std::cerr << "Error in readDataFromFile() [7]: Header file not formated correctly." << std::endl;
    		exit(1);
    	}
    	eatline(fin);
        

        if (strcmp(filetype, ESRI_ASCII) == 0) { // ASCII DATA
            // 9. Read in data
            // read raster data
            for(i = 0; i < num_rows; i++) {
                for(j = 0; j < num_cols; j++) {
                    fin >> data[i][j];
                }
                eatline(fin);
            }
    
    
            // 10. close file(s)
            fin.close();
        }

        else { // BINARY DATA
            // 9. read byte order field
            fin >> byte_field;
            fin >> byte_order;
            if ((strcmp(byte_order, byt_s1A) == 0)||(strcmp(byte_order, byt_s1B) == 0)) {
                binFileIsBigEndian = true;
            }
            else if ((strcmp(byte_order, byt_s2A) == 0)||(strcmp(byte_order, byt_s2B) == 0)) {
                binFileIsBigEndian = false;
            }
        	else {
        		std::cerr << "Error in readDataFromFile() [8]: ";
        		std::cerr << "Header file not formated correctly." << std::endl;
        		exit(1);
        	}

            // 10. close header file
        	fin.close();

            // 11. set name of binary data file
            int name_len = strlen(filenm);
        	char flt_file[50];                   // char array for index file name
        	strcpy(flt_file, filenm);
        	flt_file[name_len - 3] = 'f';
        	flt_file[name_len - 2] = 'l';
        	flt_file[name_len - 1] = 't';

            // 12. Check byte order of machine and compare to file
            numf_bytes = long(FLOAT_BYTES);           // 4 bytes
            machIsBigEndian = isMachineBigEndian();   // check byte order of machine
            // compare byte order of machine to byte order of file
            if ((machIsBigEndian)&&(binFileIsBigEndian)) {
                byteOrderNotSame = false; // both machine and file data big endian
            }
            else if ((!machIsBigEndian)&&(!binFileIsBigEndian)) {
                byteOrderNotSame = false; // both machine and file data little endian
            }
            else {
                byteOrderNotSame = true; // machine and file have different byte orders, must swap
           }

            // 13. open binary data file
            std::ifstream dfin(flt_file, std::ios_base::binary);
            if (!dfin.is_open()) {
        		std::cerr << "Error in readDataFromFile()readDataFromFile [9]: Data file "<< flt_file << " could not be opened." << std::endl;
        		exit(1);
        	}
            
            // read data into array
            float tmp_float;
            for(i = 0; i < num_rows; i++) {
                for(j = 0; j < num_cols; j++) {
                    dfin.read((char *)&tmp_float, sizeof tmp_float);
                    if (byteOrderNotSame) SwapBytes(tmp_float);
                    data[i][j] = tmp_float;
                }
            }

            // 14. close data file
            dfin.close();
        }
    }
    // ERROR -- UNKNOWN FILE TYPE =============================================
    else {
        std::cerr << "Error in readDataFromFile(): Unknown file type." << std::endl;
        exit(1);
    }
}

void gridFloat::readDataFromGrid(gridFloat & fgrd) {
    long i, j;                                 // loop indices

    // 1.  Check to see if grid size has been set and if so whether if matches the size in the file
    bool gridIsSet = false;
    bool gridIsSameSize = false;
    if ((data == NULL)&&(num_rows == -1)&&(num_cols == -1)) {
        gridIsSet = false;
    }
    else {
        gridIsSet = true;
        if ((doubleValuesWithinTolerance(cell_size, fgrd.cell_size))&&(num_rows == fgrd.num_rows)&&(num_cols == fgrd.num_cols)) {
            if (grid_bbox.equals(fgrd.grid_bbox)) {
                gridIsSameSize = true;
            }
            else {
                gridIsSameSize = false;
            }
        }
        else {
            gridIsSameSize = false;
        }
    }

    // 2. Make any necessary changes to the size of the grid
    if ((gridIsSet)&&(!gridIsSameSize)) {
        // delete old data so can resize grid
        for (i = 0; i < num_rows; i++) {
            delete [] data[i];
            data[i]  = NULL;
        }
        delete [] data;
        gridIsSet = false;
    }

    if (!gridIsSet) {
        // set grid properties
        num_rows = fgrd.num_rows;
        num_cols = fgrd.num_cols;
        cell_size = fgrd.cell_size;
        grid_bbox.reset(fgrd.grid_bbox);
        // create dynamic array
        data = new float*[num_rows];
        for (i = 0; i < num_rows; i++) {
            data[i] = new float[num_cols];
        }
        gridIsSet = true;
    }

    // 3. Read in no data value
    no_dat = fgrd.no_dat;
    

    // 4. Read in data
    for(i = 0; i < num_rows; i++) {
        for(j = 0; j < num_cols; j++) {
            data[i][j] = fgrd.data[i][j];
        }
    }
}

void gridFloat::scaleFromZeroToOne() {
    bool firstNotSet = true;
    double minVal = 0.0;
    double maxVal = 0.0;
    for (long i = 0; i < num_rows; i++) {
        for (long j = 0; j < num_cols; j++) {
            if (data[i][j] != no_dat) {
                if (firstNotSet) {
                    minVal = data[i][j];
                    maxVal = data[i][j];
                }
                else {
                    if (data[i][j] > maxVal) maxVal = data[i][j];
                    if (data[i][j] < minVal) minVal = data[i][j];
                }
            }
        }
    }
    double valRange = maxVal - minVal;
    if (valRange != 0) {
        for (long i = 0; i < num_rows; i++) {
            for (long j = 0; j < num_cols; j++) {
                if (data[i][j] != no_dat) {
                    data[i][j] = (data[i][j] - minVal)/valRange;
                }
            }
        }
    }
}





// ----------------------------------------------------------------------------

// display header file info
void gridFloat::summary() {
    std::cout << "================================================================================" << std::endl;
	std::cout << "Grid ID:" << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout << rasterID << std::endl;
    std::cout << "================================================================================" << std::endl;
    std::cout << "Grid Bounding Box:" << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
	grid_bbox.print();
    std::cout << std::endl;
    std::cout << "================================================================================" << std::endl;
    std::cout << "Grid Properties:" << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout.width(15);
    std::cout << row_s2;
    std::cout.width(15);
    std::cout << num_rows << std::endl;

    std::cout.width(15);
    std::cout << col_s2;
    std::cout.width(15);
    std::cout << num_cols << std::endl;

    std::cout.width(15);
    std::cout << cel_s2;
    std::cout.width(15);
    std::cout << cell_size << std::endl;

    std::cout.width(15);
    std::cout << nan_s3;
    std::cout.width(15);
    std::cout << no_dat << std::endl;
    std::cout << std::endl;
}

void gridFloat::display() {
    summary();
    long i, j;
    std::cout << "================================================================================" << std::endl;
    std::cout << "Grid Data:" << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    for (i = 0; i < num_rows; i++) {
        for (j = 0; j < num_cols; j++) {
            std::cout << data[i][j];
            if (j < (num_cols - 1)) {
                std::cout << "\t";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void gridFloat::printESRIascii(std::string filen) {
    // 1.  open file
    std::ofstream fout;
    fout.open(filen);
    if (!fout.is_open()) {
        std::cerr << "Error in printESRIascii(): Output file "<< filen << " could not be opened." << std::endl;
        exit(1);
    }
    // 2.  number of colums
    fout << std::setw(14) << std::left << col_s2  << num_cols << std::endl;
    // 3.  number of rows
    fout << std::setw(14) << std::left << row_s2 << num_rows << std::endl;
    // 4.  xll corner
    fout << std::setw(14) << std::left << xll_s1B << grid_bbox.get_xMin() << std::endl;
    // 5.  yll corner
    fout << std::setw(14) << std::left << yll_s1B << grid_bbox.get_yMin() << std::endl;
    // 6.  cell size
    fout << std::setw(14) << std::left << cel_s2 << cell_size << std::endl;
    // 7.  no data value
    fout << std::setw(14) << std::left << nan_s1 << no_dat << std::endl;
    // 8.  grid data
    long i, j;
    for (i = 0; i < num_rows; i++) {
        for (j = 0; j < num_cols; j++) {
            fout << data[i][j];
            if (j < (num_cols - 1)) {
                fout << " ";
            }
        }
        fout << std::endl;
    }
    fout << std::endl;
    // 9.  close file
    fout.close();
}

void gridFloat::printESRIbinary(char * filen) {
    // HDR file
    // 1.  open file
    // should we check extention on output file name?
    int name_len = strlen(filen);
    char flt_file[50];                   // char array for index file name
    if ((filen[name_len-4] == '.')&&(filen[name_len-3] == 'h')&&(filen[name_len-2] == 'd')&&(filen[name_len-1] == 'r')) {
        // set binary output file name
        strcpy(flt_file, filen);
        flt_file[name_len - 3] = 'f';
        flt_file[name_len - 2] = 'l';
        flt_file[name_len - 1] = 't';
    }
    else {
        std::cerr << "Error in printESRIbinary(): Output file "<< filen << " does not end in .hdr extension." << std::endl;
        exit(1);
    }
    
    // now, open file
    std::ofstream fout;
    fout.open(filen);
    if (!fout.is_open()) {
        std::cerr << "Error in printESRIbinary(): Output file "<< filen << " could not be opened." << std::endl;
        exit(1);
    }
    // 2.  number of colums
    fout << std::setw(14) << std::left << col_s2 << num_cols << std::endl;
    // 3.  number of rows
    fout << std::setw(14) << std::left << row_s2 << num_rows << std::endl;
    // 4.  xll corner
    fout << std::setw(14) << std::left << xll_s1B << grid_bbox.get_xMin() << std::endl;
    // 5.  yll corner
    fout << std::setw(14) << std::left << yll_s1B << grid_bbox.get_yMin() << std::endl;
    // 6.  cell size
    fout << std::setw(14) << std::left << cel_s2 << cell_size << std::endl;
    // 7.  no data value
    fout << std::setw(14) << std::left << nan_s1 << no_dat << std::endl;
    // 8.  byte order
    if (isMachineBigEndian()) {
        fout << std::setw(14) << std::left << byt_s2 << byt_s1A << std::endl;
    }
    else {
        fout << std::setw(14) << std::left << byt_s2 << byt_s2A << std::endl;
    }
    // 9. close file
    fout.close();
    // FLT FILE
    // 1.  open file
    std::ofstream datout(flt_file, std::ios::out | std::ios::app | std::ios::binary);
    if (!datout.is_open()) {
        std::cerr << "Error in printESRIbinary(): Output file "<< filen << " could not be opened." << std::endl;
        exit(1);
    }
    // 2.  write grid data
    long i, j;
    float tmp_float;
    for (i = 0; i < num_rows; i++) {
        for (j = 0; j < num_cols; j++) {
            tmp_float = data[i][j];
            datout.write((char *)(&tmp_float), sizeof(tmp_float));
        }
    }
    // 3.  close file
    datout.close();
}

/*------------------------------------------------------------------------------
* Constructor to make a 3d raster of no data values                            *
------------------------------------------------------------------------------*/
gridFloat3D::gridFloat3D(const std::vector<double> &xgrid, const std::vector<double> &ygrid, const std::vector<double> &zgrid) {
    xmin = xgrid[0];
    xmax = xgrid[xgrid.size() - 1];
    ymin = ygrid[0];
    ymax = ygrid[ygrid.size() - 1];
    zmin = zgrid[0];
    zmax = zgrid[zgrid.size() - 1];
    xsize = xgrid[1] - xgrid[0];
    ysize = ygrid[1] - ygrid[0];
    zsize = zgrid[1] - zgrid[0];
    xnum = xgrid.size();
    ynum = ygrid.size();
    znum = zgrid.size();

    int array_num = 0;
    for (double i = 0; i < znum; i++) {
        gridFloat *rst = new gridFloat(array_num++, ynum, xnum, xmin, ymin, xsize);
        rst->setAllCellsToZero(true);
        this->xy_grids.push_back(rst);
    }
}

gridFloat3D::~gridFloat3D() {
    for (int i = 0; i < xy_grids.size(); i++) {
        delete xy_grids[i];
    }
}

float gridFloat3D::getGridValue(double eX, double eY, double zcoord) {
    /*
    if(getGridZ(eZ) < zmin || getGridZ(eZ) > zmax) {
        std::cout << "uh hello?" << std::endl;
        return no_dat;
    }
    else {
        return xy_grids[getGridZ(eZ)]->getGridValue(eX, eY);
    } */
    if ((zcoord >= zmin)&&(zcoord <= zmax)) {
        return xy_grids[getGridZ(zcoord)]->getGridValue(eX, eY);
    }
    return no_dat;
}

void gridFloat3D::setGridValue(double eX, double eY, double zcoord, float val) {
    if ((zcoord >= zmin)&&(zcoord <= zmax)) {
        if (xy_grids[getGridZ(zcoord)]->getGridValue(eX, eY) != no_dat)
        xy_grids[getGridZ(zcoord)]->setGridValue(eX, eY, val);
    }
    // else, do nothing
}

void gridFloat3D::addValueToGridCell(double eX, double eY, double zcoord, float val) {
    if ((zcoord >= zmin)&&(zcoord <= zmax)) {
        if (xy_grids[getGridZ(zcoord)]->getGridValue(eX, eY) != no_dat) {
            xy_grids[getGridZ(zcoord)]->addValueToGridCell(eX, eY, val);
        }
        // no data is still no data
    }
    // else, do nothing
}

long gridFloat3D::getGridZ(double zcoord) {
    // check to ensure point is in bounds
    if ((zcoord >= zmin)&&(zcoord <= zmax)) {
        return (long)(floor((zcoord - zmin)/zsize));
    } else {
        std::cerr << "Error in gridFloat3D::getGridZ(): query point out of bounds." << std::endl;
        std::cerr << "Bounding box z min = " << zmin << "; bounding box z max = " << zmax << std::endl;
        std::cerr << "Arg z coordinate = " << zcoord << std::endl;
        exit(1);
    }
}
