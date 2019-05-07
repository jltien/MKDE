#include "ioFunctions.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <math.h>
#include <random>
using namespace std;
const char separator = ',';

/*****************************************************************************
 * delegates to the correct fileRead.
 *****************************************************************************/
unordered_map<string, AnimalData *> * fileRead(string in_filename) {
    int rows = 0;
    ifstream infile(in_filename);
    string line;
    getline(infile, line);
    istringstream ss(line);
    string next;
    while (getline(ss, next, separator)) {
    	rows++;
    }

    if (rows == 9) {
        return fileRead3d(in_filename);
    }
    else if (rows == 6) {
        return fileRead2d(in_filename);
    }
    else {
        cout << "Wrong file format" << endl;
        return 0;
    }
}

/*****************************************************************************
 * Takes an animal data file with an id(string), timestamp, x(double), y(double),
 * z(double), ObsVarXY(double), ObsVarZ(double), MoveVarXY(double),
 * MoveVarZ(double), and parses the data into a vector of vector of the data. The outer
 * vector's elements are unique animals and the inner vector has a single
 * animal's data entries. Returns true if successful, false if unsuccessful.
 * MUST HAVE 9 ARGUMENTS IN THAT ORDER!!!
 *****************************************************************************/
unordered_map<string, AnimalData *> * fileRead3d(string in_filename) {
    unordered_map<string, AnimalData *> *animals = new unordered_map<string, AnimalData *>();
    ifstream infile(in_filename);   // Initialize the file stream
    bool have_header = true;
    AnimalData *new_animal;
    static int num_grid = 0;

    // Keep reading lines until the end of file is reached
    while (infile) {
        string s;

        if (!getline(infile, s)) break;
        // Skip the file header
        if (have_header) {
            have_header = false;
            continue;
        }

        istringstream ss(s);
        vector<string> record;

        // Parses each individual line of data for each data entry
        while (ss) {
            string next;

            if (!getline(ss, next, separator)) break;
            record.push_back(next);
            if (record.size() != 9) {
                continue;
            }
            string id = record[0];
            string date_string = record[1];
            struct tm tm;
            const char *date_char = date_string.c_str();
            strptime(date_char, "%m/%d/%Y %H:%M", &tm);
            time_t epoch = mktime(&tm);
            double x = stod(record[2]);
            double y = stod(record[3]);
            double z = stod(record[4]);
            double obs_var_xy = stod(record[5]);
            double obs_var_z = stod(record[6]);
            double mov_var_xy = stod(record[7]);
            double mov_var_z = stod(record[8]);
            std::unordered_map<std::string, AnimalData *>::const_iterator exists = animals->find(id);

            // A new animal has been encountered
            if (exists == animals->end()) {
                new_animal = new AnimalData(id);
                animals->insert({id, new_animal});
                exists = animals->find(id);
            }

            new_animal = exists->second;
            new_animal->x.push_back(x);
            new_animal->y.push_back(y);
            new_animal->z.push_back(z);
            new_animal->tm.push_back(tm);
            new_animal->epochSeconds.push_back(epoch);
            new_animal->obsVarXY.push_back(obs_var_xy);
            new_animal->obsVarZ.push_back(obs_var_z);
            new_animal->moveVarXY.push_back(mov_var_xy);
            new_animal->moveVarZ.push_back(mov_var_z);
        }
    }

    // Failed to read the file
    if (!infile.eof()) {
        cerr << "FAILED TO READ " << in_filename << endl;
        return nullptr;
    }

    infile.close();

    // find the earliest first observation time
    time_t first_obs = numeric_limits<time_t>::max();
    for (auto it = animals->begin(); it != animals->end(); ++it) {
        if (it->second->epochSeconds[0] < first_obs) {
            first_obs = it->second->epochSeconds[0];
        }
    }

    // set fields for animals
    for (auto it = animals->begin(); it != animals->end(); ++it) {

        // Populate the use vector for all the animals, with all trues except the last element
        for (int i = 0; i < it->second->x.size() - 1; ++i) {
            it->second->use.push_back(true);
        }
        it->second->use.push_back(false);

        // calculate seconds/minutes elapsed since first observation, fill T and use
        updateTime(it->second, 9, first_obs);
        withinBounds(it->second, 43200); // use if within 4320 minutes


        // Find the min and max for x and y
        double xmin = numeric_limits<double>::max();
        double ymin = numeric_limits<double>::max();
        double zmin = numeric_limits<double>::max();
        time_t tmin = numeric_limits<time_t>::max();
        double xmax = 0;
        double ymax = 0;
        double zmax = 0;
        time_t tmax = 0;

        for (int i = 0; i < it->second->x.size(); ++i) {
            if (it->second->x[i] < xmin) {
                xmin = it->second->x[i];
            }
            if (it->second->x[i] > xmax) {
                xmax = it->second->x[i];
            }
            if (it->second->y[i] < ymin) {
                ymin = it->second->y[i];
            }
            if (it->second->y[i] > ymax) {
                ymax = it->second->y[i];
            }
            if (it->second->z[i] < zmin) {
                zmin =it->second->z[i];
            }
            if (it->second->z[i] > zmax) {
                zmax = it->second->z[i];
            }
            if (it->second->epochSeconds[i] < tmin) {
                tmin = it->second->epochSeconds[i];
            }
            if (it->second->epochSeconds[i] > tmax) {
                tmax = it->second->epochSeconds[i];
            }
        }
        it->second->xmin = xmin;
        it->second->xmax = xmax;
        it->second->ymin = ymin;
        it->second->ymax = ymax;
        it->second->zmin = zmin;
        it->second->zmax = zmax;
        it->second->tmin = tmin;
        it->second->tmax = tmax;
    }
    return animals;
}

/*****************************************************************************
 * Takes an animal data file with an id(string), timestamp, x(double), y(double),
 * ObsVarXY(double), MoveVarXY(double), and parses the data into a vector of 
 * vector of the data. The outer vector's elements are unique animals and the 
 * inner vector has a single animal's data entries. Returns true if successful,
 * false if unsuccessful. MUST HAVE 6 ARGUMENTS IN THAT ORDER!!!
 *****************************************************************************/

unordered_map<string, AnimalData *> * fileRead2d(string in_filename) {

    unordered_map<string, AnimalData *> *animals = new unordered_map<string, AnimalData *>();
    ifstream infile(in_filename);   // Initialize the file stream
    bool have_header = true;
    AnimalData *new_animal;
    static int num_grid = 0;

    // Keep reading lines until the end of file is reached
    while (infile) {
        string s;

        if (!getline(infile, s)) break;
        // Skip the file header
        if (have_header) {
            have_header = false;
            continue;
        }

        istringstream ss(s);
        vector<string> record;

        // Parses each individual line of data for each data entry
        while (ss) {
            string next;

            if (!getline(ss, next, separator)) break;
            record.push_back(next);
            if (record.size() != 6) {
                continue;
            }
            string id = record[0];
            string date_string = record[1];
            struct tm tm;
            const char *date_char = date_string.c_str();
            strptime(date_char, "%m/%d/%Y %H:%M", &tm);
            time_t epoch = mktime(&tm);
            double x = stod(record[2]);
            double y = stod(record[3]);
            double obs_var_xy = stod(record[4]);
            double mov_var_xy = stod(record[5]);
            std::unordered_map<std::string, AnimalData *>::const_iterator exists = animals->find(id);

            // A new animal has been encountered
            if (exists == animals->end()) {
                new_animal = new AnimalData(id);
                animals->insert({id, new_animal});
                exists = animals->find(id);
            }

            new_animal = exists->second;
            new_animal->x.push_back(x);
            new_animal->y.push_back(y);
            new_animal->tm.push_back(tm);
            new_animal->epochSeconds.push_back(epoch);
            new_animal->obsVarXY.push_back(obs_var_xy);
            new_animal->moveVarXY.push_back(mov_var_xy);
        }
    }

    // Failed to read the file
    if (!infile.eof()) {
        cerr << "FAILED TO READ " << in_filename << endl;
        return nullptr;
    }

    infile.close();

    // find the earliest first observation time
    time_t first_obs = numeric_limits<time_t>::max();
    for (auto it = animals->begin(); it != animals->end(); ++it) {
        if (it->second->epochSeconds[0] < first_obs) {
            first_obs = it->second->epochSeconds[0];
        }
    }

    for (auto it = animals->begin(); it != animals->end(); ++it) {

        // Populate the use vector for all the animals, with all trues except the last element
        for (int i = 0; i < it->second->x.size() - 1; ++i) {
            it->second->use.push_back(true);
        }
        it->second->use.push_back(false);

        // calculate seconds/minutes elapsed since first observation, fill T and use
        updateTime(it->second, 6, first_obs);

        withinBounds(it->second, 4320); // use if within 4320 minutes


        // Find the min and max for x and y
        double xmin = numeric_limits<double>::max();
        double ymin = numeric_limits<double>::max();
        time_t tmin = numeric_limits<time_t>::max();
        double xmax = 0;
        double ymax = 0;
        time_t tmax = 0;

        for (int i = 0; i < it->second->x.size(); ++i) {
            if (it->second->x[i] < xmin) {
                xmin = it->second->x[i];
            }
            if (it->second->x[i] > xmax) {
                xmax = it->second->x[i];
            }
            if (it->second->y[i] < ymin) {
                ymin = it->second->y[i];
            }
            if (it->second->y[i] > ymax) {
                ymax = it->second->y[i];
            }
            if (it->second->epochSeconds[i] < tmin) {
                tmin = it->second->epochSeconds[i];
            }
            if (it->second->epochSeconds[i] > tmax) {
                tmax = it->second->epochSeconds[i];
            }
        }
        it->second->xmin = xmin;
        it->second->xmax = xmax;
        it->second->ymin = ymin;
        it->second->ymax = ymax;
        it->second->tmin = tmin;
        it->second->tmax = tmax;
        it->second->zmin = -1;
    }

    return animals;
}

/*****************************************************************************
 * Helper functions to check whether the time is within the bounds we're
 * interested in
 *****************************************************************************/
void withinBounds(AnimalData *animal, long minutes) {
    int bound = minutes * 60;
    for (int i = 1; i < animal->t.size(); i++) {
        if (animal->epochSeconds[i] - animal->epochSeconds[i - 1] >= bound) {
            animal->use[i] = false;
            animal->use[i - 1] = false;
        }
    }
}

/*****************************************************************************
 * Helper functions to adjust time relative to the first observation time.
 *****************************************************************************/
void updateTime(AnimalData *animal, int num_args, time_t first_obs) {
    for (int i = 0; i < animal->epochSeconds.size(); i++) {
        animal->epochSeconds[i] = difftime(animal->epochSeconds[i], first_obs);
        animal->t.push_back(animal->epochSeconds[i]/60.0);
    }
	
    for (int i = 0; i < animal->x.size(); i++) {

        // 3d version
    	if (num_args == 9) {
        	animal->xyz.push_back(pointIn3D(animal->x[i], animal->y[i], animal->z[i], animal->t[i]));
        }

        // 2d version
        if (num_args == 6) {
     //    	animal->xyz.push_back(pointIn3D(animal->x[i], animal->y[i], 0, animal->t[i]));
        }
    }
}

/*****************************************************************************
 * Writes MKDE3D to a VTK file format with the name of the file being fname
 *****************************************************************************/
void writeMKDE3DtoVTK(const vector<double> &xgrid, const vector<double> &ygrid, const vector<double> &zgrid,
                      gridFloat3D * rst3d, string fname, string descr) {
    int nX = xgrid.size();
    int nY = ygrid.size();
    int nZ = zgrid.size();

    long nPoints = (long)nX*nY*nZ, ijk = 0;
/*    char * fnm = new char[strlen(fname)+1];
    fnm[strlen(fname)] = 0;
    memcpy(fnm, fname, strlen(fname)); */

    double densTmp;

    std::ofstream vtkfile;
    vtkfile.open (fname);
    // write header
    vtkfile << "# vtk DataFile Version 3.0" << std::endl;
    vtkfile << descr << std::endl;
    vtkfile << "ASCII" << std::endl;
    vtkfile << "DATASET STRUCTURED_GRID" << std::endl;
    vtkfile << "DIMENSIONS " << nX << " " << nY << " " << nZ << std::endl;
    // write points
    vtkfile << "POINTS " << nPoints << " double" << std::endl;
    vtkfile << std::scientific;

    for (int k = 0; k < nZ; k++) {
        for (int j = 0; j < nY; j++) {
            for (int i = 0; i < nX; i++) {
                vtkfile << xgrid[i] << " " << ygrid[j] << " " << zgrid[k] << std::endl;
            }
        }
    }
    // write data
    vtkfile << std::endl << "POINT_DATA " << nPoints << std::endl;
    vtkfile << "SCALARS density double 1" << std::endl;
    vtkfile << "LOOKUP_TABLE densityTable" << std::endl;

    for (double eZ = rst3d->zmin; eZ <= rst3d->zmax; eZ = eZ + rst3d->zsize) {
        for (double eY = rst3d->ymin; eY <= rst3d->ymax; eY = eY + rst3d->ysize) {
            for (long eX = rst3d->xmin; eX <= rst3d->xmax; eX = eX + rst3d->xsize) {
                if (rst3d->getGridValue(eX, eY, eZ) == FLT_NO_DATA_VALUE || rst3d->getGridValue(eX, eY, eZ) <= 0) {
                    vtkfile << "0.000000e+00" << std::endl;
                } else {
                        vtkfile << rst3d->getGridValue(eX, eY, eZ) << std::endl;
                }
            }
        }
    }

    // write lookup table
    vtkfile << std::endl << "LOOKUP_TABLE densityTable 9" << std::endl;
    vtkfile << "255 255 204 0.3" << std::endl;
    vtkfile << "255 237 160 0.4" << std::endl;
    vtkfile << "254 217 118 0.5" << std::endl;
    vtkfile << "254 178 76 0.6" << std::endl;
    vtkfile << "253 141 60 0.7" << std::endl;
    vtkfile << "252 78 42 0.8" << std::endl;
    vtkfile << "227 26 28 0.9" << std::endl;
    vtkfile << "189 0 38 0.9" << std::endl;
    vtkfile << "128 0 38 1.0" << std::endl;
    // done!
    vtkfile.close();
    return;
}



// write to GRASS GIS 3D ASCII raster file
void writeMKDE3DtoGRASS(const vector<double> &xgrid, const vector<double> &ygrid, const vector<double> &zgrid,
                        gridFloat3D * rst3d, string fnm, string nv) {
    int nX = xgrid.size();
    int nY = ygrid.size();
    int nZ = zgrid.size();
    double xSz = xgrid[1] - xgrid[0];
    double ySz = ygrid[1] - ygrid[0];
    double zSz = zgrid[1] - zgrid[0];

    // open file
    std::ofstream r3file;
    r3file.open (fnm);
    r3file << std::setprecision(12);

    // may need to set precision first
    // header
    r3file << "north: " << ygrid[nY - 1] + 0.5 * ySz << std::endl; // north: y.max
    r3file << "south: " << ygrid[0] - 0.5 * ySz << std::endl; // south: y.min
    r3file << "east: " << xgrid[nX - 1] + 0.5 * xSz << std::endl; // east: x.max
    r3file << "west: " << xgrid[0] - 0.5 * xSz << std::endl; // west: x.min
    r3file << "top: " << zgrid[nZ - 1] + 0.5 * zSz << std::endl; // top: z.max
    r3file << "bottom: " << zgrid[0] - 0.5 * zSz << std::endl;// bottom: z.min
    r3file << "rows: " << nY << std::endl; // rows
    r3file << "cols: " << nX << std::endl; // cols
    r3file << "levels: " << nZ << std::endl; // levels: integer

    // data
    // possibly use std::scientific or #include <iomanip> with std::setprecision(12)
    for (int eZ = rst3d->zmin; eZ <= rst3d->zmax; eZ = eZ + rst3d->zsize) {
        for (int eY = rst3d->ymin; eY <= rst3d->ymax; eY = eY + rst3d->ysize) {
            for (int eX = rst3d->xmin; eX <= rst3d->xmax; eX = eX + rst3d->xsize) {
                if (rst3d->getGridValue(eX, eY, eZ) == FLT_NO_DATA_VALUE) {
                    r3file << nv;
                } else {
                    r3file << rst3d->getGridValue(eX, eY, eZ);
                }
                if (eX == rst3d->xmax) {
                    r3file << std::endl;
                } else {
                    r3file << " ";
                }
            }
        }
    }

// done...
    r3file.close();
    return;
}




// write to XDMF file
void writeMKDE3DtoXDMF(const vector<double> &xgrid, const vector<double> &ygrid, const vector<double> &zgrid,
                        gridFloat3D * rst3d, string fnmXDMF, string fnmDAT) {

    int nX = xgrid.size();
    int nY = ygrid.size();
    int nZ = zgrid.size();
    double xSz = xgrid[1] - xgrid[0];
    double ySz = ygrid[1] - ygrid[0];
    double zSz = zgrid[1] - zgrid[0];
    double densTmp;

    int nmSize = 0;
// scan backward through char array, count chars, break when hit "/"; should count \0
    for (int i = fnmDAT.size(); i >= 0; i--) {
        if (fnmDAT[i] == '/') {
            break;
        } else {
            nmSize++;
        }
    }
    char * binName = new char[nmSize];
    int j = 0;
    for (int i = 0; i < nmSize; i++) {
        j = fnmDAT.size() + 1 - nmSize + i; // CHECK THIS!!!!
        binName[i] = fnmDAT[j];
    }

// now copy name to string

// write XML wrapper
    std::ofstream xmffile;
    xmffile.open (fnmXDMF);
    xmffile << std::setprecision(12);

    xmffile << "<?xml version=\"1.0\" ?>" << std::endl;
    xmffile << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
    xmffile << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">" << std::endl;
    xmffile << "<Domain>" << std::endl;
    xmffile << "    <Grid Name=\"Mesh\" GridType=\"Uniform\">" << std::endl;
    xmffile << "        <Topology name=\"topo\" TopologyType=\"3DCoRectMesh\"" << std::endl;
    xmffile << "            Dimensions=\"" << (nZ + 1) << " " << (nY + 1) << " " << (nX + 1) << "\">" << std::endl; // levels, rows, cols?
    xmffile << "        </Topology>" << std::endl;
    xmffile << "        <Geometry name=\"geo\" Type=\"ORIGIN_DXDYDZ\">" << std::endl;
    xmffile << "            <!-- Origin -->" << std::endl;
    xmffile << "            <DataItem Format=\"XML\" Dimensions=\"3\">" << std::endl;
    xmffile << "             " << " " << (xgrid[0] - 0.5*xSz) << " " << (ygrid[0] - 0.5*ySz) << " " << (zgrid[0] - 0.5*zSz) << std::endl;
    xmffile << "            </DataItem>" << std::endl;
    xmffile << "            <!-- DxDyDz -->" << std::endl;
    xmffile << "            <DataItem Format=\"XML\" Dimensions=\"3\">" << std::endl;
    xmffile << "             " << xSz << " " << ySz << " " << zSz <<  std::endl;
    xmffile << "            </DataItem>" << std::endl;
    xmffile << "        </Geometry>" << std::endl;
    xmffile << "        <Attribute Name=\"Density\" Center=\"Cell\">" << std::endl; // need AttributeType="Scalar" or Type="Scalar" ?
    xmffile << "            <DataItem Format=\"Binary\"" << std::endl;
    xmffile << "             DataType=\"Double\"" << std::endl;
    xmffile << "             Precision=\"8\"" << std::endl;
    if (isMachineBigEndian()) {
        xmffile << "             Endian=\"Big\"" << std::endl;
    } else {
        xmffile << "             Endian=\"Little\"" << std::endl;
    }
    xmffile << "             Dimensions=\"" << nZ << " " << nY << " " << nX << "\">" << std::endl;
    xmffile << "               " << binName << std::endl;
    xmffile << "            </DataItem>" << std::endl;
    xmffile << "        </Attribute>" << std::endl;
    xmffile << "    </Grid>" << std::endl;
    xmffile << "</Domain>" << std::endl;
    xmffile << "</Xdmf>" << std::endl;

// close XML file
    xmffile.close();

// write binary data (kji order)
    std::ofstream datfile(fnmDAT, std::ios::out | std::ios::trunc | std::ios::binary); // std::ios::out | std::ios::app | std::ios::binary
    if (!datfile.is_open()) {
        cout << "Error in writeMKDE3DtoXDMF(): Output file "<< fnmDAT << " could not be opened." << std::endl;
    } else {
        for (int eZ = rst3d->zmin; eZ <= rst3d->zmax; eZ = eZ + rst3d->zsize) {
            for (int eY = rst3d->ymin; eY <= rst3d->ymax; eY = eY + rst3d->ysize) {
                for (int eX = rst3d->xmin; eX <= rst3d->xmax; eX = eX + rst3d->xsize) {
                    densTmp = rst3d->getGridValue(eX, eY, eZ);
                    if (densTmp == FLT_NO_DATA_VALUE || densTmp <= 0) {
                        densTmp = 0.0;
                    }
                    datfile.write((char *)(&densTmp), sizeof(densTmp));
                }
            }
        }
        datfile.close();
    }

// done...
    return;
}



void printPoints(gridFloat * grid, std::string filename) {
	std::ofstream file;
    file.open (filename);

	for (double eY = grid->ymin; eY <= grid->ymax; eY = eY + grid->ysize) {
		for (double eX = grid->xmin; eX <= grid->xmax; eX = eX + grid->xsize) {
			if (grid->getGridValue(eX, eY) > 0) {
				file << eX << " " << eY << " " << grid->getGridValue(eX, eY) << std::endl;
			}
		}
	}
	file.close();
}