#ifndef IOFUNCTIONS_H
#define IOFUNCTIONS_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <math.h>
#include <random>
#include "AnimalData.h"
#include "raster.h"

using namespace std;

unordered_map<string, AnimalData *> * fileRead(string in_filename);
unordered_map<string, AnimalData *> * fileRead3d(string in_filename);
unordered_map<string, AnimalData *> * fileRead2d(string in_filename);

void withinBounds(AnimalData *animal, long minutes);
void updateTime(AnimalData *animal, int num_args, time_t first_obs);
void writeMKDE3DtoVTK(const vector<double> &xgrid, const vector<double> &ygrid, const vector<double> &zgrid,
                      gridFloat3D * rst3d, string fname, string descr); 
void writeMKDE3DtoGRASS(const vector<double> &xgrid, const vector<double> &ygrid, const vector<double> &zgrid,
                        gridFloat3D * rst3d, string fnm, string nv);
void writeMKDE3DtoXDMF(const vector<double> &xgrid, const vector<double> &ygrid, const vector<double> &zgrid,
                        gridFloat3D * rst3d, string fnmXDMF, string fnmDAT);
void printPoints(gridFloat * grid, string filename);
#endif
