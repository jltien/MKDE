#ifndef _MKDE_H
#define _MKDE_H

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_num_threads() 1
    #define omp_get_thread_num() 0
#endif

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <math.h>
#include <random>
#include "Raster/raster.cpp"
#include "AnimalData.h"
#include "Raster/ioFunctions.cpp"

const double MY_PI = 3.141592653589793;
const double MY_EPS = 0.00000001;


gridFloat * mkde2D(const vector<double> &T, const vector<double> &X,
                      const vector<double> &Y, const vector<bool> &use,
                      vector<double> &grid_x, vector<double> &grid_y,
                      vector<double> &move_var, vector<double> &obs_var,
                      double t_step, double pdf_thresh);

gridFloat3D * mkde3dGridv02(vector<double> &T, vector<double> &X,
                             vector<double> &Y, vector<double> &Z, vector<bool> &use,
                             vector<double> &xgrid, vector<double> &ygrid, vector<double> &zgrid,
                             gridFloat zMin, gridFloat zMax, vector<double> &msig2xy,
                             vector<double> &msig2z, vector<double> &osig2xy, vector<double> &osig2z,
                             vector<double> &t_step, vector<double> &pdf_thresh);

int coordToIndex(double x, double minGridCoord, double cellSz);
double univariateNormalDensityThreshold(double p, double sigma2);
long getLinearIndex(long row, long col, long level, long nR, long nC);
int getLowerCellIndex(double minZ, double minGridCoord, double cellSz);
int getUpperCellIndex(double maxZ, double minGridCoord, double cellSz);
double indexToCellCenterCoord(int i, double minGridCoord, double cellSz);
double integrateNormal(double x0, double x1, double mu, double sigma);
double pnorm(double x, double mu, double sigma);
bool isMachineBigEndian();
inline bool doubleEquals(double a, double b);
void withinBounds(AnimalData * animal, long minutes);
void updateTime(AnimalData * animal);

inline bool doubleEquals(double a, double b) {
    return fabs(a - b) < MY_EPS;
}


// IO Functions
unordered_map<string, AnimalData *> *fileRead(const char *in_filename);
void writeMKDE3DtoVTK(vector<double> xgrid, vector<double> ygrid, vector<double> zgrid,
                      gridFloat3D * rst3d, char * filename, char * description);
#endif //MKDE_H
