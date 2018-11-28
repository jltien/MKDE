//
// Created by Joyce Tien on 10/5/18.
//

#ifndef MKDE_PAIRWISEFUNCTIONS_H
#define MKDE_PAIRWISEFUNCTIONS_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <math.h>
#include <random>

struct pointIn3D {
    double x;
    double y;
    double z;
    double alpha_j;
    double time;
    pointIn3D(double xval, double yval, double zval, double t, double a_j) {
        x = xval;
        y = yval;
        z = zval;
        time = t;
        alpha_j = a_j;
    }
    pointIn3D(double xval, double yval, double zval, double t) {
            x = xval;
            y = yval;
            z = zval;
            time = t;
            alpha_j = -1;
    }
};

std::vector<double> assignLocationIndexToTimeGrid(const std::vector<time_t> grid_times,
                                             const std::vector<time_t> location_times, double dt_max);

std::vector<pointIn3D> interpolateCoordinateOnTimeGrid(const std::vector<time_t> & grid_times, const std::vector<double> & location_indexes,
                                                   const std::vector<time_t> & location_times, const std::vector<pointIn3D> & location_xyz);
std::vector<double> euclideanDistance(const std::vector<pointIn3D> & xyz0, const std::vector<pointIn3D> & xyz1, bool use_z);

void printInterpolateAndEuclidean(const std::string filename, const std::vector<pointIn3D> res_xyz, const std::vector<double> res_dist);

void printEuclidean(const std::string filename, const std::vector<double> res_dist);

#endif //MKDE_PAIRWISEFUNCTIONS_H
