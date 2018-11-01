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

struct Tuple {
    double x;
    double y;
    double z;

    Tuple(double xval, double yval, double zval) {
        x = xval;
        y = yval;
        z = zval;
    }
};

std::vector<double> assignLocationIndexToTimeGrid(std::vector<double> & grid_times,
                                             std::vector<double> & location_times, double dt_max);

std::vector<Tuple> interpolateCoordinateOnTimeGrid(std::vector<double> & grid_times, std::vector<double> & location_indexes,
                                                   std::vector<double> & location_times, std::vector<Tuple> & location_xyz);
std::vector<double> euclideanDistance(std::vector<Tuple> & xyz0, std::vector<Tuple> & xyz1, bool use_z = false);

#endif //MKDE_PAIRWISEFUNCTIONS_H
