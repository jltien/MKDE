//
// Created by Joyce Tien on 11/20/17.
//

#include "pairwiseFunctions.h"
#include <iostream>
#include <vector>

struct AnimalData {
public:
    std::string id;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
    time_t tmin;        // the earliest time, epoch
    time_t tmax;        // the latest time, epoch
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<pointIn3D> xyz;
    std::vector<double> t; // minutes elapsed
    std::vector<struct tm> tm;
    std::vector<time_t> epochSeconds;
    std::vector<bool> use;
    std::vector<double> obsVarXY;
    std::vector<double> obsVarZ;
    std::vector<double> moveVarXY;
    std::vector<double> moveVarZ;


    AnimalData(std::string id_name) {
        id = id_name;
        std::vector<double> * x = new std::vector<double>();
        std::vector<double> * y = new std::vector<double>();
        std::vector<double> * z = new std::vector<double>();
        std::vector<pointIn3D> * xyz = new std::vector<pointIn3D>();
        std::vector<double> * t = new std::vector<double>();
        std::vector<double> * tm = new std::vector<double>();
        std::vector<double> * epochSeconds = new std::vector<double>();
        std::vector<double> * use = new std::vector<double>();
        std::vector<double> * obsVarXY = new std::vector<double>();
        std::vector<double> * obsVarZ = new std::vector<double>();
        std::vector<double> * moveVarXY = new std::vector<double>();
        std::vector<double> * moveVarZ = new std::vector<double>();
    }
};
