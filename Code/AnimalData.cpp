//
// Created by Joyce Tien on 11/20/17.
//

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
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> t; // minutes elapsed
    std::vector<struct tm> tm;
    std::vector<time_t> epoch_t;
    std::vector<bool> use;
    std::vector<double> ObsVarXY;
    std::vector<double> ObsVarZ;
    std::vector<double> MoveVarXY;
    std::vector<double> MoveVarZ;

    AnimalData(std::string id_name) {
        id = id_name;
        std::vector<double> * x = new std::vector<double>();
        std::vector<double> * y = new std::vector<double>();
        std::vector<double> * z = new std::vector<double>();
        std::vector<double> * t = new std::vector<double>();
        std::vector<double>  * tm = new std::vector<double>();
        std::vector<double> * use = new std::vector<double>();
        std::vector<double> * ObsVarXY = new std::vector<double>();
        std::vector<double> * ObsVarZ = new std::vector<double>();
        std::vector<double> * MoveVarXY = new std::vector<double>();
        std::vector<double> * MoveVarZ = new std::vector<double>();
    }
};