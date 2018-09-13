//
// Created by Joyce Tien on 11/20/17.
//

#ifndef ANIMALDATA_H
#define ANIMALDATA_H

using namespace std;

struct AnimalData {
public:
    string id;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
    vector<double> x;
    vector<double> y;
    vector<double> z;
    vector<double> t; // minutes elapsed
    vector<struct tm> tm;
    vector<time_t> epoch_t;
    vector<bool> use;
    vector<double> ObsVarXY;
    vector<double> ObsVarZ;
    vector<double> MoveVarXY;
    vector<double> MoveVarZ;

    AnimalData(string id_name) {
        id = id_name;
        vector<double> * x = new vector<double>();
        vector<double> * y = new vector<double>();
        vector<double> * z = new vector<double>();
        vector<double> * t = new vector<double>();
        vector<double>  * tm = new vector<double>();
        vector<double> * use = new vector<double>();
        vector<double> * ObsVarXY = new vector<double>();
        vector<double> * ObsVarZ = new vector<double>();
        vector<double> * MoveVarXY = new vector<double>();
        vector<double> * MoveVarZ = new vector<double>();
    }
};

#endif //ANIMALDATA_H
