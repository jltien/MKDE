//
// Created by Joyce Tien on 11/20/17.
//

#ifndef ANIMALDATA_H
#define ANIMALDATA_H

using namespace std;

struct AnimalData {
public:
    string id;
    vector<double> x;
    vector<double> y;
    vector<double> z;
    vector<double> t;
    vector<bool> use;
    vector<double> mov_err;
    vector<double> obs_err;

    AnimalData(string id_name) {
        id = id_name;
        vector<double> * x = new vector<double>();
        vector<double> * y = new vector<double>();
        vector<double> * z = new vector<double>();
        vector<double> * t = new vector<double>();
        vector<double>  * mov_err = new vector<double>();
        vector<double> * obs_err = new vector<double>();
        vector<bool> * use = new vector<bool>();
    }
};
#endif //ANIMALDATA_H
