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

struct gridFloat3D {
public:
    vector<gridFloat *> xy_grids;
    long zgridsize;
    // gridFloat zmin;
    // gridFloat zmax;
    double xmin;
    double ymin;
    double zmin;
    double xmax;
    double ymax;
    double zmax;
    double xsz;
    double ysz;
    double xsize;
    double ysize;
    double zsize;
    long no_dat = -9999;

    gridFloat3D(double xmn, double ymn, double zmn, double xmx, double ymx, double zmx, double xsz, double ysz, double zsz) {
        xmin = xmn;
        ymin = ymn;
        zmin = zmn;
        xmax = xmx;
        ymax = ymx;
        zmax = zmx;
        xsize = xsz;
        ysize = ysz;
        zsize = zsz;
    }

    float getGridValue(double eX, double eY, double eZ) {
	    if(getGridZ(eZ) < zmin || getGridZ(eZ) > zmax) {
            return no_dat;
	    }
	    else {
		    return xy_grids[getGridZ(eZ)]->getGridValue(eX, eY);
	    }
    }

    void setGridValue(double xcoord, double ycoord, double zcoord, float val) {
        if ((zcoord >= zmin)&&(zcoord <= zmax)) {
            xy_grids[getGridZ(zcoord)]->setGridValue(xcoord, ycoord, val);
        }
        // else, do nothing
    }

    void addValueToGridCell(double eX, double eY, double zcoord, float val) {
    if ((zcoord >= zmin)&&(zcoord <= zmax)) {
        if (xy_grids[getGridZ(zcoord)]->getGridValue(eX, eY) != no_dat) {
            xy_grids[getGridZ(zcoord)]->addValueToGridCell(eX, eY, val);
        }
        // no data is still no data
    }
    // else, do nothing
}

    long getGridZ(double zcoord) {
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
};
#endif //ANIMALDATA_H
