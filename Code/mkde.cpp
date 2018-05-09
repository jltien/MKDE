#include "mkde.h"

using namespace std;
int xmin = 583022;
int xmax = 722953;
int ymin = 3364040;
int ymax = 3514539;

int main() {
    unordered_map<string, AnimalData *> *animals;
    animals = fileRead("/Users/joycetien/Desktop/SDSC/MKDE/Data/CondorTestData2.txt");

    vector<double> grid_x;
    vector<double> grid_y;
    vector<double> move_var;
    vector<double> obs_var;
    double t_step = 1.0;
    double pdf_thresh = pow(10.0, -14);

    for (int i = xmin; i < xmax; i = i + 250) {
        grid_x.push_back(i);
    }
    for (int i = ymin; i < ymax; i = i + 250) {
        grid_y.push_back(i);
    }

    for (int i = 0; i < animals->begin()->second->x.size(); i++) {
        move_var.push_back(2.0);
        obs_var.push_back(0.1);
    }

    for (auto it = animals->begin(); it != animals->end(); ++it) {
        gridFloat * rst;
        updateTime(it->second);
        withinBounds(it->second, 4320);
        rst = mkde2D(it->second->t, it->second->x, it->second->y, it->second->use,
               grid_x, grid_y, it->second->MoveVarXY, it->second->ObsVarXY, t_step, pdf_thresh);

    }

    return 0;
}

/*****************************************************************************
 * Takes an animal data file with an id(string), x(double), y(double), and
 * z(double) and parses the data into a vector of vector of the data. The outer
 * vector's elements are unique animals and the inner vector has a single
 * animal's data entries. Returns true if successful, false if unsuccessful.
 *****************************************************************************/
unordered_map<string, AnimalData *> *fileRead(const char *in_filename) {
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

            if (!getline(ss, next, '\t')) break;
            record.push_back(next);
            if (record.size() != 10) {
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
            double t = stod(record[5]);
            double obs_var_xy = stod(record[6]);
            double obs_var_z = stod(record[7]);
            double mov_var_xy = stod(record[8]);
            double mov_var_z = stod(record[9]);

            unordered_map<string, AnimalData *>::const_iterator exists = animals->find(id);

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
            new_animal->t.push_back(t);
            new_animal->tm.push_back(tm);
            new_animal->epoch_t.push_back(epoch);
            new_animal->ObsVarXY.push_back(obs_var_xy);
            new_animal->ObsVarZ.push_back(obs_var_z);
            new_animal->MoveVarXY.push_back(mov_var_xy);
            new_animal->MoveVarZ.push_back(mov_var_z);
        }
    }

    // Failed to read the file
    if (!infile.eof()) {
        cerr << "FAILED TO READ " << in_filename << endl;
        return nullptr;
    }

    infile.close();

    for (auto it = animals->begin(); it != animals->end(); ++it) {

        // Populate the use vector for all the animals, with all trues except the last element
        for (int i = 0; i < it->second->x.size() - 1; ++i) {
            it->second->use.push_back(true);
        }
        it->second->use.push_back(false);

        // Find the min and max for x and y
        double xmin = numeric_limits<double>::max();
        double ymin = numeric_limits<double>::max();
        double xmax = 0;
        double ymax = 0;

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
        }
        it->second->xmin = xmin;
        it->second->xmax = xmax;
        it->second->ymin = ymin;
        it->second->ymax = ymax;
        cout << "xmin: " << xmin << "xmax: " << xmax << "ymin: " << ymin << "ymax: " << ymax << endl;
    }

    return animals;
}

/*****************************************************************************
 * Functions that calculate MKDEs for a set of grid cells
 *****************************************************************************/

// 2D CASE
// obsT, obsX, obsY: observed data
// xyTable: a 2d array with pairs of (x,y) coordinates of cell centers
// tMax, tStep, mvSig2xy, obsSig2: parameters
gridFloat * mkde2D(const vector<double> &T, const vector<double> &X,
                      const vector<double> &Y, const vector<bool> &use,
                      vector<double> &grid_x, vector<double> &grid_y,
                      vector<double> &move_var, vector<double> &obs_var,
                      double t_step, double pdf_thresh) {

    static int num_grids = 0;
    long nObs = T.size();

    // grid specs
    int nX = grid_x.size();
    int nY = grid_y.size();
    double xSz = grid_x[1] - grid_x[0];
    double ySz = grid_y[1] - grid_y[0];

    // arrays for MKDE computation
    double *ydens = (double *) malloc(nY * sizeof(double)); // to precompute y
    vector<double> mkde(nX * nY);
    for (int i = 0; i < nX * nY; ++i) {
        mkde[i] = 0.0;
    }

    gridFloat * rst = new gridFloat(num_grids++, nX, nY, xmin, ymin, xSz);
    rst->setAllCellsToZero(true);

    // set up time variables
    double t0, t1, t, tOld, dt, alpha, totalT;

    // set up tmp variables
    double eX, eY;
    double sig2xy, sig2xy_inv;
    double distMaxXY, xyDistSq, voxx, voxy, xyterm;
    int halo1, halo2, i1k, i2k;

    // start computing MKDE
    cout << "2D MKDE Computation: STARTING" << endl;
    totalT = 0.0;
    for (int j = 0; j < (nObs - 1); ++j) {
        cout << "\tProcessing move step " << (j + 1) << " of " << (nObs - 1) << endl;

        // report percent complete after each observation
        t0 = T[j];
        t1 = T[j + 1];
        dt = t1 - t0;
        t = t0;

        if (use[j] == 1) {
            bool exitLoop = false;
            bool finalLoop = false;
            totalT += dt;   // add to total time iff move step used

            // iterate over integration time steps
            while (!exitLoop) {

                // calculate fractional distance between t0 and current time
                alpha = (t - t0) / dt;

                // calculate parameters for bilinear interpolation
                sig2xy = dt * alpha * (1.0 - alpha) * move_var[j] +
                         obs_var[j] * (1.0 - alpha) * (1.0 - alpha) +
                         obs_var[j + 1] * alpha * alpha;
                sig2xy_inv = 1.0 / sig2xy;

                // get (x,y,z) coordinates of kernel origin using linear interpolation
                eX = X[j] + alpha * (X[j + 1] - X[j]);
                eY = Y[j] + alpha * (Y[j + 1] - Y[j]);

                // convert to grid indices
                i1k = coordToIndex(eX, grid_x[0], xSz);
                i2k = coordToIndex(eY, grid_y[0], ySz);

                distMaxXY = univariateNormalDensityThreshold(pdf_thresh, sig2xy);
                halo1 = (int) ceil(distMaxXY / xSz);
                halo2 = (int) ceil(distMaxXY / ySz);

                // precompute voxel density in y dimension
                for (int i2 = 0; i2 < nY; i2++) {
                    ydens[i2] = integrateNormal(grid_y[i2] - 0.5 * ySz, grid_y[i2] + 0.5 * ySz, eY,
                                                sqrt(sig2xy));
                }

                // x-dimension
                for (int i1 = max(0, i1k - halo1); i1 < min(nX, i1k + halo1); i1++) {
                    double voxx = grid_x[i1];   // voxel x
                    double xdens = integrateNormal(voxx - 0.5 * xSz, voxx + 0.5 * xSz, eX, sqrt(sig2xy));

                    for (int i2 = max(0, i2k - halo2); i2 < min(nY, i2k + halo2); i2++) {
                        double voxy = grid_y[i2];   // voxel y
                        double pXY = xdens * ydens[i2];
                        double tmpDens;

                        // first term
                        if (doubleEquals(t, t0)) {
                            tmpDens = t_step * pXY / 2.0;
                        }

                        // last term
                        else if (doubleEquals(t, t1)) {
                            tmpDens = (t - tOld) * pXY / 2.0;
                        }

                        // intermediate term
                        else {
                            tmpDens = t_step * pXY;
                        }

                        long kk;
                        kk = getLinearIndex(i1, i2, 0, nX, nY);
                        mkde[kk] = mkde[kk] + tmpDens;
                        rst->addValueToGridCell(eX, eY, tmpDens);
                    }
                }

                // update the eval time (t) and stopping conditions
                if (finalLoop) {
                    exitLoop = true;
                } else {
                    tOld = t;
                    t += t_step;
                    if (t >= t1) {
                        t = t1;
                        finalLoop = true;
                    }
                }
            }
        }
    }

    double maxDens = 0.0, sumDens = 0.0;

    for (double eX = xmin; eX < xmax; eX=eX+250) {
        for (double eY = ymin; eY < ymax; eY=eY+250) {
            if (rst->getGridValue(eX, eY) != FLT_NO_DATA_VALUE) {
                rst->setGridValue(eX, eY, rst->getGridValue(eX, eY) / totalT);
                if (rst->getGridValue(eX, eY) > maxDens) {
                    maxDens = rst->getGridValue(eX, eY);
                }
                sumDens += rst->getGridValue(eX, eY);
            }
        }
    }
    cout << "\tMaximum voxel density = " << maxDens << endl;
    cout << "\tSum of voxel densities = " << sumDens << endl;
    cout << "2D MKDE Computation: DONE" << endl;
    rst->printESRIascii("raster_test.asc");
    return rst;

}


/*
vector<double> mkde3dGridv02(vector<double> &T, vector<double> &X,
                             vector<double> &Y, vector<double> &Z, vector<bool> &use,
                             vector<double> &xgrid, vector<double> &ygrid, vector<double> &zgrid,
                             SEXP zMinMatrix, SEXP zMaxMatrix, vector<double> &msig2xy,
                             vector<double> &msig2z, vector<double> &osig2xy, vector<double> &osig2z,
                             vector<double> &t_step, vector<double> &pdf_thresh) {

    int nObs = T.size();

    // grid speces
    int nX = xgrid.size();
    int nY = ygrid.size();
    int nZ = zgrid.size();
    long nVoxels = (long) nX * nY * nZ;
    double xSz = xgrid[1] - xgrid[0];
    double ySz = ygrid[1] - ygrid[0];
    double zSz = zgrid[1] - zgrid[0];

    // pysical constraints in z-dimension
    Rcpp::NumericMatrix zMin(zMinMatrix);  // lower z-limit at each (x,y)
    Rcpp::NumericMatrix zMax(zMaxMatrix);  // upper z-limit at each (x,y)

    // arrays for MKDE computation
    double *ydens = (double *) malloc(nY * sizeof(double)); // To precompute Y
    double *zdens = (double *) malloc(nZ * sizeof(double)); // To precompute z
    vector<double> mkde(nVoxels);
    for (long i = 0; i < nVoxels; i++) {
        mkde[i] = 0.0;
    }

    // set up time variables
    double t0, t1, t, tOld, dt, alpha;

    // set up tmp variables
    double eX, eY, eZ, factor;
    double W;
    double sig2xy, sig2z, sig2xy_inv, sig2z_inv;
    double distMaxXY, distMaxZ, xyDistSq, zDistSq, tmpZ, xyterm;
    int halo1, halo2, halo3, i1k, i2k, i3k;

    // FOLLOWING FOR DEBUGGING
    int i1min = nX, i1max = 0, i2min = nY, i2max = 0, i3min = nZ, i3max = 0;

    // start computing MKDE
    cout << "3D MKDE Computation: STARTING" << endl;
    W = 0.0;
    for (int j = 0; j < (nObs - 1); j++) {
        cout << "\tProcessing move step " << (j + 1) << " of " << (nObs - 1) << endl;

        // report percent complete after each observation
        t0 = T[j];
        t1 = T[j + 1];
        dt = t1 - t0;
        t = t0;
        if (use[j] == 1) {
            bool exitLoop = false;
            bool finalLoop = false;
            while (!exitLoop) { // iterate over integration time steps

                // Calculate fractional distance between t0 and current time
                alpha = (t - t0) / dt;

                // Calculate parameters for bilinear interpolation
                sig2xy = dt * alpha * (1.0 - alpha) * msig2xy[j] +
                         osig2xy[j] * (1.0 - alpha) * (1.0 - alpha) +
                         osig2xy[j + 1] * alpha * alpha;
                sig2z = dt * alpha * (1.0 - alpha) * msig2z[j] +
                        osig2z[j] * (1.0 - alpha) * (1.0 - alpha) +
                        osig2z[j + 1] * alpha * alpha;
                // Get (x,y,z) coordinates of kernel origin using linear interpolation
                eX = X[j] + alpha * (X[j + 1] - X[j]);
                eY = Y[j] + alpha * (Y[j + 1] - Y[j]);
                eZ = Z[j] + alpha * (Z[j + 1] - Z[j]);
                // Convert to grid indices
                i1k = coordToIndex(eX, xgrid[0], xSz);
                i2k = coordToIndex(eY, ygrid[0], ySz);
                i3k = coordToIndex(eZ, zgrid[0], zSz);

                distMaxXY = univariateNormalDensityThreshold(pdf_thresh[0], sig2xy);
                halo1 = ceil(distMaxXY / xSz);
                halo2 = ceil(distMaxXY / ySz);
                distMaxZ = univariateNormalDensityThreshold(pdf_thresh[0], sig2z);
                halo3 = ceil(distMaxZ / zSz);

                // Precompute voxel density in y dimension
                for (int i2 = 0; i2 < nY; i2++) {
                    ydens[i2] = integrateNormal(ygrid[i2] - 0.5 * ySz, ygrid[i2] + 0.5 * ySz, eY, sqrt(sig2xy));
                }
                // Precompute voxel density in z dimension
                for (int i3 = 0; i3 < nZ; i3++) {
                    zdens[i3] = integrateNormal(zgrid[i3] - 0.5 * zSz, zgrid[i3] + 0.5 * zSz, eZ, sqrt(sig2z));
                }

                for (int i1 = std::max(0, i1k - halo1); i1 < std::min(nX, i1k + halo1); i1++) { // x-dimension
                    double voxx = xgrid[i1]; // voxel x
                    double xdens = integrateNormal(voxx - 0.5 * xSz, voxx + 0.5 * xSz, eX, sqrt(sig2xy));
                    for (int i2 = std::max(0, i2k - halo2); i2 < std::min(nY, i2k + halo2); i2++) { // y-dimension
                        double voxy = ygrid[i2]; // voxel y
                        // get the range of indexes and coordinates based on the physical boundaries
                        int i3lo = std::max(0, getLowerCellIndex(zMin(i1, i2), zgrid[0], zSz));
                        int i3hi = std::min(nZ, getUpperCellIndex(zMax(i1, i2), zgrid[0], zSz) +
                                                1); // add 1 because less than
                        double loZ = indexToCellCenterCoord(i3lo, zgrid[0], zSz) - 0.5 * zSz;
                        double hiZ = indexToCellCenterCoord(i3hi - 1, zgrid[0], zSz) + 0.5 * zSz;
                        // Reflect E[Z] about the lower and upper boundaries
                        double loReflZ = 2.0 * loZ - eZ;
                        double hiReflZ = 2.0 * hiZ - eZ;

                        for (int i3 = std::max(i3lo, i3k - halo3);
                             i3 < std::min(i3hi, i3k + halo3); i3++) { // z-dimension

                            int th_id = omp_get_thread_num();
                            // set up for reflection
                            // only compute if the expected location is within distance of boundary
                            double voxz = zgrid[i3];
                            double loDZ = 0.0;
                            if (std::fabs(eZ - loZ) <= distMaxZ) {
                                loDZ = integrateNormal(voxz - 0.5 * zSz, voxz + 0.5 * zSz, loReflZ, sqrt(sig2z));
                            }
                            double hiDZ = 0.0;
                            if (std::fabs(hiZ - eZ) <= distMaxZ) {
                                hiDZ = integrateNormal(voxz - 0.5 * zSz, voxz + 0.5 * zSz, hiReflZ, sqrt(sig2z));
                            }

                            // Calculate contribution of kernel to voxel
                            double pXYZ = xdens * ydens[i2] * (zdens[i3] + loDZ + hiDZ);
                            // update density (trapezoid rule)
                            double tmpDens;
                            if (doubleEquals(t, t0)) { // first term
                                tmpDens = t_step[0] * pXYZ / 2.0;
                            } else if (doubleEquals(t, t1)) { // last term
                                tmpDens = (t - tOld) * pXYZ / 2.0;
                            } else { // intermediate terms
                                tmpDens = t_step[0] * pXYZ;
                            }

                            // Add contribution to voxel (removed Kahan summation for now)
                            long kk = getLinearIndex(i1, i2, i3, nX, nY);
                            mkde[kk] += tmpDens;
                            W += tmpDens;
                        }
                    }
                }
                // end parallel for

                // update the eval time (t) and stopping conditions
                if (finalLoop) {
                    exitLoop = true;
                } else {
                    tOld = t;
                    t += t_step[0];
                    if (t >= t1) {
                        t = t1;
                        finalLoop = true;
                    }
                }
            }
        }
    }

    // divide by totalT
    double maxDens = 0.0, sumDens = 0.0;
    long kk;
    for (int i1 = 0; i1 < nX; i1++) {
        for (int i2 = 0; i2 < nY; i2++) {
            for (int i3 = 0; i3 < nZ; i3++) {
                kk = getLinearIndex(i1, i2, i3, nX, nY);
                mkde[kk] = mkde[kk] / W;
                if (mkde[kk] > maxDens) {
                    maxDens = mkde[kk];
                }
                sumDens += mkde[kk];
            }
        }
    }
    cout << "\tMaximum voxel density = " << maxDens << std::endl;
    cout << "\tSum of voxel densities = " << sumDens << std::endl;
    cout << "3D MKDE Computation: DONE" << std::endl;

    // RETURN DENSITY HERE
    return mkde;
}
*/

/*****************************************************************************
 * Helper functions to translate between cell or voxel indexes and coordinates
 *****************************************************************************/
int coordToIndex(double x, double minGridCoord, double cellSz) {
    return (int) floor((x - minGridCoord - 0.5 * cellSz) / cellSz);
}

long getLinearIndex(long row, long col, long level, long nR, long nC) {
    return row + col * nR + level * nR * nC;
}

double indexToCellCenterCoord(int i, double minGridCoord, double cellSz) {
    return minGridCoord + ((double) i) * cellSz;
}

int getLowerCellIndex(double minZ, double minGridCoord, double cellSz) {
    int i = coordToIndex(minZ, minGridCoord, cellSz);
    double cellZ = indexToCellCenterCoord(i, minGridCoord, cellSz);
    if (cellZ < minZ) {
        i++;
    }
    return i;
}

int getUpperCellIndex(double maxZ, double minGridCoord, double cellSz) {
    int i = coordToIndex(maxZ, minGridCoord, cellSz);
    double cellZ = indexToCellCenterCoord(i, minGridCoord, cellSz);
    if (cellZ > maxZ) {
        i--;
    }
    return i;
}

/*****************************************************************************
 * Helper functions to determine range of cells or voxels over which to
 * update density.
 *****************************************************************************/
double univariateNormalDensityThreshold(double p, double sigma2) {
    double z = -2.0 * sigma2 * log(p * sqrt(2.0 * sigma2 * MY_PI));
    double res = NAN;
    if (z >= 0.0) {
        res = sqrt(z);
    }
    return res;
}

/*****************************************************************************
 * Function to find the cumulative distribution function.
 *****************************************************************************/
double pnorm(double x, double mu, double sigma) {
    double err = erf(((x - mu) / (sigma * sqrt(2))));
    return 0.5 * (1 + err);
}

/*****************************************************************************
 * Helper functions for truncated normal kernel
 *****************************************************************************/
// Uses R standalone math library
// get area under normal PDF from x0 to x1 with mean mu and standard deviation sigma
double integrateNormal(double x0, double x1, double mu, double sigma) {
    double p0 = pnorm(x0, mu, sigma);
    double p1 = pnorm(x1, mu, sigma);
    return p1 - p0;
}

bool isMachineBigEndian(void) {
    union {
        long num;
        unsigned char uc[sizeof(long)];
    } u;
    u.num = 1;
    if (u.uc[0] == 1) {
        return false;
    } else {
        return true;
    }
}

/*****************************************************************************
 * Helper functions to check whether the time is within the bounds we're
 * interested in
 *****************************************************************************/
void withinBounds(AnimalData * animal, long minutes) {
    int bound = minutes * 60;
    for (int i = 1; i < animal->t.size(); i++) {
        if (animal->epoch_t[i] - animal->epoch_t[i - 1] >= bound) {
            animal->use[i] = false;
            animal->use[i - 1] = false;
        }
    }
}

/*****************************************************************************
 * Helper functions to adjust time relative to the first observation time.
 *****************************************************************************/
void updateTime(AnimalData * animal) {
    for (int i = 0; i < animal->epoch_t.size(); i++) {
        animal->epoch_t[i] = animal->epoch_t[i] - animal->epoch_t[0];
    }
}