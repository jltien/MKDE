#include "mkde.h"

using namespace std;

// add the err to struct
// get on github, readme with future stuff
// makefile


int main() {

    unordered_map<string, AnimalData *> *animals;
    animals = fileRead("/Users/joycetien/Desktop/SDSC/MKDE/birds.txt");

    /*
    * grid_x / grid_y = array[10] = {0.5, 1.5 ... 9.5}
    * move_var = array[size of x] = {2.0... }
    * obs_var = array[size of x] = {0.1... }
    * t_step = 1.0
    * pdf_thresh = 10^-14
    */

    vector<double> grid_x;
    vector<double> grid_y;
    vector<double> move_var;
    vector<double> obs_var;
    double t_step = 1.0;
    double pdf_thresh = pow(10.0, -14);
    int value = 0.0;

    for (int i = 0; i < 10; i++) {
        grid_x.push_back(value + 0.5);
        grid_y.push_back(value + 0.5);
    }
    for (int i = 0; i < animals->begin()->second->x.size(); i++) {
        move_var.push_back(2.0);
        obs_var.push_back(0.1);
    }

    for (auto it = animals->begin(); it != animals->end(); ++it) {
        mkde2D(it->second->t, it->second->x, it->second->y, it->second->use,
               grid_x, grid_y, move_var, obs_var, t_step, pdf_thresh);
    }

/*
    for (auto it = begin(*animals); it != end(*animals); ++it) {
        for (int i = 0; i < it->second->x.size(); i++) {
            vector<double> x = (it->second)->x;
            vector<double> y = (it->second)->y;
            vector<double> z = (it->second)->z;
            vector<double> t = (it->second)->t;
            cout << it->first << " " << x[i] << " "
                 << y[i] << " " << z[i] << " " << t[i] << endl;
        }
    } */

    return 0;
}

// column of mov and obs err

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

            if (!getline(ss, next, ',')) break;

            record.push_back(next);

            // 5 data entries to populate record with
            if (record.size() != 5) {
                continue;
            }

            string id = record[0];
            double x = stod(record[1]);
            double y = stod(record[2]);
            double z = stod(record[3]);
            double t = stod(record[4]);
            cout << "w " <<record[5] << endl;
            //double merr = stod(record[5]);
           // double oerr = stod(record[6]);
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
        //    new_animal->mov_err.push_back(merr);
          //  new_animal->obs_err.push_back(oerr);
        }
    }

    // Failed to read the file
    if (!infile.eof()) {
        cerr << "FAILED TO READ " << in_filename << endl;
        return nullptr;
    }

    infile.close();

    // Populate the use vector for all the animals, with all trues except the last element
    for (auto it = animals->begin(); it != animals->end(); ++it) {
        for (int i = 0; i < it->second->x.size() - 1; ++i) {
            it->second->use.push_back(true);
        }
        it->second->use.push_back(false);
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
vector<double> mkde2D(const vector<double> &T, const vector<double> &X,
                      const vector<double> &Y, const vector<bool> &use,
                      vector<double> &grid_x, vector<double> &grid_y,
                      vector<double> &move_var, vector<double> &obs_var,
                      double t_step, double pdf_thresh) {

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

                            //last term
                        else if (doubleEquals(t, t1)) {
                            tmpDens = (t - tOld) * pXY / 2.0;
                        }

                            // intermediate terms
                        else {
                            tmpDens = t_step * pXY;
                        }

                        long kk;
                        kk = getLinearIndex(i1, i2, 0, nX, nY);
                        mkde[kk] = mkde[kk] + tmpDens;
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
    long kk;
    for (int i1 = 0; i1 < nX; i1++) {
        for (int i2 = 0; i2 < nY; i2++) {
            kk = getLinearIndex(i1, i2, 0, nX, nY);
            mkde[kk] = mkde[kk] / totalT;
            if (mkde[kk] > maxDens) {
                maxDens = mkde[kk];
            }
            sumDens += mkde[kk];
        }
    }
    cout << "\tMaximum voxel density = " << maxDens << endl;
    cout << "\tSum of voxel densities = " << sumDens << endl;
    cout << "2D MKDE Computation: DONE" << endl;

    return mkde;
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
    return minGridCoord + ((double)i)*cellSz;
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

