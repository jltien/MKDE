#include "mkde.h"

using namespace std;
/*
int xmin = 583022;
int xmax = 722953;
int ymin = 3364040;
int ymax = 3514539;
*/
int main() {
    unordered_map<string, AnimalData *> *animals;
    animals = fileRead("/Users/joycetien/Desktop/SDSC/MKDE/Data/CondorTestData2.txt");

    vector<double> grid_x;
    vector<double> grid_y;
    vector<double> grid_z;
    vector<double> move_var;
    vector<double> obs_var;
    gridFloat * rst;
    gridFloat3D * rst3d;

    double t_step = 1.0;
    double pdf_thresh = pow(10.0, -14);

    vector<double> t_step3d;
    vector<double> pdf_thresh3d;

    // find the lowest and highest x, y, z for all animals
    double xmin = numeric_limits<double>::max();
    double xmax = 0;
    double ymin = numeric_limits<double>::max();
    double ymax = 0;
    double zmin = numeric_limits<double>::max();
    double zmax = 0;

    for (auto it = animals->begin(); it != animals->end(); ++it) {
        if (it->second->xmin < xmin) {
            xmin = it->second->xmin;
        }
        if (it->second->xmax > xmax) {
            xmax = it->second->xmax;
        }
        if (it->second->ymin < ymin) {
            ymin = it->second->ymin;
        }
        if (it->second->ymax > ymax) {
            ymax = it->second->ymax;
        }
        if (it->second->zmin < zmin) {
            zmin = it->second->zmin;
        }
        if (it->second->zmax > zmax) {
            zmax = it->second->zmax;
        }
    }

    // set the x, y, z grids
    for (double i = xmin - 2000; i < xmax + 2000; i = i + 250) {
        grid_x.push_back(i);
    }
    for (double i = ymin - 2000; i < ymax + 2000; i = i + 250) {
        grid_y.push_back(i);
    }

    for (double i = 0; i < zmax + 2000; i = i + 250) {
        grid_z.push_back(i);
    }

    for (int i = 0; i < animals->begin()->second->x.size(); i++) {
        move_var.push_back(2.0);
        obs_var.push_back(0.1);
        t_step3d.push_back(1.0);
        pdf_thresh3d.push_back(pow(10.0, -14));
    }

    gridFloat * high = new gridFloat(1, grid_x.size(), grid_y.size(), xmin, ymin, grid_x[1]-grid_x[0]);

    high->setAllCellsTo(zmax);
    gridFloat * low = new gridFloat(0, grid_x.size(), grid_y.size(), xmin, ymin, grid_x[1]-grid_x[0]);

    low->setAllCellsToZero(true);

  //  for (auto it = animals->begin(); it != animals->end(); ++it) {
        auto it = animals->begin();
        updateTime(it->second);
        withinBounds(it->second, 4320);
        /*rst = mkde2D(it->second->t, it->second->x, it->second->y, it->second->use,
               grid_x, grid_y, it->second->MoveVarXY, it->second->ObsVarXY, t_step, pdf_thresh);
        */

        rst3d = mkde3dGridv02(it->second->t, it->second->x, it->second->y, it->second->z, it->second->use,
                      grid_x, grid_y, grid_z, low, high, it->second->MoveVarXY,
                      it->second->MoveVarZ, it->second->ObsVarXY, it->second->ObsVarZ, t_step3d,
                      pdf_thresh3d);

        writeMKDE3DtoVTK(grid_x, grid_y, grid_z, rst3d, "rastertest.vtk", "test data");
//    }

    //rst->printESRIascii("raster_test.asc");
    return 0;
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
    double xmin = grid_x[0];
    double ymin = grid_y[0];
    double xmax = grid_x[grid_x.size() - 1];
    double ymax = grid_y[grid_y.size() - 1];

    // arrays for MKDE computation
    double *ydens = (double *) malloc(nY * sizeof(double)); // to precompute y

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

    for (double eX = xmin; eX < xmax; eX = eX + xSz) {
        for (double eY = ymin; eY < ymax; eY = eY + ySz) {
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
    return rst;

}


gridFloat3D * mkde3dGridv02(const vector<double> &T, const vector<double> &X,
                             const vector<double> &Y, const vector<double> &Z, const vector<bool> &use,
                             const vector<double> &xgrid, const vector<double> &ygrid, const vector<double> &zgrid,
                             const gridFloat * zMin, const gridFloat * zMax, const vector<double> &msig2xy,
                             const vector<double> &msig2z, const vector<double> &osig2xy, const vector<double> &osig2z,
                             const vector<double> &t_step, const vector<double> &pdf_thresh) {

    int nObs = T.size();

    // grid speces
    int nX = xgrid.size();
    int nY = ygrid.size();
    int nZ = zgrid.size();
    long nVoxels = (long) nX * nY * nZ;
    double xSz = xgrid[1] - xgrid[0];
    double ySz = ygrid[1] - ygrid[0];
    double zSz = zgrid[1] - zgrid[0];
    double xmin = xgrid[0];
    double xmax = xgrid[xgrid.size() - 1];
    double ymin = ygrid[0];
    double ymax = ygrid[ygrid.size() - 1];
    double zmin = zgrid[0];
    double zmax = zgrid[zgrid.size() - 1];

    // arrays for MKDE computation
    double *ydens = (double *) malloc(nY * sizeof(double)); // To precompute Y
    double *zdens = (double *) malloc(nZ * sizeof(double)); // To precompute z

    // Create a vector of GridFloats and initializing GridFloat3D
    gridFloat3D * rst3d = new gridFloat3D(xgrid, ygrid, zgrid);

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
                        int i3lo = std::max(0, getLowerCellIndex(/*zMin.getGridValue(xmin + (double)i1, ymin + (double)i2)*/ 0, zgrid[0], zSz));
                        int i3hi = std::min(nZ, getUpperCellIndex(/*zMax.getGridValue(xmin + (double)i1, ymin + (double)i2)*/ 5000, zgrid[0], zSz) +
                                                1); // add 1 because less than
                        double loZ = indexToCellCenterCoord(i3lo, zgrid[0], zSz) - 0.5 * zSz;
                        double hiZ = indexToCellCenterCoord(i3hi - 1, zgrid[0], zSz) + 0.5 * zSz;
                        // Reflect E[Z] about the lower and upper boundaries
                        double loReflZ = 2.0 * loZ - eZ;
                        double hiReflZ = 2.0 * hiZ - eZ;
                        for (int i3 = std::max(i3lo, i3k - halo3); i3 < std::min(i3hi, i3k + halo3); i3++) { // z-dimension
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

                            rst3d->addValueToGridCell(eX, eY, eZ, tmpDens);
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
    for (double eZ = zmin; eZ < zmax; eZ = eZ + zSz) {
        for (double eY = ymin; eY < ymax; eY = eY + ySz) {
            for (long eX = xmin ;eX < xmax; eX = eX + xSz) {
                if (rst3d->getGridValue(eX, eY, eZ) != FLT_NO_DATA_VALUE) {
                    rst3d->setGridValue(eX, eY, eZ, rst3d->getGridValue(eX, eY, eZ) / W);
                    if (rst3d->getGridValue(eX, eY, eZ) > maxDens) {
                        maxDens = rst3d->getGridValue(eX, eY, eZ);
                    }
                    sumDens += rst3d->getGridValue(eX, eY, eZ);
                }
            }
        }
    }

/*
        for (int eZ = zmin; eZ < zmax; eZ = eZ + zSz) {
            for (int eY = ymin; eY <ymax; eY = eY + ySz) {
                for (int eX = xmin; eX < xmax; eX = eX + xSz) {
                    cout << "i: " << eX << "j: " << eY << "k: " << eZ << rst3d->getGridValue(eX, eY , eZ) << endl;
                }
            }
        } */
    cout << "\tMaximum voxel density = " << maxDens << std::endl;
    cout << "\tSum of voxel densities = " << sumDens << std::endl;
    cout << "3D MKDE Computation: DONE" << std::endl;

    // RETURN DENSITY HERE
    return rst3d;
}

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