#include "mkde.h"

int main() {
/* test code for pairwise functions
    vector<double> tgrid;
    for (double i = 50.0; i <= 200.0; i++) {
        tgrid.push_back(i);
    }
    vector<double> tloc_a = {75.2, 80.1, 85.2, 90.2, 95.1, 100.3, 105.2, 115.1, 120.2, 125.1, 130.3, 140.3, 145.1,
            150.3, 155.1, 160.2, 180.1, 185.2, 190.2, 195.1, 200.2, 205.1};
    vector<double> iloc_a = assignLocationIndexToTimeGrid(tgrid, tloc_a, 10000);
    for (int i = 0; i < iloc_a.size(); i++) {
        cout << iloc_a[i] << endl;
    }
    */
/*
    unordered_map<string, AnimalData *> *animals;
    animals = fileRead("/Users/joycetien/Desktop/test/samepoints.txt");

    vector<AnimalData *> data;

    for (auto it = animals->begin(); it != animals->end(); it++) {
        vector<time_t> tgrid;
        for (time_t i = it->second->tmin; i <= it->second->tmax; i = i + 3600) {
            tgrid.push_back(i);
        }
        data.push_back(it->second);
        vector<double> ans = assignLocationIndexToTimeGrid(tgrid, it->second->epoch_seconds, 1000000);

        for (int i = 0; i < ans.size(); i++) {
            cout << "ans: " << ans[i] << endl;
        }

        vector<pointIn3D> ans1 = interpolateCoordinateOnTimeGrid(tgrid, ans, it->second->epoch_seconds, it->second->xyz);
        for (int i = 0; i < ans.size(); i++) {
            cout << "ans1: " << ans1[i].x << " " << ans1[i].y << " " << ans1[i].z << endl;
        }

    } */
    /*
    vector<time_t> tgrid;
    for (time_t i = 1; i < 12; i++) {
        tgrid.push_back(i);
    }
    vector<time_t> location;
    location.push_back(2.5);
    location.push_back(4);
    location.push_back(7);
    vector<double> ans = assignLocationIndexToTimeGrid(tgrid, location, 1000);
    for (int i = 0; i < ans.size(); i++) {
        cout << "tgrid: " << tgrid[i] <<"   ans: " << ans[i] << endl;
    } */
/*
    vector<double> dist = euclideanDistance(data[0]->xyz, data[1]->xyz, false);

    for (int i = 0; i < dist.size(); i++) {
        cout << dist[i] << endl;
    }
*/

    unordered_map<string, AnimalData *> *animals;
    animals = fileRead("/Users/joycetien/Desktop/SDSC/MKDE/Data/CondorTestData2.txt");

    vector<double> grid_x;
    vector<double> grid_y;
    vector<double> grid_z;
    vector<double> move_var;
    vector<double> obs_var;
    gridFloat * rst;
    gridFloat3D *rst3d;

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

    for (double i = 0; i < zmax + 2000; i = i + 500) {
        grid_z.push_back(i);
    }

    for (int i = 0; i < animals->begin()->second->x.size(); i++) {
        move_var.push_back(2.0);
        obs_var.push_back(0.1);
        t_step3d.push_back(1.0);
        pdf_thresh3d.push_back(pow(10.0, -14));
    }

    // set up zmin and zmax
    gridFloat *high = new gridFloat(1, grid_x.size(), grid_y.size(), xmin, ymin, grid_x[1] - grid_x[0]);
    high->setAllCellsTo(zmax);
    gridFloat *low = new gridFloat(0, grid_x.size(), grid_y.size(), xmin, ymin, grid_x[1] - grid_x[0]);
    low->setAllCellsToZero(true);

    // code for testing regular mkde2D and mkde3D
    for (auto it = animals->begin(); it != animals->end(); ++it) {
        updateTime(it->second);
        withinBounds(it->second, 4320);
//        rst = mkde2D(it->second->t, it->second->x, it->second->y, it->second->use,
//               grid_x, grid_y, it->second->MoveVarXY, it->second->ObsVarXY, t_step, pdf_thresh);


        rst3d = mkde3dGridv02(it->second->t, it->second->x, it->second->y, it->second->z, it->second->use,
                              grid_x, grid_y, grid_z, low, high, it->second->MoveVarXY,
                              it->second->MoveVarZ, it->second->ObsVarXY, it->second->ObsVarZ, t_step3d,
                              pdf_thresh3d);


        string filename = it->first + ".vtk";
        writeMKDE3DtoVTK(grid_x, grid_y, grid_z, rst3d, filename, "test data");


//        string filename = it->first + ".grass";
//        writeMKDE3DtoGRASS(grid_x, grid_y, grid_z, rst3d, filename, "0");

//        string filename = it->first + ".xdmf";
//        writeMKDE3DtoXDMF(grid_x, grid_y, grid_z, rst3d, "test", filename);

    }

    /*
     * interaction test code
     */
    vector<AnimalData *> avec;
    for (auto it = animals->begin(); it != animals->end(); ++it) {
        avec.push_back(it->second);
    }
    rst = mkde2dGridv02interact(avec[0]->t, avec[0]->x, avec[0]->y, avec[1]->x, avec[1]->y, avec[0]->use,
                                grid_x, grid_y, avec[0]->MoveVarXY, avec[1]->MoveVarXY, avec[0]->ObsVarXY, avec[1]->ObsVarXY,
                                t_step3d, pdf_thresh3d);
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

    gridFloat *rst = new gridFloat(num_grids++, nX, nY, xmin, ymin, xSz);
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
                           gridFloat *zMin, gridFloat *zMax, const vector<double> &msig2xy,
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
    gridFloat3D *rst3d = new gridFloat3D(xgrid, ygrid, zgrid);

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
                        int i3lo = std::max(0,
                                            getLowerCellIndex(/*zMin->getGridValue(xmin + (double)i1, ymin + (double)i2)*/
                                                    0, zgrid[0], zSz));
                        int i3hi = std::min(nZ,
                                            getUpperCellIndex(/*zMax.getGridValue(xmin + (double)i1, ymin + (double)i2)*/
                                                    5000, zgrid[0], zSz) +
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
    for (double eZ = zmin; eZ <= zmax; eZ = eZ + zSz) {
        for (double eY = ymin; eY <= ymax; eY = eY + ySz) {
            for (long eX = xmin; eX <= xmax; eX = eX + xSz) {
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


/* this version computes the probability that the two individuals will be in the
   same cell, not at the exact same location.  To compute the probability of
   occurrence at the same point within a cell, we will have to integrate the
   product of the two kernels over the area of the cell */
gridFloat * mkde2dGridv02interact(const vector<double> &T, const vector<double> &X0, vector<double> &Y0,
                                 const vector<double> &X1, const vector<double> &Y1, const vector<bool> &isValid,
                                 const vector<double> &xGrid, const vector<double> &yGrid,
                                 const vector<double> &msig2xy0, const vector<double> &msig2xy1,
                                 const vector<double> &osig2xy0, const vector<double> &osig2xy1,
                                 const vector<double> &stepT, const vector<double> &pdfMin) {

    int nObs = T.size();
    static int num_grids = 0;

    // grid speces
    int nX = xGrid.size();
    int nY = yGrid.size();
    long nCells = (long) nX * nY;
    double xSz = xGrid[1] - xGrid[0];
    double ySz = yGrid[1] - yGrid[0];
    double xmin = xGrid[0];
    double ymin = yGrid[0];
    double xmax = xGrid[xGrid.size() - 1];
    double ymax = yGrid[yGrid.size() - 1];

    // arrays for MKDE computation
    double *yprob0 = (double *) malloc(nY * sizeof(double)); // To PRECOMPUTE Y
    double *yprob1 = (double *) malloc(nY * sizeof(double)); // To PRECOMPUTE Y
    double *ybhattdist = (double *) malloc(nY * sizeof(double)); // To PRECOMPUTE Y

    gridFloat * mkde = new gridFloat(num_grids++, nX, nY, xmin, ymin, xSz);
    mkde->setAllCellsToZero(true);

    // set up time variables
    double t0, t1, t, tOld, dt, alpha;
    // set up tmp variables
    double eX0, eY0, eX1, eY1;
    double W0 = 0.0, W1 = 0.0;
    double totalT; // T, Ttotal
    double sig2xy0, sig2xy1;
    double distMaxXY0, distMaxXY1, haloMinX, haloMaxX, haloMinY, haloMaxY, xyDistSq, xyterm, bhattaFact;
    int halo1min, halo1max, halo2min, halo2max;

    // FOLLOWING FOR DEBUGGING
    int i1min = nX, i1max = 0, i2min = nY, i2max = 0;

    // start computing MKDE
    cout << "2D MKDE Interaction Computation: STARTING" << endl;
    totalT = 0.0;
    for (int j = 0; j < (nObs - 1); j++) {
        cout << "\tProcessing move step " << (j + 1) << " of " << (nObs - 1) << endl;
    // report percent complete after each observation
        t0 = T[j];
        t1 = T[j + 1];
        dt = t1 - t0;
        t = t0;
        if (isValid[j] == 1) {
            totalT += dt;
            bool exitLoop = false;
            bool finalLoop = false;
            while (!exitLoop) { // iterate over integration time steps
                // Calculate fractional distance between t0 and current time
                alpha = (t - t0) / dt;
                // Calculate parameters for bilinear interpolation
                sig2xy0 = dt * alpha * (1.0 - alpha) * msig2xy0[j] +
                          osig2xy0[j] * (1.0 - alpha) * (1.0 - alpha) +
                          osig2xy0[j + 1] * alpha * alpha;
                sig2xy1 = dt * alpha * (1.0 - alpha) * msig2xy1[j] +
                          osig2xy1[j] * (1.0 - alpha) * (1.0 - alpha) +
                          osig2xy1[j + 1] * alpha * alpha;
                // Get (x,y,z) coordinates of kernel origin using linear interpolation
                eX0 = X0[j] + alpha * (X0[j + 1] - X0[j]);
                eY0 = Y0[j] + alpha * (Y0[j + 1] - Y0[j]);
                eX1 = X1[j] + alpha * (X1[j + 1] - X1[j]);
                eY1 = Y1[j] + alpha * (Y1[j + 1] - Y1[j]);
                // halo
                distMaxXY0 = univariateNormalDensityThreshold(pdfMin[0], sig2xy0); // ADD
                distMaxXY1 = univariateNormalDensityThreshold(pdfMin[0], sig2xy1); // ADD
                // x
                haloMinX = std::min(eX0 - distMaxXY0, eX1 - distMaxXY1); // ADD
                haloMaxX = std::max(eX0 + distMaxXY0, eX1 + distMaxXY1); // ADD
                halo1min = coordToIndex(haloMinX, xGrid[0], xSz); // ADD
                halo1max = coordToIndex(haloMaxX, xGrid[0], xSz); // ADD
                // y
                haloMinY = std::min(eY0 - distMaxXY0, eY1 - distMaxXY1); // ADD
                haloMaxY = std::max(eY0 + distMaxXY0, eY1 + distMaxXY1); // ADD
                halo2min = coordToIndex(haloMinY, yGrid[0], ySz); // ADD
                halo2max = coordToIndex(haloMaxY, yGrid[0], ySz); // ADD

                // Precompute voxel density in y dimension
                for (int i2 = 0; i2 < nY; i2++) {
                    yprob0[i2] = integrateNormal(yGrid[i2] - 0.5 * ySz, yGrid[i2] + 0.5 * ySz, eY0, sqrt(sig2xy0));
                    yprob1[i2] = integrateNormal(yGrid[i2] - 0.5 * ySz, yGrid[i2] + 0.5 * ySz, eY1, sqrt(sig2xy1));
                    ybhattdist[i2] = integrateKernelBC(yGrid[i2] - 0.5 * ySz, yGrid[i2] + 0.5 * ySz, eY0, sqrt(sig2xy0),
                                                       eY1, sqrt(sig2xy1), pdfMin[0]);
                }
                bhattaFact = (RSQRT2PI * RSQRT2PI) / (sqrt(sig2xy0) * sqrt(sig2xy1));

                for (int i1 = std::max(0, halo1min); i1 < std::min(nX, halo1max); i1++) { // x-dimension
                    double voxx = xGrid[i1]; // voxel x
                    double xprob0 = integrateNormal(voxx - 0.5 * xSz, voxx + 0.5 * xSz, eX0, sqrt(sig2xy0));
                    double xprob1 = integrateNormal(voxx - 0.5 * xSz, voxx + 0.5 * xSz, eX1, sqrt(sig2xy1));
                    double xbhattdist = integrateKernelBC(voxx - 0.5 * xSz, voxx + 0.5 * xSz, eX0, sqrt(sig2xy0), eX1,
                                                          sqrt(sig2xy1), pdfMin[0]);
                    for (int i2 = std::max(0, halo2min); i2 < std::min(nY, halo2max); i2++) { // y-dimension
                        double voxy = yGrid[i2]; // voxel y
                        // Calculate contribution of kernel to voxel
                        double pXY0 = xprob0 * yprob0[i2];
                        double pXY1 = xprob1 * yprob1[i2];
                        double pXY = xbhattdist * ybhattdist[i2] * bhattaFact;
                        // update (trapezoid rule)
                        double tmpDens, tmpDens0, tmpDens1;
                        if (doubleEquals(t, t0)) { // first term
                            tmpDens = stepT[0] * pXY / 2.0;
                            tmpDens0 = stepT[0] * pXY0 / 2.0;
                            tmpDens1 = stepT[0] * pXY1 / 2.0;
                        } else if (doubleEquals(t, t1)) { // last term
                            tmpDens = (t - tOld) * pXY / 2.0;
                            tmpDens0 = (t - tOld) * pXY0 / 2.0;
                            tmpDens1 = (t - tOld) * pXY1 / 2.0;
                        } else { // intermediate terms
                            tmpDens = stepT[0] * pXY;
                            tmpDens0 = stepT[0] * pXY0;
                            tmpDens1 = stepT[0] * pXY1;
                        }
                        // Add contribution to voxel (removed Kahan summation for now)
                        long kk = getLinearIndex(i1, i2, 0, nX, nY);
                        mkde->setGridValue((double)i1, (double)i2, mkde->getGridValue((double)i1, (double)i2) + tmpDens);
                        W0 += tmpDens0;
                        W1 += tmpDens1;
                    }
                }
                // update the eval time (t) and stopping conditions
                if (finalLoop) {
                    exitLoop = true;
                } else {
                    tOld = t;
                    t += stepT[0];
                    if (t >= t1) {
                        t = t1;
                        finalLoop = true;
                    }
                }
            }
        }
    }
    // divide by totalT
    double maxDist = 0.0, sumDens = 0.0;
    double normConst = 1.0 / sqrt(W0 * W1); // 1.0 / totalT
    for (double eX = xmin; eX < xmax; eX = eX + xSz) {
        for (double eY = ymin; eY < ymax; eY = eY + ySz) {
            if (mkde->getGridValue(eX, eY) != FLT_NO_DATA_VALUE) {
                mkde->setGridValue(eX, eY, mkde->getGridValue(eX, eY) * normConst);
                if (mkde->getGridValue(eX, eY) > maxDist) {
                    maxDist = mkde->getGridValue(eX, eY);
                }
                sumDens += mkde->getGridValue(eX, eY);
            }
        }
    }
    cout << "\tNormalizing Constant = " << normConst << "(" << 1.0 / normConst << ")" << endl;
    cout << "\tMaximum cell Bhattacharyya coeff = " << maxDist << endl;
    cout << "\tOverall Bhattacharyya coeff = " << sumDens << endl;
    cout << "2D MKDE Interaction Computation: DONE" << endl;
    // RETURN DENSITY HERE
    return mkde;
}


/* this version computes the probability that the two individuals will be in the
   same voxel, not at the exact same location.  To compute the probability of
   occurrence at the same point within a voxel, we will have to integrate the
   product of the two kernels over the volume of the voxel */
/*
gridFloat3D * mkde3dGridv02interact(const vector<double> &T, const vector<double> &X0, const vector<double> &Y0,
                                   const vector<double> &Z0, const vector<double> &X1, const vector<double> &Y1,
                                   const vector<double> &Z1, const vector<bool> &isValid, const vector<double> &xgrid,
                                   const vector<double> &ygrid, const vector<double> &zgrid, gridFloat *zMin,
                                   gridFloat *zMax, const vector<double> &msig2xy0, const vector<double> &msig2xy1,
                                   const vector<double> &msig2z0, const vector<double> &msig2z1, const vector<double> &osig2xy0,
                                   const vector<double> &osig2xy1, const vector<double> &osig2z0, const vector<double> &osig2z1,
                                   const vector<double> &stepT, const vector<double> &pdfMin) {
    // grid specs
    int nX = xgrid.size();
    int nY = ygrid.size();
    int nZ = zgrid.size();
    int nObs = T.size();
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
    double *yprob0 = (double *) malloc(nY * sizeof(double)); // To PRECOMPUTE Y
    double *yprob1 = (double *) malloc(nY * sizeof(double)); // To PRECOMPUTE Y
    double *ybhattdist = (double *) malloc(nY * sizeof(double)); // To PRECOMPUTE Y
    double *zprob0 = (double *) malloc(nZ * sizeof(double)); // To Precompute z
    double *zprob1 = (double *) malloc(nZ * sizeof(double)); // To Precompute z
    double *zbhattdist = (double *) malloc(nZ * sizeof(double)); // To PRECOMPUTE Z

    // Create a vector of GridFloats and initializing GridFloat3D
    gridFloat3D * mkde = new gridFloat3D(xgrid, ygrid, zgrid);

    // set up time variables
    double t0, t1, t, tOld, dt, alpha;

    // set up tmp variables
    double eX0, eY0, eZ0, eX1, eY1, eZ1;
    double W0 = 0.0, W1 = 0.0;
    double sig2xy0, sig2z0, sig2xy1, sig2z1;
    double distMaxXY0, distMaxXY1, distMaxZ0, distMaxZ1, xyDistSq, zDistSq, tmpZ, xyterm, bhattaFact;
    double haloMinX, haloMaxX, haloMinY, haloMaxY, haloMinZ, haloMaxZ;
    double loReflZ1, hiReflZ1;
    int halo1min, halo1max, halo2min, halo2max, halo3min, halo3max; //, i1k, i2k, i3k;

    // FOLLOWING FOR DEBUGGING
    int i1min = nX, i1max = 0, i2min = nY, i2max = 0, i3min = nZ, i3max = 0;

    // start computing MKDE
    cout << "3D MKDE Interaction Computation: STARTING" << endl;
    for (int j = 0; j < (nObs - 1); j++) {
        cout << "\tProcessing move step " << (j + 1) << " of " << (nObs - 1) << std::endl;

        // report percent complete after each observation
        t0 = T[j];
        t1 = T[j + 1];
        dt = t1 - t0;
        t = t0;
        if (isValid[j] == 1) {
            bool exitLoop = false;
            bool finalLoop = false;
            while (!exitLoop) { // iterate over integration time steps

                // Calculate fractional distance between t0 and current time
                alpha = (t - t0) / dt;

                // Calculate parameters for bilinear interpolation
                sig2xy0 = dt * alpha * (1.0 - alpha) * msig2xy0[j] +
                          osig2xy0[j] * (1.0 - alpha) * (1.0 - alpha) +
                          osig2xy0[j + 1] * alpha * alpha;
                sig2z0 = dt * alpha * (1.0 - alpha) * msig2z0[j] +
                         osig2z0[j] * (1.0 - alpha) * (1.0 - alpha) +
                         osig2z0[j + 1] * alpha * alpha;
                sig2xy1 = dt * alpha * (1.0 - alpha) * msig2xy1[j] +
                          osig2xy1[j] * (1.0 - alpha) * (1.0 - alpha) +
                          osig2xy1[j + 1] * alpha * alpha;
                sig2z1 = dt * alpha * (1.0 - alpha) * msig2z1[j] +
                         osig2z1[j] * (1.0 - alpha) * (1.0 - alpha) +
                         osig2z1[j + 1] * alpha * alpha;

                // Get (x,y,z) coordinates of kernel origin using linear interpolation
                eX0 = X0[j] + alpha * (X0[j + 1] - X0[j]);
                eY0 = Y0[j] + alpha * (Y0[j + 1] - Y0[j]);
                eZ0 = Z0[j] + alpha * (Z0[j + 1] - Z0[j]);
                eX1 = X1[j] + alpha * (X1[j + 1] - X1[j]);
                eY1 = Y1[j] + alpha * (Y1[j + 1] - Y1[j]);
                eZ1 = Z1[j] + alpha * (Z1[j + 1] - Z1[j]);

                // halo
                distMaxXY0 = univariateNormalDensityThreshold(pdfMin[0], sig2xy0);
                distMaxXY1 = univariateNormalDensityThreshold(pdfMin[0], sig2xy1);
                distMaxZ0 = univariateNormalDensityThreshold(pdfMin[0], sig2z0);
                distMaxZ1 = univariateNormalDensityThreshold(pdfMin[0], sig2z1);

                // x
                haloMinX = std::min(eX0 - distMaxXY0, eX1 - distMaxXY1);
                haloMaxX = std::max(eX0 + distMaxXY0, eX1 + distMaxXY1);
                halo1min = coordToIndex(haloMinX, xgrid[0], xSz);
                halo1max = coordToIndex(haloMaxX, xgrid[0], xSz);

                // y
                haloMinY = std::min(eY0 - distMaxXY0, eY1 - distMaxXY1);
                haloMaxY = std::max(eY0 + distMaxXY0, eY1 + distMaxXY1);
                halo2min = coordToIndex(haloMinY, ygrid[0], ySz);
                halo2max = coordToIndex(haloMaxY, ygrid[0], ySz);

                // z
                haloMinZ = std::min(eZ0 - distMaxZ0, eZ1 - distMaxZ1);
                haloMaxZ = std::max(eZ0 + distMaxZ0, eZ1 + distMaxZ1);
                halo3min = coordToIndex(haloMinZ, zgrid[0], zSz);
                halo3max = coordToIndex(haloMaxZ, zgrid[0], zSz);

                // Precompute voxel density in y dimension
                for (int i2 = 0; i2 < nY; i2++) {
                    yprob0[i2] = integrateNormal(ygrid[i2] - 0.5 * ySz, ygrid[i2] + 0.5 * ySz, eY0, sqrt(sig2xy0));
                    yprob1[i2] = integrateNormal(ygrid[i2] - 0.5 * ySz, ygrid[i2] + 0.5 * ySz, eY1, sqrt(sig2xy1));
                    ybhattdist[i2] = integrateKernelBC(ygrid[i2] - 0.5 * ySz, ygrid[i2] + 0.5 * ySz, eY0, sqrt(sig2xy0),
                                                       eY1, sqrt(sig2xy1), pdfMin[0]);
                }

                // Precompute voxel density in z dimension
                for (int i3 = 0; i3 < nZ; i3++) {
                    zprob0[i3] = integrateNormal(zgrid[i3] - 0.5 * zSz, zgrid[i3] + 0.5 * zSz, eZ0, sqrt(sig2z0));
                    zprob1[i3] = integrateNormal(zgrid[i3] - 0.5 * zSz, zgrid[i3] + 0.5 * zSz, eZ1, sqrt(sig2z1));
                    zbhattdist[i3] = integrateKernelBC(zgrid[i3] - 0.5 * zSz, zgrid[i3] + 0.5 * zSz, eZ0, sqrt(sig2z0),
                                                       eZ1, sqrt(sig2z1), pdfMin[0]);
                }
                bhattaFact = (RSQRT2PI * RSQRT2PI * RSQRT2PI) /
                             (sqrt(sig2xy0) * sqrt(sig2xy1) * sqrt(sqrt(sig2z0) * sqrt(sig2z1)));

                for (int i1 = std::max(0, halo1min); i1 < std::min(nX, halo1max); i1++) { // x-dimension
                    double voxx = xgrid[i1]; // voxel x
                    double xprob0 = integrateNormal(voxx - 0.5 * xSz, voxx + 0.5 * xSz, eX0, sqrt(sig2xy0));
                    double xprob1 = integrateNormal(voxx - 0.5 * xSz, voxx + 0.5 * xSz, eX1, sqrt(sig2xy1));
                    double xbhattdist = integrateKernelBC(voxx - 0.5 * xSz, voxx + 0.5 * xSz, eX0, sqrt(sig2xy0), eX1,
                                                          sqrt(sig2xy1), pdfMin[0]);
                    for (int i2 = std::max(0, halo2min); i2 < std::min(nY, halo2max); i2++) { // y-dimension
                        double voxy = ygrid[i2]; // voxel y

                        // get the range of indexes and coordinates based on the physical boundaries
                        int i3lo = std::max(0, getLowerCellIndex(zMin(i1, i2), zgrid[0], zSz));
                        int i3hi = std::min(nZ, getUpperCellIndex(zMax(i1, i2), zgrid[0], zSz) +
                                                1); // add 1 because less than
                        double loZ = indexToCellCenterCoord(i3lo, zgrid[0], zSz) - 0.5 * zSz;
                        double hiZ = indexToCellCenterCoord(i3hi - 1, zgrid[0], zSz) + 0.5 * zSz;

                        // Reflect E[Z] about the lower and upper boundaries
                        double loReflZ0 = 2.0 * loZ - eZ0;
                        double hiReflZ0 = 2.0 * hiZ - eZ0;
                        double loReflZ1 = 2.0 * loZ - eZ1;
                        double hiReflZ1 = 2.0 * hiZ - eZ1;

                        for (int i3 = std::max(i3lo, halo3min); i3 < std::min(i3hi, halo3max); i3++) { // z-dimension
                            // set up for reflection: HOW RELEVANT IS THE REFLECTION METHOD FOR INTERACTION?
                            // only compute if the expected location is within distance of boundary
                            double voxz = zgrid[i3];
                            double loDZ0 = 0.0, hiDZ0 = 0.0, loDZ1 = 0.0, hiDZ1 = 0.0, loBD = 0.0, hiBD = 0.0;
                            if ((std::fabs(eZ0 - loZ) <= distMaxZ0) || (std::fabs(eZ1 - loZ) <= distMaxZ1)) {
                                loDZ0 = integrateNormal(voxz - 0.5 * zSz, voxz + 0.5 * zSz, loReflZ0, sqrt(sig2z0));
                                loDZ1 = integrateNormal(voxz - 0.5 * zSz, voxz + 0.5 * zSz, loReflZ1, sqrt(sig2z1));
                                loBD = integrateKernelBC(voxz - 0.5 * zSz, voxz + 0.5 * zSz, loReflZ0, sqrt(sig2z0),
                                                         loReflZ1, sqrt(sig2z1), pdfMin[0]);
                            }
                            if ((std::fabs(hiZ - eZ0) <= distMaxZ0) || (std::fabs(hiZ - eZ1) <= distMaxZ1)) {
                                hiDZ0 = integrateNormal(voxz - 0.5 * zSz, voxz + 0.5 * zSz, hiReflZ0, sqrt(sig2z0));
                                hiDZ1 = integrateNormal(voxz - 0.5 * zSz, voxz + 0.5 * zSz, hiReflZ1, sqrt(sig2z1));
                                hiBD = integrateKernelBC(voxz - 0.5 * zSz, voxz + 0.5 * zSz, hiReflZ0, sqrt(sig2z0),
                                                         hiReflZ1, sqrt(sig2z1), pdfMin[0]);
                            }

                            // Calculate contribution of kernel to voxel
                            double pXYZ0 = xprob0 * yprob0[i2] * (zprob0[i3] + loDZ0 + hiDZ0);
                            double pXYZ1 = xprob1 * yprob1[i2] * (zprob1[i3] + loDZ1 + hiDZ1);
                            double pXYZ = xbhattdist * ybhattdist[i2] * (zbhattdist[i3] + loBD + hiBD) * bhattaFact;

                            // update density (trapezoid rule)
                            double tmpDens, tmpDens0, tmpDens1;
                            if (doubleEquals(t, t0)) { // first term
                                tmpDens0 = stepT[0] * pXYZ0 / 2.0;
                                tmpDens1 = stepT[0] * pXYZ1 / 2.0;
                                tmpDens = stepT[0] * pXYZ / 2.0;
                            } else if (doubleEquals(t, t1)) { // last term
                                tmpDens0 = (t - tOld) * pXYZ0 / 2.0;
                                tmpDens1 = (t - tOld) * pXYZ1 / 2.0;
                                tmpDens = (t - tOld) * pXYZ / 2.0;
                            } else { // intermediate terms
                                tmpDens0 = stepT[0] * pXYZ0;
                                tmpDens1 = stepT[0] * pXYZ1;
                                tmpDens = stepT[0] * pXYZ;
                            }
                            // Add contribution to voxel (removed Kahan summation for now)
                            long kk = getLinearIndex(i1, i2, i3, nX, nY);
                            mkde[kk] += tmpDens; // mkde[kk] + tmpDens;
                            W0 += tmpDens0;
                            W1 += tmpDens1;
                        }
                    }
                }
                // update the eval time (t) and stopping conditions
                if (finalLoop) {
                    exitLoop = true;
                } else {
                    tOld = t;
                    t += stepT[0];
                    if (t >= t1) {
                        t = t1;
                        finalLoop = true;
                    }
                }
            }
        }
    // Rcpp::Rcout << "\t\tW0 = " << W0 << ", W1 = " << W1 << std::endl;
    }

    // divide by totalT
    double maxDist = 0.0, sumDens = 0.0;
    double normConst = 1.0 / sqrt(W0 * W1);
    for (double eZ = zmin; eZ <= zmax; eZ = eZ + zSz) {
        for (double eY = ymin; eY <= ymax; eY = eY + ySz) {
            for (double eX = xmin; eX <= xmax; eX = eX + xSz) {
                if (mkde->getGridValue(eX, eY, eZ) != FLT_NO_DATA_VALUE) {
                    mkde->setGridValue(eX, eY, eZ, normConst * mkde->getGridValue(eX, eY, eZ));
                    if (mkde->getGridValue(eX, eY, eZ) > maxDist) {
                        maxDist = mkde->getGridValue(eX, eY, eZ);
                    }
                    sumDens += mkde->getGridValue(eX, eY, eZ);
                }
            }
        }
    }

    cout << "\tNormalizing Constant = " << normConst << "(" << 1.0 / normConst << ")" << endl;
    cout << "\tMaximum voxel Bhattacharyya coeff = " << maxDist << endl;
    cout << "\tOverall Bhattacharyya coeff = " << sumDens << endl;
    cout << "3D MKDE Interaction Computation: DONE" << endl;

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

/*****************************************************************************
 * Helper functions to check whether the time is within the bounds we're
 * interested in
 *****************************************************************************/
void withinBounds(AnimalData *animal, long minutes) {
    int bound = minutes * 60;
    for (int i = 1; i < animal->t.size(); i++) {
        if (animal->epoch_seconds[i] - animal->epoch_seconds[i - 1] >= bound) {
            animal->use[i] = false;
            animal->use[i - 1] = false;
        }
    }
}

/*****************************************************************************
 * Helper functions to adjust time relative to the first observation time.
 *****************************************************************************/
void updateTime(AnimalData *animal) {
    for (int i = 0; i < animal->epoch_seconds.size(); i++) {
        animal->epoch_seconds[i] = animal->epoch_seconds[i] - animal->epoch_seconds[0];
    }
}

/*****************************************************************************
 * Helper functions for spatial-temporal interaction
 *****************************************************************************/

/* Evaluate the square root of the product of the two kernels at x
 * We will very likely want to optimize the hell out of this function
 */
double kernelBC(double x, double mu1, double sigma1sq, double mu2, double sigma2sq) {
    double res = exp(-0.25*(x - mu1)*(x - mu1)/sigma1sq - 0.25*(x - mu2)*(x - mu2)/sigma2sq);
    return(res);
}

/* Based on trapzd() from Press et al.
 * Meant to be used only in a function like integrateKernelBC
 */
double trapzdKernelBC(double x0, double x1, double mu1, double sigma1sq, double mu2, double sigma2sq, int n, double old_s) {
    double x, tnm, sum, del;
    double s; // this was static in Press et al, but I made it so you have to pass the old value as an arg
    int it, j;
    if (n == 1) {
        s = 0.5*(x1 - x0)*(kernelBC(x1, mu1, sigma1sq, mu2, sigma2sq) + kernelBC(x0, mu1, sigma1sq, mu2, sigma2sq));
    } else {
        for (it = 1, j = 1; j < n - 1; j++) {
            it <<= 1;
        }
        tnm = it;
        del = (x1 - x0)/tnm;
        x = x0 + 0.5*del;
        for (sum = 0.0, j = 1; j <= it; j++, x += del) {
            sum += kernelBC(x, mu1, sigma1sq, mu2, sigma2sq);
        }
        s = 0.5*(old_s + (x1 - x0)*sum/tnm); // updates s if n > 1
    }
    return s;
}

/* Based on qsimp() from Press et al. */
double integrateKernelBC(double x0, double x1, double mu1, double sigma1, double mu2, double sigma2, double pThresh) {
    double sigma1sq = sigma1*sigma1, sigma2sq = sigma2*sigma2;
    double s, st, ost = 0.0, os = 0.0; // is this a local variable s?
    for (int j = 1; j <= JMAX; j++) {
        st = trapzdKernelBC(x0, x1, mu1, sigma1sq, mu2, sigma2sq, j, ost);
        s = (4.0*st - ost)/3.0;
        if (j > 5) {
            if (fabs(s - os) < INT_EPS*fabs(os) || (s == 0.0 && os == 0.0) || (st < pThresh)) {
                return s;
            } // also bail if st is less than minimum threshold so don't waste time on negligible associations
        }
        os = s;
        ost = st;
    }
    return 0.0;
}
/*-----------------------------------------------------------------------------
* Test byte order of machine                                                  *
-----------------------------------------------------------------------------*/
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