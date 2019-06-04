/* This file contains the main logic for MKDE. It processes the user for input and
 * runs the corresponding interaction and prints the results.
 */

#include "mkde.h"
#include "Raster/ioFunctions.h"
#define padding 1000    //padding for the grids

int main() {

    //Prompts user for file with data and processes it
    unordered_map<string, AnimalData *> *animals;
    cout << "File path to animal data: ";
 	string filepath;
    cin >> filepath;
    animals = fileRead(filepath);
    if (animals == nullptr) {   //Reading the file has failed
    	cout << "Make sure the path is absolute!" << endl;
    	return 0;
    }

    //Prompts user for interaction type
    cout << "\n(1) 2-Dimension" << "\n(2) 3-Dimension"
    			<< "\n(3) 2-Dimension Interaction" 
    			<< "\n(4) 3-Dimension Interaction" << endl << "\nType of interaction: ";
    short interaction;
    cin >> interaction;


    // Prompt user for output file type
    short outfile;
    if (interaction == 1 || interaction == 3) { //2D analysis
        cout << "\n(1) ESRI ascii" <<  "\n(2) ESRI binary" <<"\nType of output file: ";
    }
    else if (interaction == 2 || interaction == 4) { //3D analysis
        cout << "\n(1) VTK" << "\n(2) GRASS" << "\n(3) XDMF" << endl << "\nType of output file: ";
    }
    cin >> outfile;

    //Set up the grids for calculation
    vector<double> grid_x;  //x grid
    vector<double> grid_y;  //y grid
    vector<double> grid_z;  //z grid
    vector<time_t> grid_t;  //time grid
    vector<double> move_var;    //vector of move variance
    vector<double> obs_var;     //vector of observation variance
    gridFloat * rst;        //raster for 2d interaction
    gridFloat3D *rst3d;     //raster for 3f interaction
    double t_step = 1.0;    //size of time step
    double pdf_thresh = pow(10.0, -14); //threshold
    vector<double> t_step3d;    //vector of time step

    //Find the lowest and highest x, y, z for all animals
    double xmin = numeric_limits<double>::max();
    double xmax = 0;
    double ymin = numeric_limits<double>::max();
    double ymax = 0;
    double zmin = numeric_limits<double>::max();
    double zmax = 0;
    double tmin = numeric_limits<time_t>::max();
    double tmax = 0;
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
        if (it->second->tmin < tmin) {
            tmin = it->second->tmin;
        }
        if (it->second->tmax > tmax) {
            tmax = it->second->tmax;
        }
    }

    //Set the x, y, z, t grids with the minimum and maximum
    for (double i = xmin - padding; i < xmax + padding; i = i + 500) {
        grid_x.push_back(i);
    }
    for (double i = ymin - padding; i < ymax + padding; i = i + 500) {
        grid_y.push_back(i);
    }
    for (double i = zmin - padding; i < zmax + padding; i = i + 500) {
        grid_z.push_back(i);
    }
    for (time_t i = tmin - padding; i < tmax + padding; i = i + 2000) {
        grid_t.push_back(i);
    }

    //Set up the move variance, observation variance, and t step vectors
    for (int i = 0; i < animals->begin()->second->x.size(); i++) {
        move_var.push_back(2.0);
        obs_var.push_back(0.1);
        t_step3d.push_back(t_step);
    }

    //Set up zmin and zmax grids
    gridFloat *high = new gridFloat(1, grid_y.size(), grid_x.size(), xmin, ymin, grid_x[1] - grid_x[0]);
    high->setAllCellsTo(zmax);
    gridFloat *low = new gridFloat(0, grid_y.size(), grid_x.size(), xmin, ymin, grid_x[1] - grid_x[0]);
    low->setAllCellsToZero(true);

    //Regular mkde 2D and mkde 3D
    if (interaction == 1 || interaction == 2) {
    	for (auto it = animals->begin(); it != animals->end(); ++it) {

        	//2D for individual animals
        	if (interaction == 1) {

        	    //Perform actual calculations
        		rst = mkde2D(it->second->t, it->second->x, it->second->y, it->second->use,
            	grid_x, grid_y, it->second->moveVarXY, it->second->obsVarXY, t_step, pdf_thresh);

        		//Prints results
                if (outfile == 1) { //ESRI ascii
                    string filename = it->first + ".asc";
                    rst->printESRIascii(filename);
                }
                if (outfile == 2) { //ESRI binary
                    string filename = it->first + ".hdr";
                    char * cstr = &filename[0u];
                    rst->printESRIbinary(cstr);
                }

         /*       printPoints(rst, "points.txt");
                vector<double> tempZ;
                tempZ.push_back(0);
                tempZ.push_back(1);

                rst3d = new gridFloat3D(rst, grid_x, grid_y);
                writeMKDE3DtoVTK(grid_x, grid_y, tempZ, rst3d, filename, "test"); */
        	}

        	// 3D for individual animals
        	if (interaction == 2) {
        	    if (it->second->zmin == 0 && it->second->zmax == 1) {   //Data did not include z column
        	        cout << "Failed to do 3D analysis without Z column" << endl;
        	        return 1;
        	    }

                //Perform actual calculations
        		rst3d = mkde3dGridv02(it->second->t, it->second->x, it->second->y, it->second->z, it->second->use,
                              grid_x, grid_y, grid_z, low, high, it->second->moveVarXY,
                              it->second->moveVarZ, it->second->obsVarXY, it->second->obsVarZ, t_step3d,
                              pdf_thresh);

                //Prints results
                if (outfile == 1) { //VTK
                    string filename = it->first + ".vtk";
                    writeMKDE3DtoVTK(grid_x, grid_y, grid_z, rst3d, filename, "results from mkde3d");
                }
                else if (outfile == 2) { //GRASS
                    string filename = it->first + ".grass";
                    writeMKDE3DtoGRASS(grid_x, grid_y, grid_z, rst3d, filename, "results from mkde3d");
                }
                else if (outfile == 3) { //XDMF
                    string filename = it->first + ".xdmf";
                    writeMKDE3DtoXDMF(grid_x, grid_y, grid_z, rst3d, filename, "results from mkde3d");
                }
        	}
    	}
	}

    //Interaction mkde 2D and mkde 3D
	else if (interaction == 3 || interaction == 4) {
    	vector<AnimalData *> avec;	//Create a vector of the animals and do interaction on pairs
    	for (auto it = animals->begin(); it != animals->end(); ++it) {
        	avec.push_back(it->second);
    	}

    	vector<double> indexes0 = assignLocationIndexToTimeGrid(grid_t, avec[0]->epochSeconds, 1000000);
    	vector<double> indexes1 = assignLocationIndexToTimeGrid(grid_t, avec[1]->epochSeconds, 1000000);
    
    	//2D interaction
    	if (interaction == 3) {
            vector<pointIn3D> alpha_0 = interpolateCoordinateOnTimeGrid2d(grid_t, indexes0, avec[0]->epochSeconds,
                    avec[0]->xyz, avec[0]->moveVarXY, avec[0]->obsVarXY);
            vector<pointIn3D> alpha_1 = interpolateCoordinateOnTimeGrid2d(grid_t, indexes1, avec[1]->epochSeconds,
                    avec[1]->xyz, avec[1]->moveVarXY, avec[1]->obsVarXY); 
    		rst = mkde2dGridv02interact(avec[0]->t, avec[0]->x, avec[0]->y, avec[1]->x, avec[1]->y, avec[0]->use,
                    grid_x, grid_y, alpha_0, alpha_1, pdf_thresh);

            //Prints results
            if (outfile == 1) { //ESRI ascii
                string filename = "interaction2d.asc";
                rst->printESRIascii(filename);
            }
            if (outfile == 2) { //ESRI binary
                string filename = "interaction2d.hdr";
                char * cstr = &filename[0u];
                rst->printESRIbinary(cstr);
            }
		}

		//3D interaction
		if (interaction == 4) {
		    if (avec[0]->zmin == 0 && avec[0]->zmax == 1) {     //Data did not include 3D column
                cout << "Failed to do 3D analysis without Z column" << endl;
                return 1;

		    }
            vector<pointIn3D> alpha_0 = interpolateCoordinateOnTimeGrid(grid_t, indexes0, avec[0]->epochSeconds,
                    avec[0]->xyz, avec[0]->moveVarXY, avec[0]->moveVarZ, avec[0]->obsVarXY, avec[0]->obsVarZ);
            vector<pointIn3D> alpha_1 = interpolateCoordinateOnTimeGrid(grid_t, indexes1, avec[1]->epochSeconds,
                    avec[1]->xyz, avec[1]->moveVarXY, avec[1]->moveVarZ, avec[1]->obsVarXY, avec[1]->obsVarZ); 
    		rst3d = mkde3dGridv02interact(avec[0]->t, avec[0]->x, avec[0]->y, avec[0]->z, avec[1]->x, avec[1]->y,
    		        avec[1]->z, avec[0]->use, grid_x, grid_y, grid_z, low, high, alpha_0, alpha_1, pdf_thresh);

            //Prints results
            if (outfile == 1) { //VTK
                string filename = "interaction3d.vtk";
                writeMKDE3DtoVTK(grid_x, grid_y, grid_z, rst3d, filename, "results from mkde3d");
            }
            else if (outfile == 2) { //GRASS
                string filename = "interaction3d.grass";
                writeMKDE3DtoGRASS(grid_x, grid_y, grid_z, rst3d, filename, "results from mkde3d");
            }
            else if (outfile == 3) { //XDMF
                string filename = "interaction3d.xdmf";
                writeMKDE3DtoXDMF(grid_x, grid_y, grid_z, rst3d, filename, "results from mkde3d");
            }
    	}
    }
    return 1;
}


/*****************************************************************************
 * Functions that calculate MKDEs for a set of grid cells
 *****************************************************************************/

/*
 * Function used to compute 2D case of MKDE. 
 * T: vector of observed time
 * X: vector of observed x 
 * Y: vector of observed y
 * use: whether we should use the corresponding point
 * grid_x: vector of x grid
 * grid_y: vector of y grid
 * move_var: move variance
 * obs_var: observation variance
 * t_step: time step to integrate over
 * pdf_thresh: threshold
 */
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

    gridFloat *rst = new gridFloat(num_grids++, nY, nX, xmin, ymin, xSz);
    rst->setAllCellsToZero(true);

    // set up time variables
    double t0, t1, t, tOld, dt, alpha;
    float totalT;

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
    for (double eY = ymin; eY <= ymax; eY = eY + ySz) {
        for (double eX = xmin; eX <= xmax; eX = eX + xSz) {
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

/*
 * Function used to compute 2D case of MKDE. 
 * T: vector of observed time
 * X: vector of observed x 
 * Y: vector of observed y
 * Z: vector of observed z
 * use: whether we should use the corresponding point
 * grid_x: vector of x grid
 * grid_y: vector of y grid
 * grid_z: vector of z grid
 * zMin: grid of the minimum value for Z
 * zMax: grid of the maxmimum value of Z
 * msig2xy: move variance for xy
 * msigz: move variance for z
 * osig2xy: observation variance for xy
 * osig2z: observation variance for z
 * t_step: time step to integrate over
 * pdf_thresh: threshold
 */
gridFloat3D * mkde3dGridv02(const vector<double> &T, const vector<double> &X,
                            const vector<double> &Y, const vector<double> &Z, const vector<bool> &use,
                            const vector<double> &xgrid, const vector<double> &ygrid, const vector<double> &zgrid,
                            gridFloat *zMin, gridFloat *zMax, const vector<double> &msig2xy,
                            const vector<double> &msig2z, const vector<double> &osig2xy, const vector<double> &osig2z,
                            const vector<double> &t_step, const double &pdf_thresh) {

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

                distMaxXY = univariateNormalDensityThreshold(pdf_thresh, sig2xy);
                halo1 = ceil(distMaxXY / xSz);
                halo2 = ceil(distMaxXY / ySz);
                distMaxZ = univariateNormalDensityThreshold(pdf_thresh, sig2z);
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
                                  const vector<pointIn3D> &alpha_0, const vector<pointIn3D> &alpha_1,
                                  const double &pdfMin) {

    int nObs = T.size();
    static int num_grids = 0;
    int stepT = T[1] - T[0];

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

    gridFloat * mkde = new gridFloat(num_grids++, nY, nX, xmin, ymin, xSz);
    mkde->setAllCellsToZero(true);

    // set up time variables
    double t0, t1, t, tOld, dt, alpha;
    // set up tmp variables
    double W0 = 0.0, W1 = 0.0;
    double totalT; // T, Ttotal
    double distMaxXY0, distMaxXY1, haloMinX, haloMaxX, haloMinY, haloMaxY, xyDistSq, xyterm, bhattaFact;
    int halo1min, halo1max, halo2min, halo2max;

    // FOLLOWING FOR DEBUGGING
    int i1min = nX, i1max = 0, i2min = nY, i2max = 0;

    // start computing MKDE
    cout << "2D MKDE Interaction Computation: STARTING" << endl;
    totalT = 0.0;
    for (int j = 0; j < T.size(); j++) {
        cout << "\tProcessing move step " << (j + 1) << " of " << T.size() << endl;
        if (isValid[1] == 1) {
            // halo
            distMaxXY0 = univariateNormalDensityThreshold(pdfMin, alpha_0[j].sig2xy); // ADD
            distMaxXY1 = univariateNormalDensityThreshold(pdfMin, alpha_1[j].sig2xy); // ADD

            // x
            haloMinX = std::min(alpha_0[j].x - distMaxXY0, alpha_1[j].x - distMaxXY1); // ADD
            haloMaxX = std::max(alpha_0[j].x + distMaxXY0, alpha_1[j].x + distMaxXY1); // ADD
            halo1min = coordToIndex(haloMinX, xGrid[0], xSz); // ADD
            halo1max = coordToIndex(haloMaxX, xGrid[0], xSz); // ADD

            // y
            haloMinY = std::min(alpha_0[j].y - distMaxXY0, alpha_1[j].y - distMaxXY1); // ADD
            haloMaxY = std::max(alpha_0[j].y + distMaxXY0, alpha_1[j].y + distMaxXY1); // ADD
            halo2min = coordToIndex(haloMinY, yGrid[0], ySz); // ADD
            halo2max = coordToIndex(haloMaxY, yGrid[0], ySz); // ADD

            // Precompute voxel density in y dimension
            for (int i2 = 0; i2 < nY; i2++) {
                yprob0[i2] = integrateNormal(yGrid[i2] - 0.5 * ySz, yGrid[i2] + 0.5 * ySz, alpha_0[j].y,
                                             sqrt(alpha_0[j].sig2xy));
                yprob1[i2] = integrateNormal(yGrid[i2] - 0.5 * ySz, yGrid[i2] + 0.5 * ySz, alpha_1[j].y,
                                             sqrt(alpha_1[j].sig2xy));
                ybhattdist[i2] = integrateKernelBC(yGrid[i2] - 0.5 * ySz, yGrid[i2] + 0.5 * ySz, alpha_0[j].y,
                                                   sqrt(alpha_0[j].sig2xy),
                                                   alpha_1[j].y, sqrt(alpha_1[j].sig2xy), pdfMin);
            }
            bhattaFact = (RSQRT2PI * RSQRT2PI) / (sqrt(alpha_0[j].sig2xy) * sqrt(alpha_1[j].sig2xy));

            for (int i1 = std::max(0, halo1min); i1 < std::min(nX, halo1max); i1++) { // x-dimension
                double voxx = xGrid[i1]; // voxel x
                double xprob0 = integrateNormal(voxx - 0.5 * xSz, voxx + 0.5 * xSz, alpha_0[j].x,
                                                sqrt(alpha_0[j].sig2xy));
                double xprob1 = integrateNormal(voxx - 0.5 * xSz, voxx + 0.5 * xSz, alpha_1[j].x,
                                                sqrt(alpha_1[j].sig2xy));
                double xbhattdist = integrateKernelBC(voxx - 0.5 * xSz, voxx + 0.5 * xSz, alpha_0[j].x,
                                                      sqrt(alpha_0[j].sig2xy), alpha_1[j].x, sqrt(alpha_1[j].sig2xy),
                                                      pdfMin);
                for (int i2 = std::max(0, halo2min); i2 < std::min(nY, halo2max); i2++) { // y-dimension
                    double voxy = yGrid[i2]; // voxel y
                    // Calculate contribution of kernel to voxel
                    double pXY0 = xprob0 * yprob0[i2];
                    double pXY1 = xprob1 * yprob1[i2];
                    double pXY = xbhattdist * ybhattdist[i2] * bhattaFact;

                    // update (trapezoid rule)
                    double tmpDens, tmpDens0, tmpDens1;
                    if (doubleEquals(t, t0)) { // first term
                        tmpDens = stepT * pXY / 2.0;
                        tmpDens0 = stepT * pXY0 / 2.0;
                        tmpDens1 = stepT * pXY1 / 2.0;
                    } else if (doubleEquals(t, t1)) { // last term
                        tmpDens = (t - tOld) * pXY / 2.0;
                        tmpDens0 = (t - tOld) * pXY0 / 2.0;
                        tmpDens1 = (t - tOld) * pXY1 / 2.0;
                    } else { // intermediate terms
                        tmpDens = stepT * pXY;
                        tmpDens0 = stepT * pXY0;
                        tmpDens1 = stepT * pXY1;
                    }

                    // Add contribution to voxel (removed Kahan summation for now)
                    double eX = indexToCellCenterCoord(i1, xmin, xSz);
                    double eY = indexToCellCenterCoord(i2, ymin, ySz);
                    mkde->addValueToGridCell(eX, eY,tmpDens);
                    W0 += tmpDens0;
                    W1 += tmpDens1;  
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
gridFloat3D * mkde3dGridv02interact(const vector<double> &T, const vector<double> &X0, const vector<double> &Y0,
                                    const vector<double> &Z0, const vector<double> &X1, const vector<double> &Y1,
                                    const vector<double> &Z1, const vector<bool> &isValid, const vector<double> &xGrid,
                                    const vector<double> &yGrid, const vector<double> &zGrid, gridFloat *zMin,
                                    gridFloat *zMax, const vector<pointIn3D> &alpha_0, const vector<pointIn3D> &alpha_1,
                                    const double &pdfMin) {
    // grid specs
    int stepT = T[1] - T[0];
    int nX = xGrid.size();
    int nY = yGrid.size();
    int nZ = zGrid.size();
    int nObs = T.size();
    long nVoxels = (long) nX * nY * nZ;
    double xSz = xGrid[1] - xGrid[0];
    double ySz = yGrid[1] - yGrid[0];
    double zSz = zGrid[1] - zGrid[0];
    double xmin = xGrid[0];
    double xmax = xGrid[xGrid.size() - 1];
    double ymin = yGrid[0];
    double ymax = yGrid[yGrid.size() - 1];
    double zmin = zGrid[0];
    double zmax = zGrid[zGrid.size() - 1];

    // arrays for MKDE computation
    double *yprob0 = (double *) malloc(nY * sizeof(double)); // To PRECOMPUTE Y
    double *yprob1 = (double *) malloc(nY * sizeof(double)); // To PRECOMPUTE Y
    double *ybhattdist = (double *) malloc(nY * sizeof(double)); // To PRECOMPUTE Y
    double *zprob0 = (double *) malloc(nZ * sizeof(double)); // To Precompute z
    double *zprob1 = (double *) malloc(nZ * sizeof(double)); // To Precompute z
    double *zbhattdist = (double *) malloc(nZ * sizeof(double)); // To PRECOMPUTE Z

    // Create a vector of GridFloats and initializing GridFloat3D
    gridFloat3D * mkde = new gridFloat3D(xGrid, yGrid, zGrid);

    // set up time variables
    double t0, t1, t, tOld, dt, alpha;

    // set up tmp variables
    double W0 = 0.0, W1 = 0.0;
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

        // halo
        distMaxXY0 = univariateNormalDensityThreshold(pdfMin, alpha_0[j].sig2xy);
        distMaxXY1 = univariateNormalDensityThreshold(pdfMin, alpha_1[j].sig2xy);
        distMaxZ0 = univariateNormalDensityThreshold(pdfMin, alpha_0[j].sig2z);
        distMaxZ1 = univariateNormalDensityThreshold(pdfMin, alpha_1[j].sig2z);

        // x
        haloMinX = std::min(alpha_0[j].x - distMaxXY0, alpha_1[j].x - distMaxXY1);
        haloMaxX = std::max(alpha_0[j].x + distMaxXY0, alpha_1[j].x + distMaxXY1);
        halo1min = coordToIndex(haloMinX, xGrid[0], xSz);
        halo1max = coordToIndex(haloMaxX, xGrid[0], xSz);

        // y
        haloMinY = std::min(alpha_0[j].y - distMaxXY0, alpha_1[j].y - distMaxXY1);
        haloMaxY = std::max(alpha_0[j].y + distMaxXY0, alpha_1[j].y + distMaxXY1);
        halo2min = coordToIndex(haloMinY, yGrid[0], ySz);
        halo2max = coordToIndex(haloMaxY, yGrid[0], ySz);

        // z
        haloMinZ = std::min(alpha_0[j].z - distMaxZ0, alpha_1[j].z - distMaxZ1);
        haloMaxZ = std::max(alpha_0[j].z + distMaxZ0, alpha_1[j].z + distMaxZ1);
        halo3min = coordToIndex(haloMinZ, zGrid[0], zSz);
        halo3max = coordToIndex(haloMaxZ, zGrid[0], zSz);

        // Precompute voxel density in y dimension
        for (int i2 = 0; i2 < nY; i2++) {
            yprob0[i2] = integrateNormal(yGrid[i2] - 0.5 * ySz, yGrid[i2] + 0.5 * ySz, alpha_0[j].y, sqrt(alpha_0[j].sig2xy));
            yprob1[i2] = integrateNormal(yGrid[i2] - 0.5 * ySz, yGrid[i2] + 0.5 * ySz, alpha_1[j].y, sqrt(alpha_1[j].sig2xy));
            ybhattdist[i2] = integrateKernelBC(yGrid[i2] - 0.5 * ySz, yGrid[i2] + 0.5 * ySz, alpha_0[j].y, sqrt(alpha_0[j].sig2xy),
                                               alpha_1[j].y, sqrt(alpha_1[j].sig2xy), pdfMin);
        }

        // Precompute voxel density in z dimension
        for (int i3 = 0; i3 < nZ; i3++) {
            zprob0[i3] = integrateNormal(zGrid[i3] - 0.5 * zSz, zGrid[i3] + 0.5 * zSz, alpha_0[j].z, sqrt(alpha_0[j].sig2z));
            zprob1[i3] = integrateNormal(zGrid[i3] - 0.5 * zSz, zGrid[i3] + 0.5 * zSz, alpha_1[j].z, sqrt(alpha_1[j].sig2z));
            zbhattdist[i3] = integrateKernelBC(zGrid[i3] - 0.5 * zSz, zGrid[i3] + 0.5 * zSz, alpha_0[j].z, sqrt(alpha_0[j].sig2z),
                                               alpha_1[j].z, sqrt(alpha_1[j].sig2z), pdfMin);
        }
        bhattaFact = (RSQRT2PI * RSQRT2PI * RSQRT2PI) /
                     (sqrt(alpha_0[j].sig2xy) * sqrt(alpha_1[j].sig2xy) * sqrt(sqrt(alpha_0[j].sig2z) * sqrt(alpha_1[j].sig2z)));

        for (int i1 = std::max(0, halo1min); i1 < std::min(nX, halo1max); i1++) { // x-dimension
            double voxx = xGrid[i1]; // voxel x
            double xprob0 = integrateNormal(voxx - 0.5 * xSz, voxx + 0.5 * xSz, alpha_0[j].x, sqrt(alpha_0[j].sig2xy));
            double xprob1 = integrateNormal(voxx - 0.5 * xSz, voxx + 0.5 * xSz, alpha_1[j].x, sqrt(alpha_1[j].sig2xy));
            double xbhattdist = integrateKernelBC(voxx - 0.5 * xSz, voxx + 0.5 * xSz, alpha_0[j].x, sqrt(alpha_0[j].sig2xy),
                                                  alpha_1[j].x, sqrt(alpha_1[j].sig2xy), pdfMin);
            for (int i2 = std::max(0, halo2min); i2 < std::min(nY, halo2max); i2++) { // y-dimension
                double voxy = yGrid[i2]; // voxel y

                // get the range of indexes and coordinates based on the physical boundaries
                int i3lo = std::max(0, getLowerCellIndex(/*zMin(i1, i2)*/0, zGrid[0], zSz));
                int i3hi = std::min(nZ, getUpperCellIndex(/*zMax(i1, i2)*/5000, zGrid[0], zSz) +
                                        1); // add 1 because less than
                double loZ = indexToCellCenterCoord(i3lo, zGrid[0], zSz) - 0.5 * zSz;
                double hiZ = indexToCellCenterCoord(i3hi - 1, zGrid[0], zSz) + 0.5 * zSz;

                // Reflect E[Z] about the lower and upper boundaries
                double loReflZ0 = 2.0 * loZ - alpha_0[j].z;
                double hiReflZ0 = 2.0 * hiZ - alpha_0[j].z;
                double loReflZ1 = 2.0 * loZ - alpha_1[j].z;
                double hiReflZ1 = 2.0 * hiZ - alpha_1[j].z;

                for (int i3 = std::max(i3lo, halo3min); i3 < std::min(i3hi, halo3max); i3++) { // z-dimension
                    // set up for reflection: HOW RELEVANT IS THE REFLECTION METHOD FOR INTERACTION?
                    // only compute if the expected location is within distance of boundary
                    double voxz = zGrid[i3];
                    double loDZ0 = 0.0, hiDZ0 = 0.0, loDZ1 = 0.0, hiDZ1 = 0.0, loBD = 0.0, hiBD = 0.0;
                    if ((std::fabs(alpha_0[j].z - loZ) <= distMaxZ0) || (std::fabs(alpha_1[j].z - loZ) <= distMaxZ1)) {
                        loDZ0 = integrateNormal(voxz - 0.5 * zSz, voxz + 0.5 * zSz, loReflZ0, sqrt(alpha_0[j].sig2z));
                        loDZ1 = integrateNormal(voxz - 0.5 * zSz, voxz + 0.5 * zSz, loReflZ1, sqrt(alpha_1[j].sig2z));
                        loBD = integrateKernelBC(voxz - 0.5 * zSz, voxz + 0.5 * zSz, loReflZ0, sqrt(alpha_0[j].sig2z),
                                                 loReflZ1, sqrt(alpha_1[j].sig2z), pdfMin);
                    }
                    if ((std::fabs(hiZ - alpha_0[j].z) <= distMaxZ0) || (std::fabs(hiZ - alpha_1[j].z) <= distMaxZ1)) {
                        hiDZ0 = integrateNormal(voxz - 0.5 * zSz, voxz + 0.5 * zSz, hiReflZ0, sqrt(alpha_0[j].sig2z));
                        hiDZ1 = integrateNormal(voxz - 0.5 * zSz, voxz + 0.5 * zSz, hiReflZ1, sqrt(alpha_1[j].sig2z));
                        hiBD = integrateKernelBC(voxz - 0.5 * zSz, voxz + 0.5 * zSz, hiReflZ0, sqrt(alpha_0[j].sig2z),
                                                 hiReflZ1, sqrt(alpha_1[j].sig2z), pdfMin);
                    }

                    // Calculate contribution of kernel to voxel
                    double pXYZ0 = xprob0 * yprob0[i2] * (zprob0[i3] + loDZ0 + hiDZ0);
                    double pXYZ1 = xprob1 * yprob1[i2] * (zprob1[i3] + loDZ1 + hiDZ1);
                    double pXYZ = xbhattdist * ybhattdist[i2] * (zbhattdist[i3] + loBD + hiBD) * bhattaFact;

                    // update density (trapezoid rule)
                    double tmpDens, tmpDens0, tmpDens1;
                    if (doubleEquals(t, t0)) { // first term
                        tmpDens0 = stepT * pXYZ0 / 2.0;
                        tmpDens1 = stepT * pXYZ1 / 2.0;
                        tmpDens = stepT * pXYZ / 2.0;
                    } else if (doubleEquals(t, t1)) { // last term
                        tmpDens0 = (t - tOld) * pXYZ0 / 2.0;
                        tmpDens1 = (t - tOld) * pXYZ1 / 2.0;
                        tmpDens = (t - tOld) * pXYZ / 2.0;
                    } else { // intermediate terms
                        tmpDens0 = stepT * pXYZ0;
                        tmpDens1 = stepT * pXYZ1;
                        tmpDens = stepT * pXYZ;
                    }
                    // Add contribution to voxel (removed Kahan summation for now)
                    double eX = indexToCellCenterCoord(i1, xmin, xSz);
                    double eY = indexToCellCenterCoord(i2, ymin, ySz);
                    double eZ = indexToCellCenterCoord(i3, zmin, zSz);
                    mkde->addValueToGridCell(eX, eY, eZ, tmpDens); // mkde[kk] + tmpDens;
                    W0 += tmpDens0;
                    W1 += tmpDens1;
                }
            }
        }
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


// Prints all the cells with nonzero densities
void printDensities(gridFloat3D * mkde, string filename) {
    std::ofstream file;
    file.open(filename);    // open the file stream
    file << "x,y,z,density" << endl;
    for (double eZ = mkde->zmin; eZ <= mkde->zmax; eZ = eZ + mkde->zsize) {
        for (double eY = mkde->ymin; eY <= mkde->ymax; eY = eY + mkde->ysize) {
            for (double eX = mkde->xmin; eX <= mkde->xmax; eX = eX + mkde->xsize) {
                if (mkde->getGridValue(eX, eY, eZ) > 0) {
                    file << eX << "," << eY << "," << eZ << "," << mkde->getGridValue(eX, eY, eZ) << endl;
                }
            }
        }
    }
    file.close();
}

