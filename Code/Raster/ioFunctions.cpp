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
        double zmin = numeric_limits<double>::max();
        double xmax = 0;
        double ymax = 0;
        double zmax = 0;

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
            if (it->second->z[i] < zmin) {
                zmin =it->second->z[i];
            }
            if (it->second->z[i] > zmax) {
                zmax = it->second->z[i];
            }
        }
        it->second->xmin = xmin;
        it->second->xmax = xmax;
        it->second->ymin = ymin;
        it->second->ymax = ymax;
        it->second->zmin = zmin;
        it->second->zmax = zmax;
    }

    return animals;
}



/*****************************************************************************
 * Writes MKDE3D to a VTK file format with the name of the file being fname
 *****************************************************************************/
void writeMKDE3DtoVTK(const vector<double> &xgrid, const vector<double> &ygrid, const vector<double> &zgrid,
                      gridFloat3D * rst3d, string fname, string descr) {
    int nX = xgrid.size();
    int nY = ygrid.size();
    int nZ = zgrid.size();

    long nPoints = (long)nX*nY*nZ, ijk = 0;
/*    char * fnm = new char[strlen(fname)+1];
    fnm[strlen(fname)] = 0;
    memcpy(fnm, fname, strlen(fname)); */

    double densTmp;

    std::ofstream vtkfile;
    vtkfile.open (fname);
    // write header
    vtkfile << "# vtk DataFile Version 3.0" << std::endl;
    vtkfile << descr << std::endl;
    vtkfile << "ASCII" << std::endl;
    vtkfile << "DATASET STRUCTURED_GRID" << std::endl;
    vtkfile << "DIMENSIONS " << nX << " " << nY << " " << nZ << std::endl;
    // write points
    vtkfile << "POINTS " << nPoints << " float" << std::endl;
    vtkfile << std::scientific;
    for (int k = 0; k < nZ; k++) {
        for (int j = 0; j < nY; j++) {
            for (int i = 0; i < nX; i++) {
                vtkfile << xgrid[i] << " " << ygrid[j] << " " << zgrid[k] << std::endl;
            }
        }
    }
    // write data
    vtkfile << std::endl << "POINT_DATA " << nPoints << std::endl;
    vtkfile << "SCALARS density float 1" << std::endl;
    vtkfile << "LOOKUP_TABLE densityTable" << std::endl;
    for (double eZ = rst3d->zmin; eZ < rst3d->zmax; eZ = eZ + rst3d->zsize) {
        for (double eY = rst3d->ymin; eY < rst3d->ymax; eY = eY + rst3d->ysize) {
            for (long eX = rst3d->xmin; eX < rst3d->xmax; eX = eX + rst3d->xsize) {
                if (rst3d->getGridValue(eX, eY, eZ) == FLT_NO_DATA_VALUE) {
                    vtkfile << "0.0000000000000" << std::endl;
                } else {
                        vtkfile << rst3d->getGridValue(eX, eY, eZ) << std::endl;
                }
            }
        }
    }

    // write lookup table
    vtkfile << std::endl << "LOOKUP_TABLE densityTable 9" << std::endl;
    vtkfile << "255 255 204 0.3" << std::endl;
    vtkfile << "255 237 160 0.4" << std::endl;
    vtkfile << "254 217 118 0.5" << std::endl;
    vtkfile << "254 178 76 0.6" << std::endl;
    vtkfile << "253 141 60 0.7" << std::endl;
    vtkfile << "252 78 42 0.8" << std::endl;
    vtkfile << "227 26 28 0.9" << std::endl;
    vtkfile << "189 0 38 0.9" << std::endl;
    vtkfile << "128 0 38 1.0" << std::endl;
    // done!
    vtkfile.close();
    return;
}


/*
// write to GRASS GIS 3D ASCII raster file
void writeMKDE3DtoGRASS(const vector<double> &xgrid, const vector<double> &ygrid, const vector<double> &zgrid,
                        gridFloat3D * rst3d, string fname, string nv) {
    int nX = xgrid.size();
    int nY = ygrid.size();
    int nZ = zgrid.size();
    long ijk = 0;
    double xSz = xgrid[1] - xgrid[0];
    double ySz = ygrid[1] - ygrid[0];
    double zSz = zgrid[1] - zgrid[0];
    double densTmp;
    char * fnm = new char[strlen(fname)+1];
    fnm[strlen(fname)] = 0;
    memcpy(fnm, fname, strlen(fname));

    // open file
    std::ofstream r3file;
    r3file.open (fnm);
    r3file << std::setprecision(12);

    // may need to set precision first
    // header
    r3file << "north: " << ygrid[nY - 1] + 0.5 * ySz << std::endl; // north: y.max
    r3file << "south: " << ygrid[0] - 0.5 * ySz << std::endl; // south: y.min
    r3file << "east: " << xgrid[nX - 1] + 0.5 * xSz << std::endl; // east: x.max
    r3file << "west: " << xgrid[0] - 0.5 * xSz << std::endl; // west: x.min
    r3file << "top: " << zgrid[nZ - 1] + 0.5 * zSz << std::endl; // top: z.max
    r3file << "bottom: " << zgrid[0] - 0.5 * zSz << std::endl;// bottom: z.min
    r3file << "rows: " << nY << std::endl; // rows
    r3file << "cols: " << nX << std::endl; // cols
    r3file << "levels: " << nZ << std::endl; // levels: integer

    // data
    // possibly use std::scientific or #include <iomanip> with std::setprecision(12)
    for (int eZ = rst3d->zmin; eZ <= rst3d->zmax; eZ = eZ + rst3d->zsize) {
        for (int eY = rst3d->ymin; eY <= rst3d->ymax; eY = eY + rst3d->ysize) {
            for (int eX = rst3d->xmin; eX <= rst3d->xmax; eX = eX + rst3d->xsize) {
                if (rst3d->getGridValue(eX, eY, eZ) == FLT_NO_DATA_VALUE) {
                    r3file << nv;
                } else {
                    r3file << rst3d->getGridValue(eX, eY, eZ);
                }
                if (eX == rst3d->xmax) {
                    r3file << std::endl;
                } else {
                    r3file << " ";
                }
            }
        }
    }

// done...
    r3file.close();
    return;
}



/*
// write to XDMF file
SEXP writeMKDE3DtoXDMF(vector<double> xgrid, vector<double> ygrid, vector<double> zgrid, gridFloat3D * density, char * filenameXDMF, char * filenameDAT) {
    Rcpp::NumericVector xGrid(xgrid); // cell centers in the x-dimension
    Rcpp::NumericVector yGrid(ygrid); // cell centers in the y-dimension
    Rcpp::NumericVector zGrid(zgrid); // cell centers in the z-dimension
    int nX = (long)xGrid.length();
    int nY = (long)yGrid.length();
    int nZ = (long)zGrid.length();
    long ijk = 0;
    double xSz = xGrid[1] - xGrid[0];
    double ySz = yGrid[1] - yGrid[0];
    double zSz = zGrid[1] - zGrid[0];
    double densTmp;
    std::vector<double> d = Rcpp::as<std::vector<double>>(density);
    std::string strXDMF = Rcpp::as<std::string>(filenameXDMF);
    std::string strDAT = Rcpp::as<std::string>(filenameDAT);
    char * fnmXDMF = new char[strXDMF.size()+1];
    fnmXDMF[strXDMF.size()] = 0;
    memcpy(fnmXDMF, strXDMF.c_str(), strXDMF.size());
    char * fnmDAT =new char[strDAT.size()+1];
    fnmDAT[strDAT.size()] = 0;
    memcpy(fnmDAT, strDAT.c_str(), strDAT.size());
    int nmSize = 0;
// scan backward through char array, count chars, break when hit "/"; should count \0
    for (int i = strDAT.size(); i >= 0; i--) {
        if (fnmDAT[i] == '/') {
            break;
        } else {
            nmSize++;
        }
    }
    char * binName = new char[nmSize];
    int j = 0;
    for (int i = 0; i < nmSize; i++) {
        j = strDAT.size() + 1 - nmSize + i; // CHECK THIS!!!!
        binName[i] = fnmDAT[j];
    }
// now copy name to string

// write XML wrapper
    std::ofstream xmffile;
    xmffile.open (fnmXDMF);
    xmffile << std::setprecision(12);

    xmffile << "<?xml version=\"1.0\" ?>" << std::endl;
    xmffile << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
    xmffile << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">" << std::endl;
    xmffile << "<Domain>" << std::endl;
    xmffile << "    <Grid Name=\"Mesh\" GridType=\"Uniform\">" << std::endl;
    xmffile << "        <Topology name=\"topo\" TopologyType=\"3DCoRectMesh\"" << std::endl;
    xmffile << "            Dimensions=\"" << (nZ + 1) << " " << (nY + 1) << " " << (nX + 1) << "\">" << std::endl; // levels, rows, cols?
    xmffile << "        </Topology>" << std::endl;
    xmffile << "        <Geometry name=\"geo\" Type=\"ORIGIN_DXDYDZ\">" << std::endl;
    xmffile << "            <!-- Origin -->" << std::endl;
    xmffile << "            <DataItem Format=\"XML\" Dimensions=\"3\">" << std::endl;
    xmffile << "             " << " " << (xGrid[0] - 0.5*xSz) << " " << (yGrid[0] - 0.5*ySz) << " " << (zGrid[0] - 0.5*zSz) << std::endl;
    xmffile << "            </DataItem>" << std::endl;
    xmffile << "            <!-- DxDyDz -->" << std::endl;
    xmffile << "            <DataItem Format=\"XML\" Dimensions=\"3\">" << std::endl;
    xmffile << "             " << xSz << " " << ySz << " " << zSz <<  std::endl;
    xmffile << "            </DataItem>" << std::endl;
    xmffile << "        </Geometry>" << std::endl;
    xmffile << "        <Attribute Name=\"Density\" Center=\"Cell\">" << std::endl; // need AttributeType="Scalar" or Type="Scalar" ?
    xmffile << "            <DataItem Format=\"Binary\"" << std::endl;
    xmffile << "             DataType=\"Double\"" << std::endl;
    xmffile << "             Precision=\"8\"" << std::endl;
    if (isMachineBigEndian()) {
        xmffile << "             Endian=\"Big\"" << std::endl;
    } else {
        xmffile << "             Endian=\"Little\"" << std::endl;
    }
    xmffile << "             Dimensions=\"" << nZ << " " << nY << " " << nX << "\">" << std::endl;
    xmffile << "               " << binName << std::endl;
    xmffile << "            </DataItem>" << std::endl;
    xmffile << "        </Attribute>" << std::endl;
    xmffile << "    </Grid>" << std::endl;
    xmffile << "</Domain>" << std::endl;
    xmffile << "</Xdmf>" << std::endl;

// close XML file
    xmffile.close();

// write binary data (kji order)
    std::ofstream datfile(fnmDAT, std::ios::out | std::ios::trunc | std::ios::binary); // std::ios::out | std::ios::app | std::ios::binary
    if (!datfile.is_open()) {
        cout << "Error in writeMKDE3DtoXDMF(): Output file "<< fnmDAT << " could not be opened." << std::endl;
    } else {
        for (int k = 0; k < nZ; k++) {
            for (int j = 0; j < nY; j++) {
                for (int i = 0; i < nX; i++) {
                    ijk = getLinearIndex(i, j, k, nX, nY); // i, k, k, ...
                    densTmp = d[ijk];
                    if (std::isnan(densTmp)) {
                        densTmp = 0.0;
                    }
                    datfile.write((char *)(&densTmp), sizeof(densTmp));
                }
            }
        }
        datfile.close();
    }

// done...
    return Rcpp::wrap(1);
}




SEXP writeRasterToXDMF(SEXP xgrid, SEXP ygrid, SEXP rast, SEXP filenameXDMF, SEXP filenameDAT) {
    Rcpp::NumericVector xGrid(xgrid); // cell centers in the x-dimension
    Rcpp::NumericVector yGrid(ygrid); // cell centers in the y-dimension
    int nX = (long)xGrid.length();
    int nY = (long)yGrid.length();
    long ijk = 0;
    double xSz = xGrid[1] - xGrid[0];
    double ySz = yGrid[1] - yGrid[0];
    double densTmp;
    std::vector<double> r = Rcpp::as<std::vector<double>>(rast);
    std::string strXDMF = Rcpp::as<std::string>(filenameXDMF);
    std::string strDAT = Rcpp::as<std::string>(filenameDAT);
    char * fnmXDMF = new char[strXDMF.size()+1];
    fnmXDMF[strXDMF.size()] = 0;
    memcpy(fnmXDMF, strXDMF.c_str(), strXDMF.size());
    char * fnmDAT =new char[strDAT.size()+1];
    fnmDAT[strDAT.size()] = 0;
    memcpy(fnmDAT, strDAT.c_str(), strDAT.size());
    int nmSize = 0;
// scan backward through char array, count chars, break when hit "/"; should count \0
    for (int i = strDAT.size(); i >= 0; i--) {
        if (fnmDAT[i] == '/') {
            break;
        } else {
            nmSize++;
        }
    }
    char * binName = new char[nmSize];
    int j = 0;
    for (int i = 0; i < nmSize; i++) {
        j = strDAT.size() + 1 - nmSize + i; // CHECK THIS!!!!
        binName[i] = fnmDAT[j];
    }
// now copy name to string

// write XML wrapper
    std::ofstream xmffile;
    xmffile.open (fnmXDMF);
    xmffile << std::setprecision(12);

    xmffile << "<?xml version=\"1.0\" ?>" << std::endl;
    xmffile << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
    xmffile << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">" << std::endl;
    xmffile << "<Domain>" << std::endl;
    xmffile << "    <Grid Name=\"Mesh\" GridType=\"Uniform\">" << std::endl;
    xmffile << "        <Topology name=\"topo\" TopologyType=\"2DCoRectMesh\"" << std::endl;
    xmffile << "            Dimensions=\"" << (nY + 1) << " " << (nX + 1) << "\">" << std::endl; // y, x or x, y ?
    xmffile << "        </Topology>" << std::endl;
    xmffile << "        <Geometry name=\"geo\" Type=\"ORIGIN_DXDY\">" << std::endl;
    xmffile << "            <!-- Origin -->" << std::endl;
    xmffile << "            <DataItem Format=\"XML\" Dimensions=\"2\">" << std::endl;
    xmffile << "             " << (xGrid[0] - 0.5*xSz) << " " << (yGrid[0] - 0.5*ySz) << std::endl; // y, x or x, y ?
    xmffile << "            </DataItem>" << std::endl;
    xmffile << "            <!-- DxDy -->" << std::endl;
    xmffile << "            <DataItem Format=\"XML\" Dimensions=\"2\">" << std::endl;
    xmffile << "             " << xSz << " " << ySz <<  std::endl; // y, x or x, y ?
    xmffile << "            </DataItem>" << std::endl;
    xmffile << "        </Geometry>" << std::endl;
    xmffile << "        <Attribute Name=\"Raster\" Center=\"Cell\">" << std::endl; // need AttributeType="Scalar" or Type="Scalar" ?
    xmffile << "            <DataItem Format=\"Binary\"" << std::endl;
    xmffile << "             DataType=\"Double\"" << std::endl;
    xmffile << "             Precision=\"8\"" << std::endl;
    if (isMachineBigEndian()) {
        xmffile << "             Endian=\"Big\"" << std::endl;
    } else {
        xmffile << "             Endian=\"Little\"" << std::endl;
    }
    xmffile << "             Dimensions=\"" << nY << " " << nX << "\">" << std::endl; // y, x or x, y ?
    xmffile << "               " << binName << std::endl;
    xmffile << "            </DataItem>" << std::endl;
    xmffile << "        </Attribute>" << std::endl;
    xmffile << "    </Grid>" << std::endl;
    xmffile << "</Domain>" << std::endl;
    xmffile << "</Xdmf>" << std::endl;

// close XML file
    xmffile.close();

// write binary data (kji order)
    std::ofstream datfile(fnmDAT, std::ios::out | std::ios::trunc | std::ios::binary); // std::ios::out | std::ios::app | std::ios::binary
    if (!datfile.is_open()) {
        Rcpp::Rcout << "Error in writeMKDE3DtoXDMF(): Output file "<< fnmDAT << " could not be opened." << std::endl;
    } else {
// this is effectively the same as the above
        for (int i = 0; i < r.size(); i++) {
            if (i%100000 == 0) {
                Rcpp::Rcout << "writing raster cell " << (i +1 ) << " of " << r.size() << " to file " << binName << std::endl;
            }
            densTmp = r[i];
            if (std::isnan(densTmp)) {
                densTmp = 0.0;
            }
            datfile.write((char *)(&densTmp), sizeof(densTmp));
        }
        datfile.close();
    }

// done...
    return Rcpp::wrap(1);
} */