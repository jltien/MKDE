#include "pairwiseFunctions.h"
// Preconditions: Location times in data are sorted in order of increasing time, no NAN values

// Associates the observed points with the interpolated coordinates
std::vector<double> assignLocationIndexToTimeGrid(const std::vector<time_t> grid_times,
                                               const std::vector<time_t> location_times, double dt_max) {
    double NA_INDEX = NAN;
    double m = grid_times.size();
    double n = location_times.size();
    double j = 0;
    std::vector<double> location_indexes(m, NA_INDEX);

    for (int i = 0; i < m; i++) {
        while (location_times[j + 1] < grid_times[i]) {
            if (j < (n - 2)) { // there has to be a location at j+1
                j++;
            } else {
                break; // can't increment j any more
            }
        }
        if ((location_times[j] <= grid_times[i]) && (location_times[j+1] > grid_times[i])) {
            if ((location_times[j+1] - location_times[j]) <= dt_max) {
                location_indexes[i] = j;
            }
        }
        else if (location_times[j] <= grid_times[i]) {
            location_indexes[i] = j + 1;
        }
    }

    return location_indexes;
}


std::vector<pointIn3D> interpolateCoordinateOnTimeGrid(const std::vector<time_t> & grid_times, const std::vector<double> & location_indexes,
                                                       const std::vector<time_t> & location_times, const std::vector<pointIn3D> & location_xyz,
                                                       const std::vector<double> moveVarXY, const std::vector<double> moveVarZ,
                                                       const std::vector<double> obsVarXY, const std::vector<double> obsVarZ) {
    double NA_INDEX = NAN;
    int m = grid_times.size();
    int n = location_times.size();
    double alpha_j, sig2xy, sig2z, dt;
    std::vector<double> res_x(m, NA_INDEX);
    std::vector<double> res_y(m, NA_INDEX);
    std::vector<double> res_z(m, NA_INDEX);
    std::vector<pointIn3D> res_xyz;

    for (int i = 0; i < grid_times.size(); i++) {
        double j = location_indexes[i];
        if (j != NA_INDEX && j <= (n-2)) {

            if (location_times[j + 1] - location_times[j] != 0) {
                alpha_j = (grid_times[i] - location_times[j]) * 1.0 / (location_times[j + 1] - location_times[j]);
            }
            else {
                alpha_j = 0;
            }
            res_x[i] = location_xyz[j].x + alpha_j * (location_xyz[j + 1].x - location_xyz[j].x);
            res_y[i] = location_xyz[j].y + alpha_j * (location_xyz[j + 1].y - location_xyz[j].y);
            res_z[i] = location_xyz[j].z + alpha_j * (location_xyz[j + 1].z - location_xyz[j].z);

            dt = grid_times[j + 1] - grid_times[j];
            sig2xy = dt * alpha_j * (1.0 - alpha_j) * moveVarXY[j] + obsVarXY[j] * (1.0 - alpha_j)
                                                                     * (1.0 - alpha_j) + obsVarXY[j + 1] * alpha_j * alpha_j;
            sig2z = dt * alpha_j * (1.0 - alpha_j) * moveVarZ[j] + obsVarZ[j] * (1.0 - alpha_j)
                                                                     * (1.0 - alpha_j) + obsVarZ[j + 1] * alpha_j * alpha_j;
        }
     //   std::cout << alpha_j << std::endl;
        res_xyz.push_back(pointIn3D(res_x[i], res_y[i], res_z[i], grid_times[i], alpha_j, sig2xy, sig2z));
    }

    // return alpha_j with each tuple
    return res_xyz;
}


std::vector<double> euclideanDistance(const std::vector<pointIn3D> & xyz0, const std::vector<pointIn3D> & xyz1, bool use_z) {

    use_z = false;
    int n0 = xyz0.size();
    int n1 = xyz1.size();
    std::vector<double> res_dist(n0, NAN);

    if (n0 == n1) {
        for (int i = 0; i < n0; i++) {
            if (use_z) {
                res_dist[i] = sqrt((xyz0[i].x - xyz1[i].x) * (xyz0[i].x - xyz1[i].x) +
                                   (xyz0[i].y - xyz1[i].y) * (xyz0[i].y - xyz1[i].y) +
                                   (xyz0[i].z - xyz1[i].z) * (xyz0[i].z - xyz1[i].z));
            }
            else {
                res_dist[i] = sqrt((xyz0[i].x - xyz1[i].x) * (xyz0[i].x - xyz1[i].x) +
                                   (xyz0[i].y - xyz1[i].y) * (xyz0[i].y - xyz1[i].y));
            }
        }
    }
    else {
        std::cout << "Error in xyz inputs!" << std::endl;
    }
    return res_dist;
}


// Writes result from Interpolate and EuclideanDistance to csv file
void printInterpolateAndEuclidean(const std::string filename, const std::vector<pointIn3D> res_xyz, const std::vector<double> res_dist) {
    std::ofstream csv;
    csv.open(filename);
    for (int i = 0; i < res_xyz.size(); i++) {
        csv << res_xyz[i].x << " " << res_xyz[i].y << " " << res_xyz[i].z
            << " " << res_xyz[i].time << " " << res_xyz[i].alpha_j << std::endl;
    }
    for (int i = 0; i < res_dist.size(); i++) {
        csv << res_dist[i] << std::endl;
    }
}


// Writes result from EuclideanDistance to csv file
void printEuclidean(const std::string filename, const std::vector<double> res_dist) {
    std::ofstream csv;
    csv.open(filename);
    for (int i = 0; i < res_dist.size(); i++) {
        csv << res_dist[i] << std::endl;
    }
}