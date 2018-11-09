//
// Created by Joyce Tien on 10/5/18.
//

#include "pairwiseFunctions.h"

std::vector<double> assignLocationIndexToTimeGrid(std::vector<double> & grid_times,
                                               std::vector<double> & location_times, double dt_max) {
    dt_max = 10.0;
    int NA_INDEX = -1;
    double m = grid_times.size();
    double n = location_times.size();
    double j = 0;
    std::vector<double> location_indexes(m, NAN);

    for (int i = 0; i < m; i++) {
        while (location_times[j + 1] < grid_times[i]) {
            if (j < (n - 2)) { // there has to be a location at j+1
                j++;
            } else {
                break; // can't increment j any more
            }
        }
        if (location_times[j] <= grid_times[i] && location_times[j+1] > grid_times[i]) {
            if ((location_times[j+1] - location_times[j]) < dt_max) {
                location_indexes[i] = j;
            }
        }
    }
    return location_indexes;
}

std::vector<Tuple> interpolateCoordinateOnTimeGrid(std::vector<double> & grid_times, std::vector<double> & location_indexes,
                                                 std::vector<double> & location_times, std::vector<Tuple> & location_xyz) {
    int NA_INDEX = -1;
    int m = grid_times.size();
    int n = location_times.size();
    std::vector<double> res_x(m, NAN);
    std::vector<double> res_y(m, NAN);
    std::vector<double> res_z(m, NAN);
    std::vector<Tuple> res_xyz;

    for (int i = 0; i < grid_times.size(); i++) {
        int j = location_indexes[i];
        if (j != NA_INDEX && j <= (n-2)) {
            double delta_j = (grid_times[i] - location_times[j])/(location_times[j+1] - location_times[j]);
            res_x[i] = location_xyz[j].x + delta_j*(location_xyz[j+1].x - location_xyz[j].x);
            res_y[i] = location_xyz[j].y + delta_j*(location_xyz[j+1].y - location_xyz[j].y);
            res_z[i] = location_xyz[j].z + delta_j*(location_xyz[j+1].z - location_xyz[j].z);
        }
    }
    for (int i = 0; i < grid_times.size(); i++) {
        res_xyz.push_back(Tuple(res_x[i], res_y[i], res_z[i]));
    }
    return res_xyz;
}

std::vector<double> euclideanDistance(std::vector<Tuple> & xyz0, std::vector<Tuple> & xyz1, bool use_z) {
    use_z = false;
    int n0 = xyz0.size();
    int n1 = xyz1.size();
    if (n0 == n1) {
        std::vector<double> res_dist(n0, NAN);
        for (int i = 0; i < n0; i++) {
            if (use_z) {
                res_dist[i] = sqrt((xyz0[i].x - xyz1[i].x)*(xyz0[i].x - xyz1[i].x) +
                                           (xyz0[i].y - xyz1[i].y)*(xyz0[i].y - xyz1[i].y) +
                                           (xyz0[i].z - xyz1[i].z)*(xyz0[i].z - xyz1[i].z));
            }
            else {
                res_dist[i] = sqrt((xyz0[i].x - xyz1[i].x)*(xyz0[i].x - xyz1[i].x) +
                                           (xyz0[i].y - xyz1[i].y)*(xyz0[i].y - xyz1[i].y));
            }
        }
        return res_dist;
    }
    return {};
}