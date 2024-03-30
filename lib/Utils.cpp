#include "Utils.hpp"
#include <vector>
#include <array>
#include <stdexcept>
#include <cmath> 
#include <algorithm>
#include <numeric>
#include <fstream>
#include <fftw3.h>

using std::fstream;
using std::vector;
using std::string;

void SaveToFile(fftw_complex* density_map, const size_t n_cells, const string &filename)
{
    //Write the file header
    fstream image_file;
    image_file.open(filename, fstream::out);
    if(!image_file)
    {
        throw std::runtime_error("File failed to open");
    }
    image_file << "P3\n" << n_cells << " " << n_cells << "\n255\n";

    vector<double> density_xy(n_cells*n_cells);

    for(size_t i = 0; i < n_cells*n_cells; i++)
    {
        density_xy[i] = 0;
    }
    for(size_t i = 0; i < n_cells; i++)
    {
        for(size_t j = 0; j < n_cells; j++)
        {
            for(size_t k = 0; k < n_cells; k++)
            {
                density_xy[i*n_cells + j] += density_map[k + n_cells*(j + n_cells*i)][0];
            }
        }
    }
    auto max = std::max_element(density_xy.begin(), density_xy.end());
    double mean = std::accumulate(density_xy.begin(), density_xy.end(), 0.0) / (n_cells*n_cells);
    double norm = 255/mean;
    for(size_t i = 0; i < n_cells*n_cells; i++)
    {
        density_xy[i] *= norm;
    }

    for (size_t i = 0; i < n_cells; i++)
    {
        for (size_t j = 0; j < n_cells; j++)
        {
            // r
            image_file << std::min(static_cast<int>(density_xy[i * n_cells + j]), 255) << " ";
            // g
            image_file << std::min(std::max(static_cast<int>(density_xy[i * n_cells + j] - 255), 0), 255) << " ";
            // b
            image_file << std::min(std::max(static_cast<int>(density_xy[i * n_cells + j] - 550), 0), 255) << " ";

            image_file << "\n";
        }
    }
}

vector<double> correlationFunction(vector<array<double,3>> positions, int n_bins)
{
    if(n_bins <= 0)
    {
        throw std::runtime_error("Correlation function requires a positive definite number of bins.");
    }

    // Container for correlation function
    std::vector<double> CR(n_bins, 0.0);

    auto shortestDistance = [](double x1, double x2)
    {
        double d = std::abs(x1 - x2);
        return d < 0.5 ? d : (1 - d);
    };
    
    // Only take a limited sample of positions if there are too many 
    int N = std::min(1000, int(positions.size()));
    for (int i = 0; i < N; i += 1)
    {
        for (int j = i; j < N; j += 1)
        {

            double dx = shortestDistance(positions[i][0], positions[j][0]);
            double dy = shortestDistance(positions[i][1], positions[j][1]);
            double dz = shortestDistance(positions[i][2], positions[j][2]);
            double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            if(r < 0.5)  // within the 0.5 radius sphere to avoid edge effects from cube
            {
                int idx = static_cast<int>(r * n_bins * 2);
                CR[idx] += 1 / (N * 4 * M_PI * r * r);
            }
        }
    }
    for (auto &c : CR)
    {
        c = std::log(c);
    }
    return CR;
}

void PotentialSavetoTxt(std::vector<double> potential_vec, std::vector<double> real_vec, std::string &filename){

    // Open a file in write mode
    std::ofstream outFile(filename);
    if (!outFile) {
        throw std::runtime_error("Error opening file for writing.");
    }

    // Assuming all vectors are of the same length
    size_t numElements = potential_vec.size(); // or std::min({vector1.size(), vector2.size(), vector3.size()}) for safety

    for (size_t i = 0; i < numElements; ++i) {
        // Write an element from each vector to the file, separated by a space
        outFile << potential_vec[i] << " " << real_vec[i]  << "\n";
    }

    // Clean up
    outFile.close();
}

void TrajectorySavetoTxt(std::vector<double> pos_x, std::vector<double> pos_y, std::vector<double> pos_z, std::vector<double> vel_x, std::vector<double> vel_y, std::vector<double> vel_z, std::string &filename){

    // Open a file in write mode
    std::ofstream outFile(filename);
    if (!outFile) {
        throw std::runtime_error("Error opening file for writing.");
    }

    // Assuming all vectors are of the same length
    size_t numElements = pos_x.size(); // or std::min({vector1.size(), vector2.size(), vector3.size()}) for safety

    for (size_t i = 0; i < numElements; ++i) {
        // Write an element from each vector to the file, separated by a space
        outFile << pos_x[i] << " " << pos_y[i] << " " << pos_z[i] << " " << vel_x[i] << " " << vel_y[i] << " " << vel_z[i] << "\n";
    }

    // Clean up
    outFile.close();
}