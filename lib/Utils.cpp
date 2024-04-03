#include "Utils.hpp"
#include "particle.hpp"
#include <vector>
#include <array>
#include <stdexcept>
#include <cmath> 
#include <algorithm>
#include <numeric>
#include <fstream>
#include <fftw3.h>
#include <iomanip>

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

vector<double> correlationFunction(particle_group particles, int n_bins)
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
    int N = std::min(1000, int(particles.num_particles));
    for (int i = 0; i < N; i += 1)
    {
        for (int j = i; j < N; j += 1)
        {

            double dx = shortestDistance(particles.particles[i].position[0], particles.particles[j].position[0]);
            double dy = shortestDistance(particles.particles[i].position[1], particles.particles[j].position[1]);
            double dz = shortestDistance(particles.particles[i].position[2], particles.particles[j].position[2]);
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

void Save_Correlations_csv(const std::vector<std::vector<double>>& data, const std::vector<std::string>& columnLabels, const std::string& filename){
    std::ofstream file(filename);

    // Check if the file is open
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open the file.");
        return;
    }

    // Output column labels
    for (size_t i = 0; i < columnLabels.size(); ++i) {
        file << columnLabels[i];
        if (i < columnLabels.size() - 1) file << ",";
    }
    file << "\n";

    size_t numRows = data[0].size();

    // Iterate over each row
    for (size_t row = 0; row < numRows; ++row) {
        // Iterate over each column in the row
        for (size_t col = 0; col < data.size(); ++col) {
            // Check if the current column has enough rows; if not, output an empty string
            if (row < data[col].size()) {
                file << data[col][row];
            } else {
                file << "";
            }

            // If not the last column, add a comma separator
            if (col < data.size() - 1) file << ",";
        }
        // End of the row
        file << "\n";
    }

    file.close();
}



void PotentialSavetoTxt(std::vector<double>& potential_vec, std::vector<double>& real_vec, std::string &filename){

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

void TrajectorySavetoTxt(std::vector<double>& pos_x, std::vector<double>& pos_y, std::vector<double>& pos_z, std::vector<double>& vel_x, std::vector<double>& vel_y, std::vector<double>& vel_z, std::string &filename){

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


std::string findsigfig(double number){
    uint idx = 0;
    std::string result = "";
    std::string num_as_string = std::to_string(number);
    for (uint i = num_as_string.size() - 1; i > 0; i--){
        if (num_as_string[i] != '0'){
            idx = i;
            break;
        }
    }
    for (uint i = 0; i < num_as_string.size(); i++){
        result += num_as_string[i];
        if (i == idx){
            break;
        }
    }
    return result;
}

std::string formatREALToNDecimalPlaces(double value, uint dp) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(dp) << value;
    return oss.str();
}

/**
 * @brief: Rounds value to decimal place if value is the same as if it was rounded to one more decimal place. e.g. 2.50003 is rounded to 2.5. 2.05 is rounded to 2. 2.34340 is 2.3434. Therefore cuts off values from the first zero to the left.
 * @returns Formatted string with trailing zeros cut off from the left.
*/
std::string removeTrailingDecimalPlaces(double value, uint max_dp){
    std::optional<uint> idx;
    for (uint i = 0; i < max_dp - 1; i++){
        double rd1 =  std::stod(formatREALToNDecimalPlaces(value,i)); // converting to double for validation purpose
        double rd2 = std::stod(formatREALToNDecimalPlaces(value,i + 1));
        if (rd1 == rd2){
            idx.emplace(i);
            break;
        }
    }
    if (!(idx)){
        idx.emplace(max_dp - 1);
    }
    return formatREALToNDecimalPlaces(value,idx.value());
}