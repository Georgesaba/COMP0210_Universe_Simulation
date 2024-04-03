#pragma once

#include <vector>
#include <array>
#include <string>
#include <fftw3.h>
#include <random>
#include <optional>
#include "particle.hpp"

using std::vector;
using std::array;

/**
 * @brief Takes a buffer of fftw_complex values and outputs and image
 * Densities are integrated over the z axis to convert to 2D
 * @param density_map density values of type fftw_comlex. Imaginary component ignored.
 * @param n_cells size of buffer in each dimension; total size is n_cells*n_cells*n_cells
 * @param filename image output file path
 */
void SaveToFile(fftw_complex* density_map, const size_t n_cells, const std::string &filename);

/**
 * @brief Calculates a log radial correlation for coordinates 0 <= r < 0.5
 * Calculates pair-wise distances and counts how many fall into radial bins
 * @param positions a list of 3D positions of particles within a unit cube
 * @param n_bins the resolution of the histogram
 * @return vector<double> log of radial correlation function evenly spaced from r = 0 to 0.5
 */
vector<double> correlationFunction(particle_group particles, int n_bins);

void Save_Correlations_csv(const std::vector<std::vector<double>>& data, const std::vector<std::string>& columnLabels, const std::string& filename);

void PotentialSavetoTxt(std::vector<double>& potential_vec, std::vector<double>& real_vec, std::string &filename);

void TrajectorySavetoTxt(std::vector<double>& pos_x, std::vector<double>& pos_y, std::vector<double>& pos_z, std::vector<double>& vel_x, std::vector<double>& vel_y, std::vector<double>& vel_z, std::string &filename);


/**
 * @brief: Removes trailing zeros in values starting from the most amount of decimal places. Cuts off values from the first zero to the right.
 * @returns: Formatted string with trailing zeros cut off from the right.
*/
std::string findsigfig(double number);
std::string formatREALToNDecimalPlaces(double value, uint dp = 3);
std::string removeTrailingDecimalPlaces(double value, uint max_dp = 3);