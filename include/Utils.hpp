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

/**
 * @brief: Saves log radial correlation for coordinates 0 <= r < 0.5 (output of correlationFunction)
 * Saves results to csv file.
 * @param data: 2D std::vector that contains the radial correlation function for each bin.
 * @param columnLabels: vector of strings that contain the expansion factor the correlation function was computed for.
 * @param filename: string of the file path and name that the csv will be saved to.
*/
void Save_Correlations_csv(const std::vector<std::vector<double>>& data, const std::vector<std::string>& columnLabels, const std::string& filename);

/**
 * @brief: Function that saves two vectors of doubles to a txt file. This was used to plot the potential evalutated by the simulation to check that it had suitable shape.
 * @param filename: String that contains the path the txt file is saved to.
*/
void PotentialSavetoTxt(std::vector<double>& potential_vec, std::vector<double>& real_vec, std::string &filename);

/**
 * @brief: Function that saves particle trajectories to a txt file. This was used to check that the cpp tests were evaluating the trajectories properly when searching for a 0 mean particle trajectory.
*/
void TrajectorySavetoTxt(std::vector<double>& pos_x, std::vector<double>& pos_y, std::vector<double>& pos_z, std::vector<double>& vel_x, std::vector<double>& vel_y, std::vector<double>& vel_z, std::string &filename);


/**
 * @brief: Removes trailing zeros in values starting from the most amount of decimal places. Cuts off values from the first zero to the right.
 * Keeps decimal point however useful to prevent small numbers from being rounded to 0 as in the removeTrailingDecimalPlace function.
 * @returns: Formatted string with trailing zeros cut off from the right however keeps decimal point.
*/
std::string findsigfig(double number);

/**
 * @brief: Helper function to removeTrailingDecimalPlaces that formats a double to a string with a set number of decimal places.
 * @param dp: Number of decimal places that the string is rounded to.
*/
std::string formatREALToNDecimalPlaces(double value, uint dp = 3);

/**
 * @brief: Rounds value to decimal place if value is the same as if it was rounded to one more decimal place. e.g. 2.50003 is rounded to 2.5. 2.05 is rounded to 2. 2.34340 is 2.3434. 
 * Therefore cuts off values from the first zero to the left. An example of failure is rounding 0.01 to 0.
 * @param value: double that needs to be converted to a rounded string.
 * @param max_dp: Threshold for maximum number of decimal places that the value will be rounded to.
 * @returns Formatted string with trailing zeros cut off from the left.
*/
std::string removeTrailingDecimalPlaces(double value, uint max_dp = 3);