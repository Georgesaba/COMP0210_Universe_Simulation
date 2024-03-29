#pragma once

#include <vector>
#include <array>
#include <random>

/**
 * @brief: Class designed to hold collection of particle objects.
*/
class particle_group
{
public:
    particle_group(double width, double mass, uint num_particles, uint random_seed);
    particle_group(double width, double mass, uint num_particles, const std::vector<std::array<double,3>> &positions);
    double mass;
    uint num_particles;
    std::vector<particle> particles;
};


/**
 * @brief: Class designed to hold position and velocity data for single particle.
*/
class particle
{
public:
    particle(const std::array<double, 3> &initial_position);
    std::array<double, 3> position;
    std::array<double, 3> velocity = {0, 0, 0};
};