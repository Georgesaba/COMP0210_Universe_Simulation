#pragma once

#include <vector>
#include <array>
#include <random>

/**
 * @brief: Class designed to hold position and velocity data for single particle.
*/
class particle
{
public:
    /**
     * @brief: Constructor for particle class assigning particle positions manually Contains error handling for out of range inputs.
     * @param initial_position: 1D array of length 3 to hold coordinates in x, y and z directions with respect to the unit cube describing the box.
    */
    particle(const std::array<double, 3> &initial_position);
    std::array<double, 3> position;
    std::array<double, 3> velocity = {0, 0, 0};
};

/**
 * @brief: Class designed to hold collection of particle objects.
*/
class particle_group
{
    public:

    /**
     * @brief: Constructor for particle_group class allowing for uniform random initialisation of particle positions.
     * @param mass: Mass of each particle.
     * @param num_particles: Number of particles to be created in the group.
     * @param random_seed: Random seed that will be applied to the STL standard library default random number generator following the uniform distribution.
    */
    particle_group(double mass, uint num_particles, uint random_seed);
    
    /**
     * @brief: Constructor for particle_group class allowing for manual assignment of particle positions. Contains error handling to check if inputted number of particles value is correct
     * @param mass: Mass of each particle.
     * @param num_particles: Number of particles to be created in the group.
     * @param positions: Vector of length 3 arrays that contain the coordinates in the unit cube in all 3 directions of cartesian space.
    */
    particle_group(double mass, uint num_particles, const std::vector<std::array<double,3>> &positions);

    size_t get_num_particles();

    double mass;
    std::vector<particle> particles;

    private:
    uint num_particles;
};
