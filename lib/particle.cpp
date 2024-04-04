#include "particle.hpp"
#include <random>
#include <stdexcept>
#include <iostream>



particle::particle(const std::array<double, 3> &initial_position){
    for (double pos: initial_position){
        if (pos > 1 || pos < 0){
            throw std::range_error("Error - Element in vector " + std::to_string(pos) + " is outside of boundary conditions!");
        }
        else if (pos == 1){
            pos = 0; //For human input, almost zero change of rng 1 value being obtained
        }
    }
    position = initial_position;
}


particle_group::particle_group(double mass, uint num_particles, const std::vector<std::array<double,3>> &positions) : 
                            mass(mass), num_particles(num_particles) 
{
    if (mass <= 0){
        throw std::invalid_argument("Error - The particle masses must be larger than 0!");
    }
    if (num_particles != positions.size()){
        throw std::invalid_argument("Error - The number of particles does not match the size of the given position vector!");
    }
    for (uint i = 0; i < num_particles; i++){
        particles.push_back(particle(positions[i]));
    }
}


particle_group::particle_group(double mass, uint num_particles, uint random_seed) :
                            mass(mass), num_particles(num_particles)
{
    if (mass <= 0){
        throw std::invalid_argument("Error - The particle masses must be larger than 0!");
    }
    if (num_particles > 10000000000){
        std::cerr << "Warning - More than 10,000,000,000 particles have been generated! This may negatively impact performance." << std::endl;
    }
    std::default_random_engine generator(random_seed);
    std::uniform_real_distribution<double> initial_dist(0, 1);
    std::array<double, 3> initial_position;
    for (uint i = 0; i < num_particles; i++){
        for (uint j = 0; j < 3; j++){
            initial_position[j] = initial_dist(generator);
        }
        particles.push_back(particle(initial_position));
    }
}

size_t particle_group::get_num_particles(){
    num_particles = particles.size();
    return num_particles;
}