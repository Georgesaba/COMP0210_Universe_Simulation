#pragma once
#include "particle.hpp"
#include <fftw3.h>

class Simulation
{
public:
    Simulation(double t_max, double t_step, particle_group collection, double W, uint num_cells, double e_factor);  
    
    /**
     * @brief Run a particle mesh simulation from t=0 to t_max
     * 
     */
    void run();

    void fill_density_buffer();
    void fill_potential_buffer();
    
    const fftw_complex * get_density_buffer() const;
    const fftw_complex * get_potential_buffer() const;

    ~Simulation();

    private:
    std::vector<std::vector<std::vector<uint>>> bin_particles();

    double time_max;
    double time_step;
    particle_group particle_collection;
    double box_width;
    uint number_of_cells;
    double expansion_factor;

    fftw_complex * density_buffer;
    fftw_complex * potential_buffer;
    fftw_complex * k_space_buffer;
    fftw_plan forward_plan;
    fftw_plan backward_plan;
};