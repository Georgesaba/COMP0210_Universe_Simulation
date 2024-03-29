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
private:
    double time_max;
    double time_step;
    particle_group particle_collection;
    double box_width;
    uint number_of_cells;
    double expansion_factor;
};