#pragma once
#include "particle.hpp"
#include <fftw3.h>
#include <vector>
#include <optional>

class Simulation
{
public:
    Simulation(double t_max, double t_step, particle_group collection, double W, uint num_cells, double e_factor);  
    
    /**
     * @brief Run a particle mesh simulation from t=0 to t_max
     * 
     */
    void run(std::optional<std::string> output_folder = std::nullopt);

    void fill_density_buffer();
    void fill_potential_buffer();
    std::vector<std::vector<std::vector<std::array<double, 3>>>> calculate_gradient(const fftw_complex * potential);
    void update_particles();
    
    const fftw_complex * get_density_buffer() const;
    const fftw_complex * get_potential_buffer() const;
    const particle_group & get_particle_collection() const;
    void box_expansion();

    ~Simulation();

    private:
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