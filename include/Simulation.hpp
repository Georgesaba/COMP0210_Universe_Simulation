#pragma once
#include "particle.hpp"
#include <fftw3.h>
#include <vector>
#include <optional>

/**
 * @brief: Class that takes an initial distribution of particles and then uses the particle mesh method to simulate the trajectories of N bodies due to the resultant gravitational field.
 * Calculates the gravitational potential at each point in the cubic mesh and then evaluates the acceleration due to gravity for each cell. Updates particle positions based on this gravity.
*/
class Simulation
{
public:
    /**
     * @brief Constructor for Simulation class. Allocates memory in heap for FFT plans for forward and backwards fast fourier transform and respective buffers. Initialises member variables of class.
     * @param t_max: Time at which Simulation terminates.
     * @param t_step: Timestep which separates each moment that the Simulation evaluates particle positions for.
     * @param collection: Particle_group instance that contains the initial distribution of particles to be passed to the Simulation.
     * @param num_cells: Number of cells per length of the cubic box the Simulation runs in.
     * @param e_factor: Expansion factor - Factor by which the simulation is scaled by every iteration.
    */
    Simulation(double t_max, double t_step, particle_group collection, double W, uint num_cells, double e_factor);  
    
    /**
     * @brief Run a particle mesh simulation from t=0 to t_max in slices separated by dt.
     * @param output_folder string containing the output folder that the simulation images will be saved to. Optional argument that defaults to a std::nullopt object and results in no saved plots.
     */
    void run(std::optional<std::string> output_folder = std::nullopt);

    /**
     * @brief: Calculates the density of every cell in the cubic box. Stores in the density buffer array with type fftw_complex.
    */
    void fill_density_buffer();

    /**
     * @brief: Evaluates the gravitational potential of every cell in the cubic box. Stores in the density buffer array with type fftw_complex.
     * Evaluates Fast Fourier Transform of density buffer, applies factors and performs back transformation.
    */
    void fill_potential_buffer();
    std::vector<std::vector<std::vector<std::array<double, 3>>>> calculate_gradient(const fftw_complex * potential);
    
    /**
     * @brief: Given cell graviational potential calculates the acceleration due to gravity in every direction in each cell of the box.
     * Applies acceleration to constant acceleration equations of motion to evaluate updated velocities and acceleration.
    */
    void update_particles();
    
    /**
     * @brief: Applies expansion factor to width of box and velocity of every particle.
    */
    void box_expansion();

    /**
     * @brief: Destructor deallocates fftw_complex * c array memory in heap and deallocates memory used to store FFT plans.
    */
    ~Simulation();

    const fftw_complex * get_density_buffer() const;
    const fftw_complex * get_potential_buffer() const;
    const particle_group & get_particle_collection() const;

    private:
    double time_max;
    double time_step;
    particle_group particle_collection;
    double box_width;
    uint number_of_cells;
    double expansion_factor;

    fftw_complex * density_buffer; // buffers and plans
    fftw_complex * potential_buffer;
    fftw_complex * k_space_buffer;
    fftw_plan forward_plan;
    fftw_plan backward_plan;
};