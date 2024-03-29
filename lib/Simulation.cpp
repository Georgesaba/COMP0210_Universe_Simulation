#include "Simulation.hpp"
#include "Utils.hpp"
#include "particle.hpp"

Simulation::Simulation(double t_max, double t_step, particle_group collection, double W, uint num_cells, double e_factor) : 
                        time_max(t_max), time_step(t_step), particle_collection(collection), box_width(W), number_of_cells(num_cells),
                         expansion_factor(e_factor){}

void Simulation::run()
{

}