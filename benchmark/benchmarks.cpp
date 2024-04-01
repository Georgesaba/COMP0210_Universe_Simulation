#include <iostream>
#include <vector>
#include <string>
#include <omp.h>
#include <chrono>
#include "Simulation.hpp"

class BenchmarkData
{
public:
    BenchmarkData(std::string benchmark_name, int threads) : name(benchmark_name), num_threads(threads) {}
    double time;
    std::string name;
    int num_threads;
    //use for additional context such as number of particles / cells 
    std::string info;

    void start()
    {
        t1 = std::chrono::high_resolution_clock::now();
    }

    void finish()
    {
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        time = (t2 - t1).count() / 1e9;
    }

    std::chrono::high_resolution_clock::time_point t1;
};

std::ostream& operator<<(std::ostream &os, const BenchmarkData& b)
{
    std::cout << "Benchmarking " << b.name << " with " << b.num_threads << " threads." << std::endl;
    std::cout << "Time = " << b.time << std::endl;
    std::cout << "Info: " << b.info << std::endl;
    return os;
}

int main()
{
    uint average_particles_per_cell = 10;
    uint num_cells = 201;
    uint num_particles = num_cells * num_cells * num_cells * average_particles_per_cell;
    double mass = 10.0 * 10.0 * 10.0 * 10.0 * 10.0/num_particles;
    particle_group particles(mass, num_particles, 42); // initialising particle group as is passed by value to the constructor
                                                       // does not need to be re-initiated
    
    uint max_threads =  omp_get_max_threads(); // my CPU has 16 threads and 8 cores

    std::string info = "The number of cells per length of the box is " + std::to_string(num_cells) + " and the number of particles is " + std::to_string(num_particles) + ".";
    
    std::vector<BenchmarkData> density_benches;
    std::vector<BenchmarkData> potential_benches;
    std::vector<BenchmarkData> gradient_benches;
    std::vector<BenchmarkData> update_par_benches;
    std::vector<BenchmarkData> expansion_benches;

    for (uint i = 1; i <= max_threads; i++){
        omp_set_num_threads(i);
        double width = 100.0;
    
        Simulation sim(1.5, 0.01, particles, width, num_cells, 1.02); // initialise new simulation each time so runs are identical.
        
        // density evaluation
        BenchmarkData density_bench("Density Calculation", i);
        density_bench.start();
        sim.fill_density_buffer();
        density_bench.finish();
        density_bench.info = info;
        density_benches.push_back(density_bench);

        // potential evaluation
        BenchmarkData potential_bench("Potential Calculation", i);
        potential_bench.start();
        sim.fill_potential_buffer();
        potential_bench.finish();
        potential_bench.info = info;
        potential_benches.push_back(potential_bench);

        // calculate gradient
        BenchmarkData gradient_bench("Gradient Calc", i); // gradient calc and particle update are joined
        gradient_bench.start();
        sim.calculate_gradient(sim.get_potential_buffer());
        gradient_bench.finish();
        gradient_bench.info = info;
        gradient_benches.push_back(gradient_bench);

        // update particles
        BenchmarkData update_par_bench("Particle Update and Gradient Calc", i); // gradient calc and particle update are joined
        update_par_bench.start();
        sim.update_particles();
        update_par_bench.finish();
        update_par_bench.info = info;
        update_par_benches.push_back(update_par_bench);

        // box_expansion
        BenchmarkData expansion_bench("Expansion Calculation", i);
        expansion_bench.start();
        sim.box_expansion();
        expansion_bench.finish();
        expansion_bench.info = info;
        expansion_benches.push_back(expansion_bench);
    }
    for (uint i = 0; i < density_benches.size(); i++){
        std::cout << density_benches[i] << std::endl;
    }
    for (uint i = 0; i < potential_benches.size(); i++){
        std::cout << potential_benches[i] << std::endl;
    }
    for (uint i = 0; i < gradient_benches.size(); i++){
        std::cout << gradient_benches[i] << std::endl;
    }
    for (uint i = 0; i < update_par_benches.size(); i++){
        std::cout << update_par_benches[i] << std::endl;
    }
    for (uint i = 0; i < expansion_benches.size(); i++){
        std::cout << expansion_benches[i] << std::endl;
    }
    return 0;
}