#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <cmath>
#include "Simulation.hpp"

using namespace Catch::Matchers;

TEST_CASE("Test density calculation function for no particles","[Density_Calc]"){
    double mass = 0.01;
    double width = 1;
    uint number_particles = 0;
    particle_group particles(mass, number_particles, {});
    uint num_cells = 100;
    double cell_width = width/num_cells;
    Simulation sim(10, 0.1, particles, width, num_cells, 2);
    sim.fill_density_buffer();
    const fftw_complex* density_buffer = sim.get_density_buffer();
    for (uint i = 0; i < num_cells * num_cells * num_cells; i++){
        CHECK(density_buffer[i][1] == 0);
        CHECK_THAT(density_buffer[i][0], WithinRel(0,1e-10));
    }
}


TEST_CASE("Test density calculation function for a single particle","[Density_Calc]"){
    double mass = 0.01;
    double width = 1;
    uint number_particles = 1;
    particle_group particles(mass, number_particles, {{0.45, 0.45, 0.45}});
    uint num_cells = 100;
    double cell_width = width/num_cells;
    Simulation sim(10, 0.1, particles, width, num_cells, 2);
    sim.fill_density_buffer();
    const fftw_complex* density_buffer = sim.get_density_buffer();
    for (uint i = 0; i < num_cells * num_cells * num_cells; i++){
        CHECK(density_buffer[i][1] == 0);
        if (45 + num_cells * (45 + num_cells * 45) == i){
            CHECK_THAT(density_buffer[i][0], WithinRel(mass/(cell_width * cell_width * cell_width),1e-10));
        }
        else{
            CHECK_THAT(density_buffer[i][0], WithinRel(0,1e-10));
        }
    }
}

TEST_CASE("Test density calculation function for multiple particles","[Density_Calc]"){
    double mass = 0.01;
    double width = 1;
    uint number_particles = 10;
    uint num_cells = 10;
    double cell_width = width/num_cells;
    std::vector<std::array<double, 3>> particle_pos = {{0.45,0.45,0.45},{0.42,0.44,0.43},{0.45,0.49,0.41},{0.1,0.2,0.3},{0.6,0.6,0.6},
    {0.9,0.9,0.9},{0.8,0.8,0.8},{0.5,0.5,0.5},{0.2,0.2,0.2},{0.3,0.2,0.1}};
    particle_group particles(mass, number_particles, particle_pos);
    Simulation sim(10, 0.1, particles, width, num_cells, 2);
    sim.fill_density_buffer();
    const fftw_complex* density_buffer = sim.get_density_buffer();
    for (uint i = 0; i < num_cells * num_cells * num_cells; i++){
        // set conditions for single particle cells
        bool condition1 = (1 + num_cells * (2 + num_cells * 3) == i);
        bool condition2 = (6 + num_cells * (6 + num_cells * 6) == i);
        bool condition3 = (9 + num_cells * (9 + num_cells * 9) == i);
        bool condition4 = (8 + num_cells * (8 + num_cells * 8) == i);
        bool condition5 = (5 + num_cells * (5 + num_cells * 5) == i);
        bool condition6 = (2 + num_cells * (2 + num_cells * 2) == i);
        bool condition7 = (3 + num_cells * (2 + num_cells * 1) == i);
        
        CHECK(density_buffer[i][1] == 0);
        if (4 + num_cells * (4 + num_cells * 4) == i){
            CHECK_THAT(density_buffer[i][0], WithinRel(3 * mass/(cell_width * cell_width * cell_width),1e-10));
        }
        else if ( condition1|| condition2 || condition3 || condition4 || condition5 || condition6 || condition7)
        {
            CHECK_THAT(density_buffer[i][0], WithinRel(mass/(cell_width * cell_width * cell_width),1e-10));
        }
        else{
            CHECK_THAT(density_buffer[i][0], WithinRel(0,1e-10));
        }
    }
}




/**
 * @brief Fill out this test function by filling in the TODOs
 * Tests the calculation of the gravitational potential due to a single particle
 * Potential is examined along the x axis throught the centre of the cube
 */
// TEST_CASE("Test potential function for single particle", "Potential Tests")
// {
//     //One particle in the centre of the box with mass 0.01
//     //TODO: Set up single particle with position {0.5, 0.5, 0.5}
//     double mass = 0.01;

//     // These properties should apply to your box
//     double width = 100;
//     int ncells = 101;

//     //Declare simulation with particle setup, width, and number of cells. 
//     //TODO: Simulation sim(...);
//     //TODO: calculate density
//     //TODO: calculate potential

//     double w_c = width / ncells;
//     // Look at the potential function along the x axis
//     for (int i = 0; i < 101; i++)
//     {
//         int j = 50;
//         int k = 50;
//         if (i == 50 && j == 50 && k == 50)
//         {
//             //ignore the centre where point potential would be singular
//             continue;
//         }
//         else
//         {
//             // ignore imaginary component
//             // approximate potential for multiple sources due to periodicity
//             double dx1 = std::abs(i - 50) * w_c;
//             double dx2 = std::abs(i - 151) * w_c;
//             double dx3 = std::abs(i + 51) * w_c;
//             double pot; //TODO = your potential function at indices (i,j,k)
//             double expected_pot =  -mass *(1/ dx1 + 1/dx2 + 1/dx3);
//             double diff = pot - expected_pot;
//             REQUIRE_THAT(pot, WithinRel(expected_pot, 0.3));
//         }
//     }
// }