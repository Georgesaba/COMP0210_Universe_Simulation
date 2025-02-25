#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <cmath>
#include "Simulation.hpp"
#include "Utils.hpp"
#include <iostream>
#include <algorithm>

using namespace Catch::Matchers;

TEST_CASE("Test Simulation constructor for error handling with invalid arguments", "[Simulation_Constructor]"){
    double mass = 0.01;
    double width = 1;
    uint number_particles = 0;
    particle_group particles(mass, number_particles, {});
    uint num_cells = 100;
    double cell_width = width/num_cells;
    REQUIRE_THROWS(Simulation(-1, 0.1, particles, width, num_cells, 2));
    REQUIRE_THROWS(Simulation(1, -0.1, particles, width, num_cells, 2));
    REQUIRE_THROWS(Simulation(1, 0.1, particles, width, num_cells,-2));
    REQUIRE_THROWS(Simulation(1, 0.1, particles, width,-1 * num_cells,2));
    REQUIRE_THROWS(Simulation(1, 0.1, particles, -1 * width, num_cells,2));
}

TEST_CASE("Test particle constructor for error handling with invalid arguments", "[particle_constructor]"){
    REQUIRE_THROWS(particle({-1,-1,-1}));
    REQUIRE_THROWS(particle({2,2,2}));
}

TEST_CASE("Test particle group constructor for error handling with invalid arguments","[particle_constructor]"){
    uint number_particles = 0;
    REQUIRE_THROWS(particle_group(-1, number_particles, {}));
    REQUIRE_THROWS(particle_group(1, number_particles, {{1,1,1}}));
}


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
TEST_CASE("Test potential function for single particle", "[Potential_Calc]")
{
    //One particle in the centre of the box with mass 0.01
    //TODO: Set up single particle with position {0.5, 0.5, 0.5}
    double mass = 0.01;
    particle_group particles(mass, 1, {{0.5, 0.5, 0.5}});

    // These properties should apply to your box
    double width = 100;
    int ncells = 101;

    //Declare simulation with particle setup, width, and number of cells. 
    Simulation sim(10, 0.1, particles, width, ncells, 1);
    sim.fill_density_buffer();
    sim.fill_potential_buffer();

    double w_c = width / ncells;

    // save to txt file
    // std::vector<double> pot_store; //commented out usually as just extra sanity check
    // std::vector<double> expected_pot_store;

    // Look at the potential function along the x axis
    for (int i = 0; i < 101; i++)
    {
        int j = 50;
        int k = 50;
        if (i == 50 && j == 50 && k == 50)
        {
            //ignore the centre where point potential would be singular
            continue;
        }
        else
        {
            // ignore imaginary component
            // approximate potential for multiple sources due to periodicity
            double dx1 = std::abs(i - 50) * w_c;
            double dx2 = std::abs(i - 151) * w_c;
            double dx3 = std::abs(i + 51) * w_c;
            double pot = sim.get_potential_buffer()[k + ncells * (j + ncells * i)][0]; //TODO = your potential function at indices (i,j,k)
            double expected_pot =  -mass *(1/ dx1 + 1/dx2 + 1/dx3);
            double diff = pot - expected_pot;
            REQUIRE_THAT(pot, WithinRel(expected_pot, 0.3));
            // pot_store.push_back(pot);
            // expected_pot_store.push_back(expected_pot);
        }
    }
    // std::string filename = "test_potential/test.txt";
    // PotentialSavetoTxt(pot_store, expected_pot_store, filename);
}

TEST_CASE("Test gradient function for periodic f(x) = sin(x) + cos(y) + sin(z)", "[Gradient_Function]"){
    uint num_cells = 100;
    double width = 2*M_PI;

    uint buffer_length = num_cells * num_cells * num_cells;
    fftw_complex * test_func_buffer = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * buffer_length);
    std::memset(test_func_buffer, 0, sizeof(fftw_complex) * buffer_length);
    std::vector<double> func(buffer_length, 0);
    std::vector<std::vector<std::vector<std::array<double, 3>>>> test_grad;
    
    for (uint i = 0; i < num_cells; i++){
        std::vector<std::vector<std::array<double, 3>>> test_grad_jk;
        for (uint j = 0; j < num_cells; j++){
            std::vector<std::array<double, 3>> test_grad_k;
            for (uint k = 0; k < num_cells; k++){
                double x = width * (i + 0.5)/num_cells;
                double y = width * (j + 0.5)/num_cells;
                double z = width * (k + 0.5)/num_cells;
                test_func_buffer[k + num_cells * (j + num_cells * i)][0] = std::sin(x) + std::cos(y) + std::sin(z);
                std::array<double, 3> test_grad_section;
                test_grad_section[0] = std::cos(x);
                test_grad_section[1] = - std::sin(y);
                test_grad_section[2] =  std::cos(z);
                
                test_grad_k.push_back(test_grad_section);
            }
            test_grad_jk.push_back(test_grad_k);
        }
        test_grad.push_back(test_grad_jk);
    }
    particle_group particles(1, 1, {{0.5, 0.5, 0.5}});
    Simulation sim(10, 1, particles, width, num_cells, 3);

    std::vector<std::vector<std::vector<std::array<double, 3>>>> grad = sim.calculate_gradient(test_func_buffer);
    
    for (uint i = 0; i < num_cells; i++){
        for (uint j = 0; j < num_cells; j++){
            for (uint k = 0; k < num_cells; k++){
                REQUIRE_THAT(grad[i][j][k][0], WithinRel(test_grad[i][j][k][0], 1e-3));
                REQUIRE_THAT(grad[i][j][k][1], WithinRel(test_grad[i][j][k][1], 1e-3));
                REQUIRE_THAT(grad[i][j][k][2], WithinRel(test_grad[i][j][k][2], 1e-3));
            }
        }
    }
    
    fftw_free(test_func_buffer);
}


TEST_CASE("Test gradient function for periodic f(x) = sin(x) + sin(y) + sin(z)", "[Gradient_Function]"){
    uint num_cells = 100;
    double width = 2*M_PI;

    uint buffer_length = num_cells * num_cells * num_cells;
    fftw_complex * test_func_buffer = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * buffer_length);
    std::memset(test_func_buffer, 0, sizeof(fftw_complex) * buffer_length);
    std::vector<double> func(buffer_length, 0);
    std::vector<std::vector<std::vector<std::array<double, 3>>>> test_grad;
    
    for (uint i = 0; i < num_cells; i++){
        std::vector<std::vector<std::array<double, 3>>> test_grad_jk;
        for (uint j = 0; j < num_cells; j++){
            std::vector<std::array<double, 3>> test_grad_k;
            for (uint k = 0; k < num_cells; k++){
                double x = width * (i + 0.5)/num_cells;
                double y = width * (j + 0.5)/num_cells;
                double z = width * (k + 0.5)/num_cells;
                test_func_buffer[k + num_cells * (j + num_cells * i)][0] = std::sin(x) + std::sin(y) + std::sin(z);
                std::array<double, 3> test_grad_section;
                test_grad_section[0] = std::cos(x);
                test_grad_section[1] = std::cos(y);
                test_grad_section[2] =  std::cos(z);
                
                test_grad_k.push_back(test_grad_section);
            }
            test_grad_jk.push_back(test_grad_k);
        }
        test_grad.push_back(test_grad_jk);
    }
    particle_group particles(1, 1, {{0.5, 0.5, 0.5}});
    Simulation sim(10, 1, particles, width, num_cells, 3);

    std::vector<std::vector<std::vector<std::array<double, 3>>>> grad = sim.calculate_gradient(test_func_buffer);
    
    for (uint i = 0; i < num_cells; i++){
        for (uint j = 0; j < num_cells; j++){
            for (uint k = 0; k < num_cells; k++){
                REQUIRE_THAT(grad[i][j][k][0], WithinRel(test_grad[i][j][k][0], 1e-3));
                REQUIRE_THAT(grad[i][j][k][1], WithinRel(test_grad[i][j][k][1], 1e-3));
                REQUIRE_THAT(grad[i][j][k][2], WithinRel(test_grad[i][j][k][2], 1e-3));
            }
        }
    }
    
    fftw_free(test_func_buffer);
}


TEST_CASE("Test gradient function for periodic f(x) = cos^2(x) + sin^2(y) + cos(x)", "[Gradient_Function]"){
    uint num_cells = 100;
    double width =2 * M_PI;

    uint buffer_length = num_cells * num_cells * num_cells;
    fftw_complex * test_func_buffer = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * buffer_length);
    std::memset(test_func_buffer, 0, sizeof(fftw_complex) * buffer_length);
    std::vector<double> func(buffer_length, 0);
    std::vector<std::vector<std::vector<std::array<double, 3>>>> test_grad;
    
    for (uint i = 0; i < num_cells; i++){
        std::vector<std::vector<std::array<double, 3>>> test_grad_jk;
        for (uint j = 0; j < num_cells; j++){
            std::vector<std::array<double, 3>> test_grad_k;
            for (uint k = 0; k < num_cells; k++){
                double x = width * (i + 0.5)/num_cells;
                double y = width * (j + 0.5)/num_cells;
                double z = width * (k + 0.5)/num_cells;
                test_func_buffer[k + num_cells * (j + num_cells * i)][0] = std::cos(x) * std::cos(x) + std::sin(y) * std::sin(y) + std::cos(z);
                std::array<double, 3> test_grad_section;
                test_grad_section[0] = -2 * std::sin(x) * std::cos(x);
                test_grad_section[1] = 2 * std::sin(y) * std::cos(y);
                test_grad_section[2] =  -std::sin(z);
                
                test_grad_k.push_back(test_grad_section);
            }
            test_grad_jk.push_back(test_grad_k);
        }
        test_grad.push_back(test_grad_jk);
    }
    particle_group particles(1, 1, {{0.5, 0.5, 0.5}});
    Simulation sim(10, 1, particles, width, num_cells, 3);

    std::vector<std::vector<std::vector<std::array<double, 3>>>> grad = sim.calculate_gradient(test_func_buffer);
    
    for (uint i = 0; i < num_cells; i++){
        for (uint j = 0; j < num_cells; j++){
            for (uint k = 0; k < num_cells; k++){
                REQUIRE_THAT(grad[i][j][k][0], WithinRel(test_grad[i][j][k][0], 0.3));
                REQUIRE_THAT(grad[i][j][k][1], WithinRel(test_grad[i][j][k][1], 0.3));
                REQUIRE_THAT(grad[i][j][k][2], WithinRel(test_grad[i][j][k][2], 0.3));
            }
        }
    }
    
    fftw_free(test_func_buffer);
}

TEST_CASE("Ensure two particles approach each other:","[Update_Particle]"){
    double mass = 0.1;
    double width = 1;
    uint number_particles = 2;
    particle_group particles(mass, number_particles, {{0.4, 0.4, 0.4}, {0.8, 0.8, 0.8}});
    uint num_cells = 100;
    double cell_width = width/num_cells;
    Simulation sim(10, 0.1, particles, width, num_cells, 2);
    
    double prev_distance = std::sqrt(0.4 * 0.4 * 3);

    for (uint i = 0; i < 10; i++){
        sim.fill_density_buffer();
        sim.fill_potential_buffer();
        sim.update_particles();
        const particle_group particle_collection = sim.get_particle_collection(); // testing particles approaching each other
        
        // given limited timesteps the particles must approach each other. Calculating total distance
        double new_distance = std::pow((particle_collection.particles[0].position[0] - particle_collection.particles[1].position[0]), 2);
        new_distance += std::pow((particle_collection.particles[0].position[1] - particle_collection.particles[1].position[1]), 2);
        new_distance += std::pow((particle_collection.particles[0].position[2] - particle_collection.particles[1].position[2]), 2);
        new_distance = std::sqrt(new_distance);
        REQUIRE(new_distance < prev_distance);
        prev_distance = new_distance;
    }
}

TEST_CASE("Ensure displacement is centered around 0", "[Update_Particle]"){
    double mass = 0.1;
    double width = 1;
    uint number_particles = 2;
    particle_group particles(mass, number_particles, {{0.30, 0.30, 0.30}, {0.70, 0.70, 0.70}});
    uint num_cells = 10;
    double cell_width = width/num_cells;
    Simulation sim(10, 0.01, particles, width, num_cells, 2);
    
    std::vector<double> distances_x = {0.4};
    std::vector<double> distances_y = {0.4};
    std::vector<double> distances_z = {0.4};

    std::vector<double> velocities_x = {0};
    std::vector<double> velocities_y = {0};
    std::vector<double> velocities_z = {0};

    for (uint i = 0; i < 200000; i++){
        sim.fill_density_buffer();
        sim.fill_potential_buffer();
        sim.update_particles();
        const particle_group particle_collection = sim.get_particle_collection(); // testing particles approaching each other
        
        // given limited timesteps the particles must approach each other. Calculating total distance

        double distance_x = (particle_collection.particles[0].position[0] - particle_collection.particles[1].position[0]);
        double distance_y = (particle_collection.particles[0].position[1] - particle_collection.particles[1].position[1]);
        double distance_z = (particle_collection.particles[0].position[2] - particle_collection.particles[1].position[2]);

        distances_x.push_back(distance_x);
        distances_y.push_back(distance_y);
        distances_z.push_back(distance_z);

        double velocity_x = (particle_collection.particles[0].velocity[0] - particle_collection.particles[1].velocity[0]);
        double velocity_y = (particle_collection.particles[0].velocity[1] - particle_collection.particles[1].velocity[1]);
        double velocity_z = (particle_collection.particles[0].velocity[2] - particle_collection.particles[1].velocity[2]);

        velocities_x.push_back(velocity_x);
        velocities_y.push_back(velocity_y);
        velocities_z.push_back(velocity_z);
    }

    double mean_dist_x = std::accumulate(distances_x.begin(), distances_x.end(), 0.0)/distances_x.size();
    double mean_dist_y = std::accumulate(distances_y.begin(), distances_y.end(), 0.0)/distances_y.size();
    double mean_dist_z = std::accumulate(distances_z.begin(), distances_z.end(), 0.0)/distances_z.size();

    double mean_vel_x = std::accumulate(velocities_x.begin(), velocities_x.end(), 0.0)/velocities_x.size();
    double mean_vel_y = std::accumulate(velocities_y.begin(), velocities_y.end(), 0.0)/velocities_y.size();
    double mean_vel_z = std::accumulate(velocities_z.begin(), velocities_z.end(), 0.0)/velocities_z.size();

    REQUIRE_THAT(mean_dist_x, WithinAbs(0,0.1)); // use hypothesis that displacmement and difference of velocity between particles is 0 centered                                      
    REQUIRE_THAT(mean_dist_y, WithinAbs(0,0.1)); 
    REQUIRE_THAT(mean_dist_z, WithinAbs(0,0.1));


    REQUIRE_THAT(mean_vel_x, WithinAbs(0,0.05));
    REQUIRE_THAT(mean_vel_y, WithinAbs(0,0.05));
    REQUIRE_THAT(mean_vel_z, WithinAbs(0,0.05));

    // std::string filename("test_potential/trajectories.txt");
    // TrajectorySavetoTxt(distances_x, distances_y, distances_z, velocities_x, velocities_y, velocities_z,filename); //for analysis purposes
}

TEST_CASE("Ensure particles remain stationary when distance between them through boundaries and directly is the same","[Update_Particle]"){
    double mass = 0.1;
    double width = 1;
    uint number_particles = 2;
    particle_group particles(mass, number_particles, {{0.25, 0.25, 0.25}, {0.75, 0.75, 0.75}});
    uint num_cells = 10;
    double cell_width = width/num_cells;
    Simulation sim(10, 0.01, particles, width, num_cells, 2);
    


    for (uint i = 0; i < 2000; i++){
        sim.fill_density_buffer();
        sim.fill_potential_buffer();
        sim.update_particles();
        const particle_group particle_collection = sim.get_particle_collection(); 
        
        // particles must remain stationary as gravitational force in either direction is the same. Total field at particle coordinates is uniform

        REQUIRE_THAT(particle_collection.particles[0].position[0], WithinRel(0.25,1e-6));
        REQUIRE_THAT(particle_collection.particles[0].position[1], WithinRel(0.25,1e-6));
        REQUIRE_THAT(particle_collection.particles[0].position[2], WithinRel(0.25,1e-6));
        REQUIRE_THAT(particle_collection.particles[1].position[0], WithinRel(0.75,1e-6));
        REQUIRE_THAT(particle_collection.particles[1].position[1], WithinRel(0.75,1e-6));
        REQUIRE_THAT(particle_collection.particles[1].position[2], WithinRel(0.75,1e-6));
        
        REQUIRE_THAT(particle_collection.particles[0].velocity[0], WithinAbs(0,1e-6));
        REQUIRE_THAT(particle_collection.particles[0].velocity[1], WithinAbs(0,1e-6));
        REQUIRE_THAT(particle_collection.particles[0].velocity[2], WithinAbs(0,1e-6));
        REQUIRE_THAT(particle_collection.particles[1].velocity[0], WithinAbs(0,1e-6));
        REQUIRE_THAT(particle_collection.particles[1].velocity[1], WithinAbs(0,1e-6));
        REQUIRE_THAT(particle_collection.particles[1].velocity[2], WithinAbs(0,1e-6));
    }

}