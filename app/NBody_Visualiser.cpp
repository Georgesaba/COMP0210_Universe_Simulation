#include "Simulation.hpp"
#include <iostream>
#include <optional>
#include <filesystem>
#include <memory>
#include "Utils.hpp"

/**
 * @brief: This function prints a help message for the NBody_Visualiser application
*/
void HelpMessage(){
    std::cout << "This program visualises a developing universe through modelling the graviational fields of multiple particles with the same mass using the particle mesh method.\n\nBrief instructions can be found below." << std::endl;
    std::cout << "Usage: NBody_Visualiser -nc <number_of_cells> -np <average_particles_per_cell> -t <total_time> -dt <time_step> -F <expansion_factor> -o <output_folder> -s <random_seed>\n"
              << "Options:\n"
              << "  -h                                       Show this help message\n"
              << "  -nc <number_of_cells>                    Number of cells wide the equal sided box has\n"
              << "  -np <average_particles_per_cell>         Average number of particles each cell has\n"
              << "  -t  <total_time>                         Total time the simulation will run for\n"
              << "  -dt <time_step>                          Amount of time that is incremented each propagation\n"
              << "  -F  <expansion_factor>                   Factor that the absolute value of the box expands\n"
              << "  -o  <output_folder>                      Folder that output images are sent to\n"
              << "  -s  <random_seed>                        Seed that is used to generate initial randomised positions" << std::endl;
}

int main(int argc, char** argv)
{
    
    std::string output_folder;
    uint num_cells;
    uint random_seed;
    double average_particles_per_cell;
    double time_step;
    double expansion_factor;
    double max_time;

    bool output_folder_set = false;
    bool num_cells_set = false;
    bool average_particle_per_cell_set = false; // track if values have been set already.
    bool time_step_set = false;
    bool expansion_factor_set = false;
    bool random_seed_set = false;
    bool max_time_set = false;
    
    for (uint i = 1; i < argc; i+=2){
        std::string arg(argv[i]);
        if (arg == "-h"){
            HelpMessage();
            return 0;
        }
        else if (arg == "-o") {
            if (output_folder_set){
                std::cerr << "Error - Cannot set output folder twice!" << std::endl; // check for clashes caused by setting the same value to different conditions. 
                HelpMessage();
                return 1;
            }
            std::string arg1(argv[i+1]);
            output_folder = arg1;
            output_folder_set = true;
        }
        else if (arg == "-nc"){
            if (num_cells_set){
                std::cerr << "Error - Cannot set the number of cells twice!" << std::endl;
                HelpMessage();
                return 1;
            }
            std::string arg1(argv[i+1]);
            num_cells = std::atoi(arg1.c_str());
            if (num_cells > 220){
                std::cerr << "Warning - Process may be killed as the number of cells exceeds 220! Reduce the -np or -nc settings if this happens!" << std::endl;
            }
            num_cells_set = true;
        }
        else if (arg == "-np"){
            if (average_particle_per_cell_set){
                std::cerr << "Error - Cannot set the number of particles per cell twice!" << std::endl;
                HelpMessage();
                return 1;
            }
            std::string arg1(argv[i+1]);
            average_particles_per_cell = std::stod(arg1.c_str()); // double so more flexibility is possible. Performance impact is negligible
            average_particle_per_cell_set = true;
        } 
        else if (arg == "-t"){ 
            if (max_time_set){
                std::cerr << "Error - the maximum time has already been set!" << std::endl;
                HelpMessage();
                return 1;
            }
            std::string arg1(argv[i+1]);
            max_time = std::stod(arg1.c_str());
            max_time_set = true;
        }
        else if (arg == "-dt"){
            if (time_step_set){
                std::cerr << "Error - the timestep has already been set!" << std::endl;
                HelpMessage();
                return 1;
            }
            std::string arg1(argv[i+1]);
            time_step = std::stod(arg1.c_str());
            time_step_set = true;
        }
        else if (arg == "-F"){
            if (expansion_factor_set){
                std::cerr << "Error - the expansion factor has already been set!" << std::endl;
                HelpMessage();
                return 1;
            }
            std::string arg1(argv[i + 1]);
            expansion_factor = std::stod(arg1.c_str());
            expansion_factor_set = true;
        }
        else if (arg == "-s"){
            if (random_seed_set){
                std::cerr << "Error - the expansion factor has already been set!" << std::endl;
                HelpMessage();
                return 1;
            }
            std::string arg1(argv[i + 1]);
            random_seed = std::atoi(arg1.c_str());
            random_seed_set = true;
        }
        else{ // extra error handling
            std::cerr << "Invalid Flag Detected: " << arg << std::endl;
            HelpMessage();
            return 1;
        }
    }
    
    if (!(output_folder_set && num_cells_set && average_particle_per_cell_set && time_step_set && expansion_factor_set && random_seed_set && max_time_set)){
        std::cerr << "Please Input the Required Flags!" << std::endl;
        HelpMessage();
        return 1;
    }

    double width = 100.0;
    uint num_particles = num_cells * num_cells * num_cells * average_particles_per_cell;
    double mass = 10.0 * 10.0 * 10.0 * 10.0 * 10.0/num_particles;
    
    std::unique_ptr<Simulation> Simulation_ptr;

    try{
        particle_group particles(mass, num_particles, random_seed);
        Simulation_ptr = std::make_unique<Simulation>(max_time, time_step, particles, width, num_cells, expansion_factor);
    }
    catch (const std::bad_alloc &e){
        std::cerr << "Error - Memory Overflow: Please use smaller values for -nc <number_of_cells> or -np <average_number_particles_per_cell> arguments!" << std::endl;
        HelpMessage();
        return 1;
    }
    catch(const std::exception &e){
        std::cerr << e.what() << std::endl;
        HelpMessage();
        return 1;
    }
    output_folder += "/" +  removeTrailingDecimalPlaces(random_seed);
    Simulation_ptr->run(output_folder);
    
    return 0;
}