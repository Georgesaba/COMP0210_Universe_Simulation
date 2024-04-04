# UniverseInABox
The intention of this project is to write a simple Particle-Mesh gravitational simulation that allows us to simulate the motion of N bodies. This consists of four applications: `TestSimulation`, `BenchmarkSimulation`, `NBody_Comparison` and `NBody_Visualiser`.

This project is compiled using CMake so compiling requires cmake version 3.16 and C++17 at a minimum.

In the same level in the directory as this README.md file, run `cmake -B build` to configure the project and create the build directory. To compile the programs run `cmake --build build`. Now you should be able to find `TestSimulation`, `BenchmarkSimulation`, `NBody_Comparison` and `NBody_Visualiser` in the `/build/bin/` folders. To run a program type `./build/bin/{program_name}`. `TestSimulation` just contains unit tests for the different functions, classes and algorithms used in this project and `BenchmarkSimulation` contains code to print out benchmark times for different functions using different numbers of threads.

###  NBody_Visualiser

This application is built to generate a cubic grid with periodic boundary conditions of length $n_c$ with total number of cells $N_c$ and then generate an amount of particles $n_p$ determined by a user given average number of particles per cell value. The particle mesh method is then used to calculate a graviation field in this grid given a box length and an expansion factor is applied to the box to simulate an expanding universe. This program can be run using a command similar to the one shown below:

```
./build/bin/NBody_Visualiser -nc 101.0 -np 12 -t 1.5 -dt 0.01 -F 1.02 -o Images -s 42
```

The flag `-h` when used displays a help message shows run instructions an explains the required flags used to run the program. All the flags are required and no optional commands exist. `-s` is the random seed that is used with the std::default_random_engine generator from the STL `<random>` library in c++. `-o` is the output folder that is the images are outputted time. `-F` is the factor by which the box is scaled with, `-dt` is the time-step for each iteration in the simulation and `-t` is the total time elapsed. 

```
This program visualises a developing universe through modelling the graviational fields of multiple particles with the same mass using the particle mesh method.

Brief instructions can be found below.
Usage: NBody_Visualiser -nc <number_of_cells> -np <average_particles_per_cell> -t <total_time> -dt <time_step> -F <expansion_factor> -o <output_folder> -s <random_seed>
Options:
  -h                                       Show this help message
  -nc <number_of_cells>                    Number of cells wide the equal sided box has
  -np <average_particles_per_cell>         Average number of particles each cell has
  -t  <total_time>                         Total time the simulation will run for
  -dt <time_step>                          Amount of time that is incremented each propagation
  -F  <expansion_factor>                   Factor that the absolute value of the box expands
  -o  <output_folder>                      Folder that output images are sent to
  -s  <random_seed>                        Seed that is used to generate initial randomised positions
```

This will then output `.pbm` images to the directory `<output_folder>/<seed>/<Expansion_Factor>/`. It should be noted that all values that are used in naming conventions that are not restricted to integers will that at least a decimal `.` following the number even if it is whole. The file naming convention is `UniverseSim_dt_<time_step>_time_<current_time_simulation>_num_cells_<number_of_cells>_ppc_<average_particles_per_cell>.pbm` where `<current-time_simulation>` is the value of the time at the timestep the image of the particle density distribution was captured at. 

### NBody_Comparison

This application runs $x$ different simulations in parallel using distributed memory and then outputs radial correlation statistics of said simulations to a user specified folder in `.csv` format. Each simulation is ran with a different expansion factor. The user needs to input four arguments: The number of simulations $x$, the output folder that the results are saved to, the maximum expansion factor and the minimum expansion factor. The $x$ simulations are generated with expansion equally spaced expansion factors that range between the maximum and minimum ones specified. The program can be run using the below command format:

```
mpirun -np <number_simulations> ./build/bin/NBody_Comparison -o <output_folder> -emin <minimum_expansion_factor> -emax <maximum_expansion_factor>

mpirun -np 4 ./build/bin/NBody_Comparison -o Correlation -emin 1 -emax 1.04
```
Here `mpirun` is used to distribute the program across the specified nodes in order to commence the parallel computation. The `-np` flag is used to specify the number of parallel proccesses that will be used to run independent simulations. The `-o` flag is used to specify the output folder that the binned radial correlations will be outputted to, the `-emin` flag is used to specify the minimum expansion factor that will be used and `-emax` represents the maximum expansion factor.

The file naming convention of the output `.csv` file is `Comparison_<number_simulations>_<minimum_expansion_factor>_<maximum_expansion_factor>.csv`.
