# Responses

## Section 1.2

In `include/particle.hpp` and `lib/particle.cp`p I have outlined two classes. `particle` is responsible for defining characteristics of a single particle, and `particle_group` is reponsible for describing a collection of particles. The `particle` class has a size 3 velocity std::array member variable which is set to all zero values and a size 3 position std::array member variable. Its constructor sets an initial position for it. All member variables are public so that individual coordinates can be changed independently without having to replace the whole position or velocity array every time it is accessed. The `particle_group` class has particle mass and the number of particles stored as member variables as well as an std::vector of particle objects to store all the particle specific particle information. I have chosen this format so that I can store the features that every particle has only once and then bundle the features each particle has in a single object for ease of use. The std::vector is used to store the particle objects as it stores data on the heap making it suitable for extremely large numbers of particles. This also has performance benefits as all the qualities of a particle can be accessed with a single memory access of a vector instead of indexing vectors multiple times. This is especially true of larger vectors where the data in question may be stored outside of the cache. I decided to make constructors that work with different argument combinations by taking advantage of constructor overloading. The first constructor allows for the position of each particle to be set individually through passing an std::vector of length 3 std::arrays containing the x, y and z coordinates of the particle and the second constructor gives the particles random positions between 0 and 1 for each x, y and z particle coordinates. Both constructors assign the particle mass and the number of particles for the collection of particles.


## Section 1.10 - Complexity

#### Density Calculation:

The density calculation has time complexity $O(N_c) + O(n_p)$ as it used the memset operation from the standard library which initialises the density_buffer to all zero values due to the fact that it accesses each element of an array of size $N_c$. The main loop that follows iterates through every particle and its complexity is determined by the number of particles $n_p$ as the other indexing operations are constant time so do not contribute due to having $O(1)$ complexity.

#### Potential Calculation:

The fast fourier transform has time complexity $O(N_c log(N_c))$. Every element in the buffer is then iterated over individually making it have time complexity of $O(N_c)$. The inverse fast fourier transform has time complexity $O(N_c log(N_c))$ as well. Therefore the effective total time complexity of the operation is $O(N_c log(N_c)) + O(N_c)$ however as $N_c$ gets larger the time complexity reduces to $O(N_c log(N_c))$.

#### Gradient of Potential Calculation:

The initialisation of the three dimensional vector of length three arrays all filled with zeros has a time complexity of $O(N_c)$ as each element in the $n_c \times n_c \times n_c$ size vector is accessed. The nested for loop also has time complexity $O(N_c)$ as the partial derivative with respect to $x$, $y$ and $z$ is calculated $N_c$ times with each iteration performing a constant amount of work. The total time complexity is therefore $O(N_c)$.

#### Particle Update Function

Ignoring the usage of the gradient function in the particle update function, the time compexity of the particle update function is $O(n_p)$ as the operations in each loop have constant time complexity assuming that the boundary conditions are exceeded a limited number of times in one go and the algorithm loops $n_p$ times making it linearly dependent on $n_p$.

#### Conclusion

In order to achieve a smooth density function the number of particles must always be larger than the number of particles in the mesh. The total time complexity of the mesh method is $O(n_p) + O(N_c) + O(N_c log(N_c)) + O(N_c) + O(N_c) + O(n_p)$ which can be reduced to $O(n_p) + O(N_clog(N_c))$ which for larger simulations with large numbers of particles is less than $O(n_p^2)$. Therefore the Particle-Mesh method will scale much better than the Particle-Particle method for these simulations.


## Section 1.11: Running Simulation

File storage format is `Images/seed_num/expansion_factor` with the other parameters included in the labelling. 4 different expansion factors were used. 0.98, 1, 1.02 and 1.1. A dt value of 0.01 with a max time of 1.5 was used and the average particles per cell value that was used was 12. The fixed random seed was always 42. The number of cells wide the box was is 101 although 201 was also tested for more granular results.

## Section 2.1: Identifying Simulation functions for Parallelisation

In the Simulation class the member functions fill_density_buffer(), fill_potential_buffer(), calculate_gradient(fftw_complex * potential), update_particles() and box_expansion() can all be parallelised as the results of operations in the loop do not depend on other processes. The `#pragma omp atomic` was required in the fill_density_buffer as the particles were being handled by separate threads although the density array that was being incremented had each index corresponding to a cell in the cube grid which meant that this could have caused a data race. The atomic command was used to keep the process atomic and prevent this.

#### Function Timings

##### Run 1

Benchmarking Density Calculation with 1 threads.
Time = 0.407183
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Density Calculation with 2 threads.
Time = 0.202571
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Density Calculation with 3 threads.
Time = 0.143954
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Density Calculation with 4 threads.
Time = 0.101066
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Density Calculation with 5 threads.
Time = 0.0843676
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Density Calculation with 6 threads.
Time = 0.110924
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Density Calculation with 7 threads.
Time = 0.0895665
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Density Calculation with 8 threads.
Time = 0.0877549
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Density Calculation with 9 threads.
Time = 0.0558973
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Density Calculation with 10 threads.
Time = 0.0602805
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Density Calculation with 11 threads.
Time = 0.055638
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Density Calculation with 12 threads.
Time = 0.0501531
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Density Calculation with 13 threads.
Time = 0.0488205
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Density Calculation with 14 threads.
Time = 0.0414185
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Density Calculation with 15 threads.
Time = 0.0496775
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Density Calculation with 16 threads.
Time = 0.0491985
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Potential Calculation with 1 threads.
Time = 0.086307
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Potential Calculation with 2 threads.
Time = 0.105839
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Potential Calculation with 3 threads.
Time = 0.0830604
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Potential Calculation with 4 threads.
Time = 0.0796004
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Potential Calculation with 5 threads.
Time = 0.0808205
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Potential Calculation with 6 threads.
Time = 0.0825749
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Potential Calculation with 7 threads.
Time = 0.0821773
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Potential Calculation with 8 threads.
Time = 0.0827839
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Potential Calculation with 9 threads.
Time = 0.0917558
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Potential Calculation with 10 threads.
Time = 0.08282
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Potential Calculation with 11 threads.
Time = 0.0803666
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Potential Calculation with 12 threads.
Time = 0.0813482
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Potential Calculation with 13 threads.
Time = 0.0832948
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Potential Calculation with 14 threads.
Time = 0.0812494
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Potential Calculation with 15 threads.
Time = 0.0853361
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Potential Calculation with 16 threads.
Time = 0.0897332
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Gradient Calc with 1 threads.
Time = 0.0095149
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Gradient Calc with 2 threads.
Time = 0.00515912
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Gradient Calc with 3 threads.
Time = 0.00424612
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Gradient Calc with 4 threads.
Time = 0.00379641
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Gradient Calc with 5 threads.
Time = 0.00372309
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Gradient Calc with 6 threads.
Time = 0.00362644
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Gradient Calc with 7 threads.
Time = 0.00380239
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Gradient Calc with 8 threads.
Time = 0.00397423
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Gradient Calc with 9 threads.
Time = 0.00410529
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Gradient Calc with 10 threads.
Time = 0.00387656
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Gradient Calc with 11 threads.
Time = 0.00391536
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Gradient Calc with 12 threads.
Time = 0.00410103
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Gradient Calc with 13 threads.
Time = 0.00377883
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Gradient Calc with 14 threads.
Time = 0.00391388
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Gradient Calc with 15 threads.
Time = 0.00395912
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Gradient Calc with 16 threads.
Time = 0.00715385
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Particle Update and Gradient Calc with 1 threads.
Time = 0.716457
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Particle Update and Gradient Calc with 2 threads.
Time = 0.321261
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Particle Update and Gradient Calc with 3 threads.
Time = 0.241704
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Particle Update and Gradient Calc with 4 threads.
Time = 0.157543
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Particle Update and Gradient Calc with 5 threads.
Time = 0.160937
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Particle Update and Gradient Calc with 6 threads.
Time = 0.150223
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Particle Update and Gradient Calc with 7 threads.
Time = 0.147559
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Particle Update and Gradient Calc with 8 threads.
Time = 0.132156
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Particle Update and Gradient Calc with 9 threads.
Time = 0.118514
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Particle Update and Gradient Calc with 10 threads.
Time = 0.10748
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Particle Update and Gradient Calc with 11 threads.
Time = 0.0990091
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Particle Update and Gradient Calc with 12 threads.
Time = 0.0880951
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Particle Update and Gradient Calc with 13 threads.
Time = 0.0919515
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Particle Update and Gradient Calc with 14 threads.
Time = 0.0736343
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Particle Update and Gradient Calc with 15 threads.
Time = 0.0885005
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Particle Update and Gradient Calc with 16 threads.
Time = 0.0867701
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Expansion Calculation with 1 threads.
Time = 0.0382308
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Expansion Calculation with 2 threads.
Time = 0.0211765
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Expansion Calculation with 3 threads.
Time = 0.0234968
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Expansion Calculation with 4 threads.
Time = 0.0221637
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Expansion Calculation with 5 threads.
Time = 0.0238748
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Expansion Calculation with 6 threads.
Time = 0.022937
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Expansion Calculation with 7 threads.
Time = 0.0248403
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Expansion Calculation with 8 threads.
Time = 0.0246648
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Expansion Calculation with 9 threads.
Time = 0.0232902
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Expansion Calculation with 10 threads.
Time = 0.024
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Expansion Calculation with 11 threads.
Time = 0.0240411
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Expansion Calculation with 12 threads.
Time = 0.0238106
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Expansion Calculation with 13 threads.
Time = 0.0237971
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Expansion Calculation with 14 threads.
Time = 0.0236321
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Expansion Calculation with 15 threads.
Time = 0.0240168
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

Benchmarking Expansion Calculation with 16 threads.
Time = 0.0241153
Info: The number of cells per length of the box is 101 and the number of particles is 12363612.

##### Run 2

Benchmarking Density Calculation with 1 threads.
Time = 1.72939
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Density Calculation with 2 threads.
Time = 1.13353
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Density Calculation with 3 threads.
Time = 0.712041
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Density Calculation with 4 threads.
Time = 0.66141
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Density Calculation with 5 threads.
Time = 0.520107
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Density Calculation with 6 threads.
Time = 0.467166
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Density Calculation with 7 threads.
Time = 0.417774
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Density Calculation with 8 threads.
Time = 0.38871
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Density Calculation with 9 threads.
Time = 0.408177
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Density Calculation with 10 threads.
Time = 0.408422
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Density Calculation with 11 threads.
Time = 0.359446
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Density Calculation with 12 threads.
Time = 0.401336
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Density Calculation with 13 threads.
Time = 0.387255
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Density Calculation with 14 threads.
Time = 0.431734
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Density Calculation with 15 threads.
Time = 0.406753
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Density Calculation with 16 threads.
Time = 0.391926
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Potential Calculation with 1 threads.
Time = 0.372276
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Potential Calculation with 2 threads.
Time = 0.369512
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Potential Calculation with 3 threads.
Time = 0.350012
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Potential Calculation with 4 threads.
Time = 0.358653
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Potential Calculation with 5 threads.
Time = 0.355189
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Potential Calculation with 6 threads.
Time = 0.339495
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Potential Calculation with 7 threads.
Time = 0.343342
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Potential Calculation with 8 threads.
Time = 0.334725
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Potential Calculation with 9 threads.
Time = 0.341988
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Potential Calculation with 10 threads.
Time = 0.359727
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Potential Calculation with 11 threads.
Time = 0.336278
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Potential Calculation with 12 threads.
Time = 0.3428
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Potential Calculation with 13 threads.
Time = 0.335189
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Potential Calculation with 14 threads.
Time = 0.341599
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Potential Calculation with 15 threads.
Time = 0.342395
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Potential Calculation with 16 threads.
Time = 0.339855
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Gradient Calc with 1 threads.
Time = 0.0619885
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Gradient Calc with 2 threads.
Time = 0.0379624
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Gradient Calc with 3 threads.
Time = 0.0329366
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Gradient Calc with 4 threads.
Time = 0.0353653
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Gradient Calc with 5 threads.
Time = 0.0330857
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Gradient Calc with 6 threads.
Time = 0.0315926
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Gradient Calc with 7 threads.
Time = 0.0300559
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Gradient Calc with 8 threads.
Time = 0.0314777
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Gradient Calc with 9 threads.
Time = 0.029237
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Gradient Calc with 10 threads.
Time = 0.0385783
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Gradient Calc with 11 threads.
Time = 0.0332832
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Gradient Calc with 12 threads.
Time = 0.0298546
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Gradient Calc with 13 threads.
Time = 0.0302792
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Gradient Calc with 14 threads.
Time = 0.0304599
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Gradient Calc with 15 threads.
Time = 0.0317404
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Gradient Calc with 16 threads.
Time = 0.0319498
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Particle Update and Gradient Calc with 1 threads.
Time = 3.26482
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Particle Update and Gradient Calc with 2 threads.
Time = 1.75527
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Particle Update and Gradient Calc with 3 threads.
Time = 1.27501
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Particle Update and Gradient Calc with 4 threads.
Time = 0.97757
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Particle Update and Gradient Calc with 5 threads.
Time = 0.877254
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Particle Update and Gradient Calc with 6 threads.
Time = 0.759223
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Particle Update and Gradient Calc with 7 threads.
Time = 0.692309
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Particle Update and Gradient Calc with 8 threads.
Time = 0.604826
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Particle Update and Gradient Calc with 9 threads.
Time = 0.566183
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Particle Update and Gradient Calc with 10 threads.
Time = 0.639536
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Particle Update and Gradient Calc with 11 threads.
Time = 0.660705
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Particle Update and Gradient Calc with 12 threads.
Time = 0.467521
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Particle Update and Gradient Calc with 13 threads.
Time = 0.447448
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Particle Update and Gradient Calc with 14 threads.
Time = 0.438275
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Particle Update and Gradient Calc with 15 threads.
Time = 0.457765
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Particle Update and Gradient Calc with 16 threads.
Time = 0.573821
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Expansion Calculation with 1 threads.
Time = 0.162249
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Expansion Calculation with 2 threads.
Time = 0.0824112
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Expansion Calculation with 3 threads.
Time = 0.0802064
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Expansion Calculation with 4 threads.
Time = 0.0787748
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Expansion Calculation with 5 threads.
Time = 0.0762822
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Expansion Calculation with 6 threads.
Time = 0.0786626
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Expansion Calculation with 7 threads.
Time = 0.0795726
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Expansion Calculation with 8 threads.
Time = 0.0779164
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Expansion Calculation with 9 threads.
Time = 0.0778634
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Expansion Calculation with 10 threads.
Time = 0.0789523
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Expansion Calculation with 11 threads.
Time = 0.0795112
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Expansion Calculation with 12 threads.
Time = 0.0792901
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Expansion Calculation with 13 threads.
Time = 0.0795983
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Expansion Calculation with 14 threads.
Time = 0.0803269
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Expansion Calculation with 15 threads.
Time = 0.0805623
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

Benchmarking Expansion Calculation with 16 threads.
Time = 0.0798672
Info: The number of cells per length of the box is 151 and the number of particles is 41315412.

##### Run 3

Benchmarking Density Calculation with 1 threads.
Time = 7.1269
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Density Calculation with 2 threads.
Time = 1.92785
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Density Calculation with 3 threads.
Time = 1.44682
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Density Calculation with 4 threads.
Time = 1.02325
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Density Calculation with 5 threads.
Time = 0.994784
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Density Calculation with 6 threads.
Time = 0.890907
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Density Calculation with 7 threads.
Time = 0.77327
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Density Calculation with 8 threads.
Time = 0.70953
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Density Calculation with 9 threads.
Time = 0.65449
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Density Calculation with 10 threads.
Time = 0.673552
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Density Calculation with 11 threads.
Time = 0.592343
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Density Calculation with 12 threads.
Time = 0.618197
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Density Calculation with 13 threads.
Time = 0.554707
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Density Calculation with 14 threads.
Time = 0.529864
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Density Calculation with 15 threads.
Time = 0.543532
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Density Calculation with 16 threads.
Time = 0.563541
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Potential Calculation with 1 threads.
Time = 0.96153
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Potential Calculation with 2 threads.
Time = 0.769982
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Potential Calculation with 3 threads.
Time = 0.755752
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Potential Calculation with 4 threads.
Time = 0.792984
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Potential Calculation with 5 threads.
Time = 0.790276
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Potential Calculation with 6 threads.
Time = 0.76409
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Potential Calculation with 7 threads.
Time = 0.781921
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Potential Calculation with 8 threads.
Time = 0.754889
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Potential Calculation with 9 threads.
Time = 0.79205
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Potential Calculation with 10 threads.
Time = 0.809308
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Potential Calculation with 11 threads.
Time = 0.786385
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Potential Calculation with 12 threads.
Time = 0.790811
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Potential Calculation with 13 threads.
Time = 0.774962
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Potential Calculation with 14 threads.
Time = 0.781886
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Potential Calculation with 15 threads.
Time = 0.799603
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Potential Calculation with 16 threads.
Time = 0.793272
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Gradient Calc with 1 threads.
Time = 0.0919779
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Gradient Calc with 2 threads.
Time = 0.0568788
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Gradient Calc with 3 threads.
Time = 0.0490351
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Gradient Calc with 4 threads.
Time = 0.0482513
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Gradient Calc with 5 threads.
Time = 0.0528761
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Gradient Calc with 6 threads.
Time = 0.0502631
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Gradient Calc with 7 threads.
Time = 0.0524614
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Gradient Calc with 8 threads.
Time = 0.049492
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Gradient Calc with 9 threads.
Time = 0.0525983
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Gradient Calc with 10 threads.
Time = 0.0540007
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Gradient Calc with 11 threads.
Time = 0.0540036
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Gradient Calc with 12 threads.
Time = 0.0539652
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Gradient Calc with 13 threads.
Time = 0.0584742
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Gradient Calc with 14 threads.
Time = 0.054549
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Gradient Calc with 15 threads.
Time = 0.0688302
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Gradient Calc with 16 threads.
Time = 0.0685671
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Particle Update and Gradient Calc with 1 threads.
Time = 6.93234
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Particle Update and Gradient Calc with 2 threads.
Time = 3.49144
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Particle Update and Gradient Calc with 3 threads.
Time = 2.48071
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Particle Update and Gradient Calc with 4 threads.
Time = 1.86549
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Particle Update and Gradient Calc with 5 threads.
Time = 1.50419
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Particle Update and Gradient Calc with 6 threads.
Time = 1.32143
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Particle Update and Gradient Calc with 7 threads.
Time = 1.23134
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Particle Update and Gradient Calc with 8 threads.
Time = 1.14441
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Particle Update and Gradient Calc with 9 threads.
Time = 1.09315
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Particle Update and Gradient Calc with 10 threads.
Time = 0.998685
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Particle Update and Gradient Calc with 11 threads.
Time = 1.01025
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Particle Update and Gradient Calc with 12 threads.
Time = 0.955928
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Particle Update and Gradient Calc with 13 threads.
Time = 0.888577
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Particle Update and Gradient Calc with 14 threads.
Time = 0.91668
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Particle Update and Gradient Calc with 15 threads.
Time = 0.928762
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Particle Update and Gradient Calc with 16 threads.
Time = 0.866622
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Expansion Calculation with 1 threads.
Time = 0.316385
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Expansion Calculation with 2 threads.
Time = 0.155196
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Expansion Calculation with 3 threads.
Time = 0.155455
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Expansion Calculation with 4 threads.
Time = 0.150995
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Expansion Calculation with 5 threads.
Time = 0.151513
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Expansion Calculation with 6 threads.
Time = 0.155578
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Expansion Calculation with 7 threads.
Time = 0.156343
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Expansion Calculation with 8 threads.
Time = 0.156845
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Expansion Calculation with 9 threads.
Time = 0.154952
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Expansion Calculation with 10 threads.
Time = 0.153485
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Expansion Calculation with 11 threads.
Time = 0.155504
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Expansion Calculation with 12 threads.
Time = 0.156205
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Expansion Calculation with 13 threads.
Time = 0.156022
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Expansion Calculation with 14 threads.
Time = 0.154474
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Expansion Calculation with 15 threads.
Time = 0.154125
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

Benchmarking Expansion Calculation with 16 threads.
Time = 0.157658
Info: The number of cells per length of the box is 201 and the number of particles is 81206010.

## Section 2.2: Comparing Universes with Distributed Memory

#### Changes to CorrelationFunction()

I have changed the function CorrelationFunction() in Utils.hpp/cpp to take in a particle_group object instead of a vector of positions. This is to remove the need for an extra loop that loops through all the particles to extract their positons and place them in an `std::vector<std::array<std::double>>` object. The Simulation had the hard coded parameters `num_cells = 101`; `average_particles_per_cell = 13`; `width = 100`; `random_seed = 42`; `t_max = 1.5` and `dt = 0.01`. The total mass remains $10^5$. 4 simulations and 101 bins were used for the correlation function.

When running the application I used three sets of command line arguments:

`minimum_expansion_factor = 1; maximum_expansion_factor = 1.05`<br>
`minimum_expansion_factor = 1; maximum_expansion_factor = 1.04`<br>
`minimum_expansion_factor = 1.02; maximum_expansion_factor = 1.04`<br>
`minimum_expansion_factor = 1.03; maximum_expansion_factor = 1.04`<br>

