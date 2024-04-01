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
