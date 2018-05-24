---
layout: post
title: Initial Benchmarks of Distance Calculations
---

*Using Jupyter Notebooks for the first time, and they are awesome!!!*

This is the first post after the start of GSoC. A reflection on the past few weeks and summary of the things I have tried, failed and learned should be a good start to this blog.

# Journey so far
---------------

Community bonding period started as soon as the good news arrived, but it ended quickly before I could notice it. Although exams and few other chores also had their fair share in consuming the time. With the remaining time at hand, I focussed on the first part of the project. 
To give a brief introduction, the project focusses on speeding up the distance calculations, which is one of the fundamental requirement in the MDAnalysis package and Molecular Dynamics Analysis packages in general. More information on the project can be found [here](https://summerofcode.withgoogle.com/projects/#5050592943144960).
My first part of the project include 
1. Benchmarking the current methods in MDAnalysis
To evaluate the distances efficiently, multiple methods are already implemented in MDAnalysis. While every method has its own advantages and disadvantages, it is essential to quantify the applicability of these methods. Below is a brief description of these methods and their relevant application in reference to this project. It should be mentioned that the applications mentioned below donot cover the complete breadth of the method.
   1. **Brute Force** :-  This is naive pairwise distance calculation, where distance between every possible pair of particles is calculated individually and stored in a huge matrix. The advantage of the method is its simplicity and easy parallelization. Furthermore, it is easy to implement periodic boundary conditions using the minimum image convention, where the minimum distance among all the neighbouring images is stored in the distance matrix.
   2. **KDTree** :-  KDtree is a space partitioning data structure for organizing points in a k dimensional space. Typically, KDtrees are used for range searches and near neighbour searches in k dimensional space. For the requirements of MDAnalysis, KDtree has been implememented to generate tree structure in 3 dimensions. 
   
   ![alt text](/images/kdtree.gif "KDtree Data structure")
   
   Furthermore, complete tree is usually not built since time for generation of tree structure is high and becomes infeasible for analysis of MD trajectories. A bucket size of 10 is chosen as a default in MDAnalysis package, whereas a larger/smaller size can be used to optimize the calculations depending on the data points. Another major issue which takes up a significant amount of time, specific to Molecular Dynamics simulations is Periodic boundary conditions(PBC). In simple words, an atom crossing one of the edge will reappear from the other side, this is a common boundary condition used in atomic simulations. While no direct treatment is available to handle PBC using KDtree, a common trick is to construct a tree structure after bringing all the atoms to the central cell and then check the distance in all the image cells during search. For a 3D system, this operation directly increases the computation cost by 27. 
   3. **CellGrid** :- A more common approach, particularly used in MD integrators, is to subdivide the domain using a rectangular grid and assigning atoms to cells based on their location. The size of the grid can be different than the cutoff distance. Typically, the MD integrators evaluate force upto a cutoff distance, which is used as the size of the grid. MDAnalysis users typically query for varying search radius around a protein or a particle. Currently, Cellgrid uses the cutoff radius as the grid size (which has its own disadvantages, and we will discuss it later). Once the list is created, a particular search query from a point will only evaluate the distances from the particles residing in the neighbouring cell as well as the particles in the same cell. 
These methods are already in place in MDAnalysis. The task is to analyze the methods and their suitable applications. This will involve analysis of two use cases(range search and contact search) with respect to their running time, memory consumption. The range search involves selecting the group of particles surrounding a particuar particle of a group of particles within a certain cutoff distance. Contant search on the other hand involves looking for all the pair of particles in contact with each other or within a certain cutoff distance from each other. 
2. The cellgrid module in MDAnalysis is not optimized in terms of grid size and cutoff distances. The next task is to optimize the Cellgrid API and Benchmark it against the previous cases. 
3.  Third task involves search of a method which can be utilized for even faster distance calculations and its implementation. This is to ensure that the remaining time of the project is in the correct direction in implementing solution which enhances the capability of MDAnalysis in fast and efficient distance evaluations. Some of the probabe candidates are Octree, Loose Octree, Verlet List, Combination of Cellgrid and Verlet list, Combination of Cell-list and fast distance calculations using SIMD programming. More details will be posted in later blogs as the project continues. 

# Benchmark Studies
--------------------

For all the benchmark studies presented below, I am using Intel Core i5-5200 processor clocked at 2.20 GHz and RAM of 8GB. The related ipython notebooks are available [here](https://github.com/ayushsuhane/Benchmarks_Distance).

1. For the use case to identify contact pairs within a cutoff distance, all the three methods mentioned above are classified based on their execution time. The settings for the benchmarks include a contact distance of 10 units in a uniform distribution of points (which is representative of solvent particles in molecular simulations). All three methods returned the same output. timeit and memit magic methods are used to capture the execution time and memory consumption during the whole step of searching. While Brite force is difficult to beat for small number of particles, the eexecution time increases rapidly for brute force with increase in number of particles. 

![alt text](/images/230518_paircon-bm-timeall.PNG "Execution time for contact pair search")

While Brute force is not used here to benchmark larger systems,  it can be easily seen that cellgrid is clearly advantageous for dense systems. 

![alt text](/images/230518_paircon-bm-time-cgkd.PNG "Execution time for contact pair search  for dense system")

Furthermore, memory consumption using memit shows the rise in memory consumption  in brute force method once the number of particles reaches beyond 4k particles. This sudden increase in memory is not desirable, since using brute force method in a personal computer becomes unfeasible as typical number of particles involved in molecular simualtions are around 100k. The reason behind this lies in the method self_distance_array in MDAnalysis. It stores the distances for all the possible combination, which consumes most of the memory.

![alt text](/images/230518_paircon-bm-memall.PNG "Memory Consumption for contact pair search")

Cellgrid currently uses similar cutoff distance and grid size for the list creation. It can be understood that this methodology will result in huge execution time for the case when the cutoff distance is very small as well as for the case when the grid size is very high. The latter case results in generation of large list  while the former results in large number of computations between two cells. As the particle density increases, larger grid size results in increase in execution time as can be seen from the figure.   
![alt text](/images/230518_paircon-bm-cutoffcg.PNG "Execution time for contact pair search using cellgrid with different cutoff radius")

2. Similar method has been utilized to obtain the execution time for selection of particles around a group of particles using the above listed methods. For this purpose, a uniform spherical protein particle is created at the centre of the box with 100 particles distributed inside a sphere, while a uniform distribution of solvent particles is generated outside the sphere. The selection attribute is similar to 'around' keyword in MDAnalysis. The results are for the selection distance around 6 units around the spherical protein particle. As opposed to previous use case, cellgrid takes a huge amount of time whereas KDtree is the best method for selection procedure. 

Execution Time             |  Memory Consumption        | Cellgrid Cutoff Radius
:-------------------------:|:-------------------------: |:--------------------------:
![alt text](/images/230518_sphsel-bm-timeall.PNG "Execution time for particle selection")  |  ![alt text](/images/230518_sphsel-bm-memall.PNG "Memory Consumption for particle selection")  | ![alt text](/images/230518_sphsel-bm-cutoffcg.PNG "Execution time for particle selection using cellgrid with different cutoff radius") 

However, for dense system, cellgrid takes lower time as compared to KDtree. 

3. For a inhomogeneous system i.e. two different groups of particles distributed in slabs are analyzed with respect to memory consumption and execution time. A similar trend as obtained in previous selection case is reproduced in this as well, such that KDtree is the best candidate for number of particles ranging from 100 to 100k in a simulation box. 

           
:-------------------------:|:-------------------------: |:--------------------------:
![alt text](/images/230518_slabsel-bm-timeall.PNG "Execution time for particle selection")  |  ![alt text](/images/230518_slabsel-bm-time-cgkd.PNG "Execution time for particle selection (Comparison between KDtree and Cellgrid")  | ![alt text](/images/230518_slabsel-bm-cutoffcg.PNG "Execution time for particle selection using cellgrid with different cutoff radius") 


All this for now. See you in the next blog post. 
