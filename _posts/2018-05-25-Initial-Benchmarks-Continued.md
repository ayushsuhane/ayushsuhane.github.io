---
layout: post
title: Initial Benchmarks (Continued ..)
---

Last post concluded few initial benchmarks comparing KDtree, Cellgrid and Brute force method for two use case (1)Contact searches (2) Near neighbour selection. A obvious conclusion that cell list and KDtree are better alternatives as compared to brute force was quantified in terms of execution time and memory consumption. Some more aspects are dealt in this post. 

Shifting the focus to contact searches, since cellgrid is more relevant for that particular use case, it can be observed that cellgrid becomes beneficial with increase in number of particles for a chosen grid size. As bin size and cutoff radius are attached to same value as of now, it is relevant to check the variation of execution time due to different bin sizes/cutoff distance.  

Execution Time             |  Memory Consumption        
:-------------------------:|:-------------------------: 
![alt text](/images/250518_paircon-bm-cutoffcg.PNG "Execution time for particle selection using cellgrid for different binsize")  |  ![alt text](/images/250518_paircon-bm-memcutoff.PNG "Memory Consumption for particle selection using cellgrid for different binsize") 

It can be seen that small values lead to longer execution time, so does for larger cutoff radius. Furthermore, the memory consumption is also high in case for large cellsize, as it requires more number of particles i.e. more space for pair distances. The reason for long execution time for small cutoff radius is due to longer building time of cellgrid which can also be observed by evaluating the build time at different radius. The memory consumption of building is almost constant. As building requires only assigning particles to theier respective cells based on the location. This assignment scales with O(n), which results in a small slope in log scale.  Interesting thing to note is building time of KDtree is independent of the selection radius and is close to the lowest time among all the bin sizes used here for Cellgrid. 

Build Time             |  Memory Consumption (Build)        
:-------------------------:|:-------------------------: 
![alt text](/images/250518_build-bm-cutoffcg.PNG "Build  time for particle selection using cellgrid for different binsize")  |  ![alt text](/images/250518_build-bm-cutoffcg-mem.PNG "Memory Consumption (Build) for particle selection using cellgrid for different binsize") 

The next post will discuss about the optimization of cellgrid, and quantify the optimum bin size to make the build independent of the number of particles.





