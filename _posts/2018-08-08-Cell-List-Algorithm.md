---
layout: post
title: Cell Lists and advantages
---


We had discussed earlier that cell-lists have algorithmic advantage over naive brute force method as well as tree based data structures for uniform distributed datasets. To reiterate, KDTree require o(NlogN) time for construction of data structure and o(logN) for the search. Cell-list algorithm, on the contrary, require O(N) time for the construction of the data structure, and O(1) time for search for a best case. The best case for cell-list is achievable when the particles are well-dispered in the domain, and the minimum required distance between neighbors is such there exist only one particle in a cell. However, such a case ceases to exist almost always in practical applications. 

At  its core, cell-list data structure includes splitting of the domain into smaller subdomains called cells. This is followed by sorting the atoms into these cells. Once sorted, for any given query, its cell-location can be identified (``cellindex``) and naive brute force can be used to calculate the distances between query and all the particles residing in nearest neighbor cells. As can be contemplated, a competition between cell-list and KDtree should exist since increasing the number of particles per cell would eventually lead to increase in evaluations (``O(N^2)``) as compared to KDTree (``O(logN)``).

From an implementation perspective, the basic implementation of cell-list is quite straightforward and can also be found in Appendix F,  Page 552 of
``Understanding Molecular Dynamics: From Algorithm to Applications`` by Frenkel and Smit. Provided the cutoff radius, cell-size is typically chosen as the cutoff radius. All the particles are then distributed into their respective cells based on their position. Any query for a nearest neighbour / fixed radius search, identifies the cell of the query and searches all the neighboring cells only i.e. 27 cells. This reduces the computational cost tremendously as compared to naive brute force calculations. 

While its implementation is also relatively straightforward for orthogonal boxes, another challenge remains is to make cell-list work with triclinic boxes. The challenge of different systems arise due to the requirement to handle Periodic boundary conditions in Molecular dynamics. To explain in brief, a periodic boundary condition in terms of distance calculations means a minimum distance between two particles and their periodic images should be considered as the distance between the particles. Inorder to satisfy the minimum image convention, the particles must be translated in all three directions equivalent to the box vectors, and corresponding distances should be evaluated to get the minimum distance. Inorder to calculate the distances in triclinic systems, an equivalent brick shaped box is constructed and all the particles are wrapped around the brick shaped box. The dimension of the brick-shaped box are defined as the component of box vectors in respective directions i.e. ``[[ax, ay, az], [bx, by, bz], [cx, cy, cz]]`` as the box vectors, then the boundary of the brick shaped box is defined by ``[-ax/2, ax/2], [-bx/2, by/2], [-cz/2, cz/2]]``. To evaluate the minimum distances, every displacement vector is modified by the addition/substraction of box vectors such that the distance becomes less than half the box width in every direction. In python, it can be written as follows:

.. code-block ::python:
    
    for i in range(natoms):
        for m in range(DIM - 1, -1, -1):
            while bbox_coords[i, m] < 0:
                for d in range(m+1):
                    bbox_coords[i, d] += self.c_pbcbox.box[m][d]
            while bbox_coords[i, m] >= self.c_pbcbox.box[m][m]:
                for d in range(m+1):
                    bbox_coords[i, d] -= self.c_pbcbox.box[m][d]


Another practical problem is that the cell-list data structure can become huge for small cutoff radius. This is prevented by defining a ``max_gridsize`` which limits the number of cells within a domain. A class ``FastNS`` is defined which contains all the serach functions, while other helper Classes like ``PBCBox`` handles the brick-shaped box and evaluation of minimum distance for periodic/non-periodic calculations, ``NSResults`` is a class which initializes and maintains containers to access the results of search query, and ``NSGrid`` contains all the grid related functions to distribute atoms in cells, define cell index and more. For reference, the actual file ``nsgrid`` can be found [here](https://github.com/MDAnalysis/mdanalysis/blob/develop/package/MDAnalysis/lib/nsgrid.pyx). Few of the relevant snippets are shown below:

In ``NSGrid`` object, a linear mapping of ever coordinate to cellindex in defined as :

.. code-block ::python:

    cdef int coord2cellid(self, vec coord) nogil:
        return <int> (coord[ZZ] / cellsize[ZZ]) * (cell_offsets[ZZ]) +\
               <int> (coord[YY] / cellsize[YY]) * cell_offsets[YY] + \
               <int> (coord[XX] / cellsize[XX])

However, this is valid only under the assumption that all the ``coord`` are within the brick shaped box spanned by cells. Inorder to distribute the coordinates to their respective cell

.. code-block ::python:
    
    for i in range(ncoords):

        # Add bead to grid cell
        cellindex = self.cellids[i]
        self.beadids[cellindex * self.nbeads_per_cell + beadcounts[cellindex]] = i
        beadcounts[cellindex] += 1

This snippet, finds the ``cellid`` of every coordinate, and traverses the ``beadids`` array which is a list of all the atoms present in a particluar cell. A fixed size of ``nbeads_per_cell`` is used to make a uniform array, where ``nbeads_per_cell`` is the maximum number of atoms in any cell. Once the particles are distributed in the grid, search can be performed using search functions defined in ``FastNS`` class. 

Primarily, two methods are introduced in the ``FastNS`` class with different objectives. The ``search`` method searches the already populated grid with a list of queries to return all the pairs between list of query coordinates and populated coordinates within a fixed radius. ``self_search`` searches the coordinates within the grid to find all the neighbours.  While in first look, ``self_search`` might look as a subset of ``search`` method, but time can be halved directly for this special case as distance between every pair needs to be calculated just once. This small trick saves a significant amount of time from the ``search`` method. 

The ``search`` method finds the pairs and distances by creating another grid with the ``NSGrid`` object followed by searching querines from one grid against initially populated grid. In contrast, the ``self_search`` doesn't involve creating another grid, but searches within itself. This is the layout of the Grid-Search method which is supposed to increase the performance of distance evaluations. We shall compare the performance results in the next blog post. 

Stay tuned to see the performance improvements!!
