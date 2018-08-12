---
layout: post
title: Augment Coordinates for Implementation of Periodic Boundary Conditions
excerpt: "Handling Periodic boundary conditions with augmenting particle coordinates"
---

Periodic Boundary conditions are one of the essential requirements of Molecular Dynamics Simulations. Majority of the MD simulations are executed with PBC in all three directions. Therefore, all the distance computations in MDAnalysis should also support PBC in triclinic simulation boxes. In simpler terms, a particle exiting from one of the face of a box will also enter simulataneously from the corresponding opposite face. In three dimensions, this creates 26 exact replicas of any particle in the neighbourhood of the box under consideration. Inorder to calculate the distance between two particles under the influence of periodic boundary conditions, a minimum image convention is followed. Minimum image distance is the smallest possible distance between two particles in the box as well as between the images in the surrounding duplicate boxes. This naive implementation of PBC for distance calculation increases the computation cost by 27 times.

Present implementation to handle PBC using KDTree in MDAnalysis is a more advanced approach. Rather than evaluating all the images, only the relevant images are generated and checked based on the relative position of both the particles. For a given box dimensions, its reciprocal vectors are evaluated. If `[[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]]` are the dimensions of the box such that `[a, b, c]` are the three vectors of the box, the reciprocal vectors `[ra, rb, rc]` are evaluated as:

``` math
ra = b X c
rb = c X a
rc = a X b
```
The reciprocal vectors represent the normal vectors to the plane formed by the its member vectors i.e. `ra` is the plane normal vector of a plane formed by `b` and `c` vectors. Since the problem is inherently three dimensional, both the particles in consideration have a location defined by their respective three dimensional vectors. A dot product with the reciprocal vectors yield the distance of the coordinate from the corresponding face along the direction of normal (reciprocal) vector. This concept can be exploited to limit the generation of images which are within a certain distance from the box faces. The cutoff distance of `[a/2, b/2, c/2]` would be the condition of equivalency of this approach with the naive PBC implementation described above.    

However, typical problems in MDAnalysis doesnot require to go to the extreme of choosing half the box distance as the cutoff distance. A conditional evaluation of image distance only of the coordinate lie inside the cutoff distance saves a significant amount of time. The current algorithm in MDAnalysis is pure python implementation and can be described as follows:

* Evaluate the reciprocal vectors.
* Find the distance of the coordinates from each of the 6 faces.
* If any of the distance is less than the specified distance, generate images and check the distance based on the minimum image convention.
* Return the minimum distance.

This has two drawbacks though :

* Pure python implementation is slower as the evaluation of minimum distance requires multiple loops.
* Other algorithms for near neighbour searching which donot have specific implementation of PBC cannot be handled direclty.

To make the algorithm more flexible, an approach of augment coordinates is adopted. In this approach, all the particles within the cutoff distance are duplicated as per the Periodic Boundary Conditions. Another mapping of these particles with the original particles is also constructed to undo the augmentation. These duplicated particles along with the original particles can now be operated with any algorithm to evaluate the distances without any special consideratio to PBC. All the indices/coordinates can be mapped back to the original set of particles using the mapping described above. A cython implementation is implemented for triclinic boxes for fast looping based on this [notebook](http://www.richardjgowers.com/2018/06/28/make_halos.html). Furthermore, memoryviews are used here based on this [observation](https://jakevdp.github.io/blog/2012/08/08/memoryview-benchmarks/). 

Similar algorithm described above is used here for cython implementation. To evaluate the distance of particle from each of the faces, dot product operation of two vectors i.e. `[x, y, z]` (particle coordinate) and `[x, y, z] - ([a1, a2, a3] + [b1, b2, b3] + [c1, c2, c3])` i.e distance of coordinate from the opposite diagonal box is performed. Furthermore, duplicates in the diagonal images of box are generated based on the condition if the particle distance evaluated from more than one faces in less than the specified distance. For example, if the particle lies close to the xy plane but far from the yz plane, then corresponding duplicate will be generated in xy plane near the opposite plane. The function to augment the coordinates and undo the arguments can be written as

``` python
	%%cython --annotate

cimport cython
import cython

cimport numpy as np
import numpy as np

from libc.math cimport sqrt

@cython.boundscheck(False)
@cython.wraparound(False)
def make_halo(float[:, ::1] coordinates, float[:,::1] box, float[:,::1] reciprocal, float r):
    """Calculate augmented coordinate set
    
    Parameters
    ----------
    coordinates : np.ndarray
      coordinates to augment
    box : np.ndarray
      size of box
    rm : np.ndarray
      reciprocal matrix
    r : float
      thickness of halo region to buffer by
      
    Returns
    -------
    augmented : np.ndarray
      coordinates of the new augmented coordinates
    indices : np.ndarray
      original indices of the augmented coordinates
    """
    cdef bint lo_x, hi_x, lo_y, hi_y, lo_z, hi_z
    cdef int i, j, p, N
    cdef float shiftX[3]
    cdef float shiftY[3]
    cdef float shiftZ[3]
    cdef float coord[3]
    cdef float end[3]
    cdef float other[3]
    
    cdef int dim
    dim = coordinates.shape[1]
    # room for adding triclinic support by using (3,) vectors
    for i in range(dim):
        shiftX[i] = box[0, i]
        shiftY[i] = box[1, i]
        shiftZ[i] = box[2, i]
        end[i] = box[0, i] + box[1, i] + box[2,i]
    N = coordinates.shape[0]
    p = 0  # output counter
    # allocate output arrays
    # could be more conservative with this
    # or use C++ vectors + push etc
    cdef float[:, :] output = np.zeros((N, 3), dtype=np.float32)
    cdef int[:] indices = np.zeros(N, dtype=np.int32)

    for i in range(0,N):
        for j in range(3):
            coord[j] = coordinates[i, j]
            other[j] = end[j] - coordinates[i,j]
        # identify the condition 
        lo_x = _dot(&coord[0], &reciprocal[0,0]) <= r
        hi_x = _dot(&other[0], &reciprocal[0,0]) <= r
        lo_y = _dot(&coord[0], &reciprocal[1,0]) <= r
        hi_y = _dot(&other[0], &reciprocal[1,0]) <= r
        lo_z = _dot(&coord[0], &reciprocal[2,0]) <= r
        hi_z = _dot(&other[0], &reciprocal[2,0]) <= r
        
        if lo_x:
            # if X, face piece
            for j in range(3):
                # add to output
                output[p, j] = coord[j] + shiftX[j]
            # keep record of which index this augmented position was created from
            indices[p] = i
            p += 1
           
            if lo_y:
                # if X&Y, edge piece
                for j in range(3):
                    output[p, j] = coord[j] + shiftX[j] + shiftY[j]
                indices[p] = i
                p += 1
                
                
                if lo_z:
                    # if X&Y&Z, corner piece
                    for j in range(3):
                        output[p, j] = coord[j] + shiftX[j] + shiftY[j] + shiftZ[j]
                    indices[p] = i
                    p += 1
                    
                    
                elif hi_z:
                    for j in range(3):
                        output[p, j] = coord[j] + shiftX[j] + shiftY[j] - shiftZ[j]
                    indices[p] = i
                    p += 1
                    

            elif hi_y:
                for j in range(3):
                    output[p, j] = coord[j] + shiftX[j] - shiftY[j]
                indices[p] = i
                p += 1
                
                
                if lo_z:
                    for j in range(3):
                        output[p, j] = coord[j] + shiftX[j] - shiftY[j] + shiftZ[j]
                    indices[p] = i
                    p += 1
                    
                    
                elif hi_z:
                    for j in range(3):
                        output[p, j] = coord[j] + shiftX[j] - shiftY[j] - shiftZ[j]
                    indices[p] = i
                    p += 1
                    
                
            if lo_z:
                for j in range(3):
                    output[p, j] = coord[j] + shiftX[j] + shiftZ[j]
                indices[p] = i
                p += 1
                
               
            elif hi_z:
                for j in range(3):
                    output[p, j] = coord[j] + shiftX[j] - shiftZ[j]
                indices[p] = i
                p += 1
                
                
        elif hi_x:
            for j in range(3):
                output[p, j] = coord[j] - shiftX[j]
            indices[p] = i
            p += 1
            
            
            if lo_y:
                for j in range(3):
                    output[p, j] = coord[j] - shiftX[j] + shiftY[j]
                indices[p] = i
                p += 1
                
                
                if lo_z:
                    for j in range(3):
                        output[p, j] = coord[j] - shiftX[j] + shiftY[j] + shiftZ[j]
                    indices[p] = i
                    p += 1
                    
                    
                elif hi_z:
                    for j in range(3):
                        output[p, j] = coord[j] - shiftX[j] + shiftY[j] - shiftZ[j]
                    indices[p] = i
                    p += 1
                    
                

            elif hi_y:
                for j in range(3):
                    output[p, j] = coord[j] - shiftX[j] - shiftY[j]
                indices[p] = i
                p += 1
                
                
                if lo_z:
                    for j in range(3):
                        output[p, j] = coord[j] - shiftX[j] - shiftY[j] + shiftZ[j]
                    indices[p] = i
                    p += 1
                    
                elif hi_z:
                    for j in range(3):
                        output[p, j] = coord[j] - shiftX[j] - shiftY[j] - shiftZ[j]
                    indices[p] = i
                    p += 1
                    
            if lo_z:
                for j in range(3):
                    output[p, j] = coord[j] - shiftX[j] + shiftZ[j]
                indices[p] = i
                p += 1
                
            elif hi_z:
                for j in range(3):
                    output[p, j] = coord[j] - shiftX[j] - shiftZ[j]
                indices[p] = i
                p += 1
                
        if lo_y:
            for j in range(3):
                output[p, j] = coord[j] + shiftY[j]
            indices[p] = i
            p += 1
            
            if lo_z:
                for j in range(3):
                    output[p, j] = coord[j] + shiftY[j] + shiftZ[j]
                indices[p] = i
                p += 1
                
            elif hi_z:
                for j in range(3):
                    output[p, j] = coord[j] + shiftY[j] - shiftZ[j]
                indices[p] = i
                p += 1
                
        elif hi_y:
            for j in range(3):
                output[p, j] = coord[j] - shiftY[j]
            indices[p] = i
            p += 1
            

            if lo_z:
                for j in range(3):
                    output[p, j] = coord[j] - shiftY[j] + shiftZ[j]
                indices[p] = i
                p += 1
                
            elif hi_z:
                for j in range(3):
                    output[p, j] = coord[j] - shiftY[j] - shiftZ[j]
                indices[p] = i
                p += 1
                
        if lo_z:
            for j in range(3):
                output[p, j] = coord[j] + shiftZ[j]
            indices[p] = i
            p += 1
            
        elif hi_z:
            for j in range(3):
                output[p, j] = coord[j] - shiftZ[j]
            indices[p] = i
            p += 1
            
    return np.array(output[:p]), np.array(indices[:p])


@cython.boundscheck(False)
@cython.wraparound(False)
cdef float _dot(float  *a, float *b):
    """Return dot product of two sequences in range."""
    cdef ssize_t n
    cdef float sum1
    cdef ssize_t dim
    
    dim=3
    sum1 = 0.0
    for n in range(dim):
        sum1 += a[n] * b[n]
    return sum1


@cython.boundscheck(False)
@cython.wraparound(False)
def undo_augment(int[:] results, int[:] translation, int nreal):
    """Translate augmented indices back to originals
    
    Note: modifies results in place!
    
    Parameters
    ----------
    results : ndarray of ints
      indices of coordinates, including "augmented" indices
    translation : ndarray of ints
      original indices of augmented coordinates
    nreal : int
      number of real coordinates, ie values in results equal or larger than this
      need to be translated to their real counterpart
      
    Returns
    -------
    results : ndarray of ints
    """
    cdef int N
    
    N = results.shape[0]
    
    for i in range(N):
        if results[i] >= nreal:
            results[i] = translation[results[i] - nreal]
            
    return results
```    

Followed by the definition of this function, the API for distance evaluations which returns all the pairs within a specified distance using KDTree can be written as:

``` python
def augment_kdtree(coords, cutoff, box):
    """
    box of type [A,B,C,a1,a2,a3]
    """
    box = triclinic_vectors(box)
    # Define reciprocal matrix
    rm = np.zeros(9, dtype=np.float32).reshape(3, 3)
    rm[0] = np.cross(box[1], box[2])
    rm[1] = np.cross(box[2], box[0])
    rm[2] = np.cross(box[0], box[1])
    for i in range(3):
        rm[i] /= norm(rm[i])  # normalize
    aug, idx = make_halo(coords, box, rm, cutoff)
    aug_coord = np.concatenate([coords, aug])
    kdtree = spatial.cKDTree(aug_coord)
    pairs = np.array(list(kdtree.query_pairs(cutoff)), dtype=np.int32)
    if len(pairs) > 1:
        undo_augment(pairs[:, 0], idx, len(coords))
        undo_augment(pairs[:, 1], idx, len(coords))
        pairs = np.unique(np.sort(pairs), axis=0)
    return pairs
```

Similarly, other distance evaluations such as atom selections can be performed using similar methods. To benchmark the cython implementation of KDtree and already implemented PKDtree in MDAnalysis, a use case to identify all the bonds i.e. distance between two particles less than the sum of radius of each particle is chosen. Two different real case are also chosen with ~12k and ~1000k particles to check the extremeties of performance. While more detailed analysis is provided [here](https://github.com/ayushsuhane/Benchmarks_Distance/blob/master/Notebooks/Augment_Functionality.ipynb) and [here](https://github.com/ayushsuhane/Benchmarks_Distance/blob/master/Notebooks/Augment_Functionality.ipynb), major highlights of the benchmark is provided below.

For the two cases described above the guess_bonds function in `MDAnalysis.topology.guessers.py` took `2.86 sec` for smaller system, while the `PKDTree` in `MDAnalysis.lib.pkdtree` took `2.69 sec`, whereas the function defined here took `920 ms` i.e increase in performance by 3 times. For the larger system, the computation time for  ` guess_bonds` was more than `30 min` whereas the current implementation of PKDTree in MDAnalysis took around `7 mins`. For this case, the method described above took `2 mins` i.e. the performance boost of 3 times than already implemented PKDtree (but not in guess bonds), whereas the performance boost is significantly high compared to brute force approach currently implemented in MDAnalysis for guessing the bonds.

The next step is to include this functionality in MDAnalysis and replacement of Bio.KDtree, which is obsolete in BioPython, with scipy.spatial.cKDTree. 

See you next time!!!

