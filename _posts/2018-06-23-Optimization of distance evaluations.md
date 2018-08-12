---
layout: post
title: Optimization of distance evaluations
excerpt: "Inclusion of capped function"
---

In the last blog post, we discussed the advantages of gridsearch algorithm for three dimensional datasets, but it can be memory expensive based on the required cutoff radius. As KDTree and brute force are implemented and deeply rooted in MDAnalysis, an optimization of selection of algorithm i.e. choosing the algorithm based on certain empirical rules hidden from the user, would be helpful in fast distance evaluations. For this reason, the primary objective is to maintain a function which can handle different algorithms and switch between them based on the data provided to it.
The skeleton of the function included the definition of the function which takes two sets of coordinates i.e. atom positions accessed by `Universe.atoms.positions`, maximum and minimum cutoff radius where minimum cutoff radius is kept as optional, box which can be accessed from `universe.dimensions` but is kept as optional along with an optional handle provided to the user to override the automatic selection of the method. The implememtation of periodic boundary conditions are automatically considered when a box is supplied. Currently, KDTree and Bruteforce method are supported in this function, but the flexible structure allows adding other methods i.e. grid search method once it becomes available.

The main functionality is handled by `cappedfunction(reference, configuration, max_cutoff, min_cutoff=None, box=None, method=None)` which resides at  `MDAnalysis.lib.distances.py`. The function returns array of all the indices of the pairs as well as their corresponding distances which are within the distance specified by `max_cutoff` and `min_cutoff`. Unless the method is specified by the user, a function (`_determine_method`) which automatically choses the method based on the number of particles in either set of coordinates and/or search radius is called by the `cappeddistance` function. This function (`determine_method`) returns the method object which points towards an efficient algorithm. Each algorithm has a different method defined i.e. KDTree has `_pkdtree_capped` and bruteforce has `_bruteforce_capped` method defined in the function as:

``` python
def _bruteforce_capped(reference, configuration, max_cutoff, min_cutoff=None, box=None):
    
    pairs, distance = [], []

    reference = np.asarray(reference, dtype=np.float32)
    configuration = np.asarray(configuration, dtype=np.float32)

    if reference.shape == (3, ):
        reference = reference[None, :]
    if configuration.shape == (3, ):
        configuration = configuration[None, :]

    _check_array(reference, 'reference')
    _check_array(configuration, 'configuration')

    for i, coords in enumerate(reference):
        dist = distance_array(coords[None, :], configuration, box=box)[0]
        if min_cutoff is not None:
            idx = np.where((dist <= max_cutoff) & (dist > min_cutoff))[0]
        else:
            idx = np.where((dist <= max_cutoff))[0]
        for j in idx:
            pairs.append((i, j))
            distance.append(dist[j])
    return pairs, distance
```

and similarly for KDTree. To remove the circular imports within `pkdtree.py` and `distances.py`, local import to the modules are performed in the `_pkdtree_capped` function. The rules are selected such as if either of the particle set consisted of less than 5000 particles or if the number of particles is greater than 100k and the maximum cutoff distance is greater than 10% of the box size naive brute force is selected, whereas kdtree is selected in all the other cases. A notebook with the benchmarks to evaluate these rules can be found [here](https://github.com/ayushsuhane/Benchmarks_Distance/blob/master/Notebooks/Capped_Distance.ipynb). 

For any new method, it needs to be registered into the `_determine_method` function dictionary i.e 

``` python
methods = {'bruteforce': _bruteforce_capped,
            'pkdtree': _pkdtree_capped}
```

and a corresponding function needs to be defined as well. The corresponding rules can be inserted in the `_determine_method` without disturbing the main algorithm for distance evaluations.

This is the brief discussion of the method which has various advantages i.e. it can be used to replace a lot of analysis functions such as guess_bonds, radial distribution function analysis, particle selections etc in an optimized way and totally hidden from the user (along with the flexibility to give user the control to override), and resulting in a more maintainable and clear code. 

This is all there is to it now. We shall discuss it later as well, after the implementation of gridsearch method where optimized rules will be added for gridsearch method in this function as well.






