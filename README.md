# Tetrahedral_trees

We address the problem of performing efficient spatial and topological queries on large tetrahedral meshes with arbitrary topology and complex boundaries. Such meshes arise in several application domains, such as 3D Geographic Information Systems (GISs), scientific visualization, and finite element analysis. To this aim, we propose Tetrahedral trees, a family of spatial indexes based on a nested space subdivision (an octree or a kD-tree) and defined by several different subdivision criteria. We provide efficient algorithms for spatial and topological queries on Tetrahedral trees and compare to state-of-the-art approaches. Our results indicate that Tetrahedral trees are an improvement over R*-trees for querying tetrahedral meshes; they are more compact, faster in many queries, and stable at variations of construction thresholds. They also support spatial queries on more general domains than topological data structures, which explicitly encode adjacency information for efficient navigation but have difficulties with domains with a non-trivial geometric or topological shape.

### Acknowledgments ###

This work has been partially supported by the US National Science Foundation under grant number IIS-1910766 and by the University of Maryland under the 2017-2018 BSOS Dean Research Initiative Program. It has also been performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344.

### Publications ###

**Main paper**

- **Tetrahedral Trees: A Family of Hierarchical Spatial Indexes for Tetrahedral Meshes**  
R. Fellegara, L. De Floriani, P. Magillo, and K. Weiss
*ACM Transaction on Spatial Algorithms and Systems 6, 4, Article 23, 34 pages, 2020* - [doi](https://doi.org/10.1145/3385851)

**Paper describing the Stellar decomposition idea and the compressed encoding used in Tetrahedral trees library**

- **The Stellar tree: a Compact Representation for Simplicial Complexes and Beyond**  
R. Fellegara, K. Weiss, and L. De Floriani  
*Arxiv e-prints, 2019 (v2)* - [doi](https://arxiv.org/abs/1707.02211)

### Features ###

+ Three spatial indexes based on
    * point threshold (PR-T tree)
    * triangle threshold (PMR-T tree)
    * point and triangle thresholds (PM-T tree)
+ Two spatial decomposition
    * quadtree
    * kD-tree
+ Spatial queries
    * point location
    * box query
    * line query
+ Topological queries
    * batched/ranged co-boundary queries
    * batched/ranged adjacency queries
    * batched/ranged link extraction for vertices

### How to compile ###

The library requires the [boost library](http://www.boost.org/) (for dynamic_bitset class), [qmake](https://doc.qt.io/archives/qt-4.8/qmake-manual.html), and [BitMagic library](http://bitmagic.io/) installed in your system.

In Debian-based OSes (including Ubuntu and Linux Mint) these requirements can be satisfied by running the following command in the terminal:
```
#!

sudo apt install libboost-all-dev qt5-qmake bmagic
```

Once installed, execute from the command line in the root folder of the project the following command:
```
#!

qmake
```
and once configured execute:
```
#!

make
```
This latter command generates an executable in `dist` folder.

The compilation has been test on linux systems.

### Use the main library ###

In the bin folder there is the main executable file named `tetrahedral_trees` that contains the whole library. 
For a complete list of the command line options just type in the command line the following command:
```
./tetrahedral_trees
```

### Supported File Formats ###

The library supports files in `.ts` format, that are simple ASCII files containing the explicit representation of vertices and tetrahedral
```
nV mT              -  number of vertices (nV) and tetrahedra (mT)
x1 y1 z1 f1        -  x y z coordinates for each vertex and field value
x2 y2 z2 f2
.  .  .  .
xn yn zn fn

v11 v21 v31 v41    -  vertices (v1 v2 v3 v4) of the first tetrahedron
v12 v22 v32 v42
.   .   .   .
v1m v2m v3m v4m

```
