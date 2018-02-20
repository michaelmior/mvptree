# MVPTree C library 

* Version: 1.0.0
* Date: 2010/12/17
* Creator: D. Grant Starkweather 
* License: GPLv3
* Contact: starkd88@gmail.com

[![Build Status](https://travis-ci.org/michaelmior/mvptree.svg?branch=master)](https://travis-ci.org/michaelmior/mvptree)<Paste>

## Background

The MVP tree is a distance-based data structure for the storage and retrieval of n-dimensional data points.
It relies on the relative distances from selected vantage points to index the points into a tree-like hierarchy.
It thus cuts the search space into distinct 'hyper-spheres' around each vantage point.  

`libmvptree.a` is a generic implementation of the MVP tree.
It allows the user to define a distance function, the type of data and array length (e.g. its bit width for each data element - 1,2,4 and 8), as well as experiment with various tree shapes (e.g. branch factor, leaf capacity, and a path length variable to save the distances between each point and all all the vantage points).

-------------------------------------------------------------------------------

## Platforms

This release should work fine on all Linux/Unix platforms.
Successful compilation and testing
has been achieved on Windows using Cygwin.
However, msys/mingw is still a problem due to the memory mapping functions in Windows.

## Installation

1. Type `make all` to build the `libmvptree.a` library and test programs.
   Run `./testmvp` to do a basic test of the library.
   More involved tests can be done with `./testmvp2` to test it with various number of randomly simulated data points.
   Run it without arguments to see what options are available.

   Note: For the testing, a specified number of uniformly random data points are generated and added to the tree.
   Then a cluster of data points around another randomly chosen point is generated and added to a tree; each element in these data points is a Poisson distributed random variable to serve as a difference from the central cluster point's respective element.
   The point that serves as the center of the cluster is then used to retrieve k nearest neighbors (in this case, the number in the cluster) from the tree.
   For the test to be successful, all data points must be retrieved.

2. Type `make imget` to build the `imget` image indexing program.

3. `make install` to install in the target directory.
   You might want to edit the `Makefile` to change the `DESTDIR` variable from `/usr/local/lib`.

4. Run `./testmvp` to run the test program.

## API

A demo of API use exists in the `testmvp.c` file.  

## References

Bozkaya, Tolga; Ozsoyoglu, Meral 1999."Indexing Large Metric Spaces for Similarity
Search Queries". ACM Transactions in Database Systems, Vol. 24, No. 3, September 1999,
pg. 361-404.
