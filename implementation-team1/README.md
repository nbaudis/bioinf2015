mlalign [![Build Status](https://semaphoreci.com/api/v1/projects/738d699e-a7f2-48b3-9000-ae3132aa3f8c/423554/badge.svg)](https://semaphoreci.com/sgraf/mlalign)
=======

This C/C++ library and console program computes the maximum likelihood alignment 
of two sequences after [tkf91](http://link.springer.com/article/10.1007/BF02193625) 
as was part of our bioinformatics practical assignment.

How To Build
============

This project uses CMake for project file generation. These are the steps to perform a build:

  1. `$ cd /path/to/mlalign`
  2. `$ mkdir build` to have a folder where the project files can reside in an out of tree build
  3. `$ cd build`
  4. `$ cmake .. -DCMAKE_BUILD_TYPE=Release`

Now, there should be appropriate project files in the `build/` directory, depending on which target 
CMake chose for you (or you chose), with which compilation should be straightforward.

The POSIX (Linux, Mac OS X) way:
 
  5. `$ make` possibly with `-jX` to build every X objects in parallel for better multi-core utilization.
  
(e.g. `$ make` or opening it in Visual Studio).

How To Use
==========
For instructions on how to use the standalone *console application*, do `$ ./mlalign` with no arguments. 
The example invocation from the usage instructions: 

```
$ ./mlalign --sequence-1 ACGATA --sequence-2 AGTGGTA --lambda 1 --mu 2 --tau 0.1 --pa 0.25 --pc 0.25 --pg 0.25 --pt 0.25

> Found alignment of length 7 with score 4.341e-12:
> A-CGATA
> AGTGGTA
```

**Be warned however that for actually achieving good performance you should use the C API!** 
Most time will be spent on process creation and IO that could be avoided entirely that way.

For a small code sample of how to wire up with the API in header `mlalign.h` and CMake 
target `libmlalign`, dig into src/main.cpp:main and inspect the bottommost lines for how 
to call the only exported function `mlalign_TKF91`.

Adding `libmlalign` to your dependencies is probably done the easiest via a CMake sub-project.

To run the benchmarks, do
  1. `$ cd tests`
  2. `$ ../build/benchmark`

Ignore the baseline measurements, there is just no other way to tell Celero to benchmark your code.