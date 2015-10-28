  1. `$ cd crlibm-11.0beta4`
  2. `$ ./configure`
  3. `$ sudo make install`
  4. `$ cd ..`
  5. `$ mkdir build` to have a folder where the project files can reside in an out of tree build
  6. `$ cd build`
  7. `$ cmake .. -DCMAKE_BUILD_TYPE=Release`
  8. `$ make` possibly with `-jX` to build every X objects in parallel for better multi-core utilization.
  9. `$ ./TKF91Sequential -h` for instructions on how to use this.

