#!/bin/bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --target satellitePRJ -j $(nproc)