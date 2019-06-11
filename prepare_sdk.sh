#!/bin/bash
sudo apt-get update
sudo apt-get install gcc
sudo apt-get install g++
sudo apt-get install make
sudo apt-get install cmake
sudo apt-get install libboost-all-dev
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
