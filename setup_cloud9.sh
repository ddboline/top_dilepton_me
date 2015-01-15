#!/bin/bash

sudo apt-get install -y gfortran libroot-bindings-python-dev libroot-graf2d-postscript-dev libroot-core-dev libroot-math-physics-dev
sudo apt-get install -y libgsl0-dev liblhapdf-dev make
mkdir -p obj/
make -f bin/Makefile
