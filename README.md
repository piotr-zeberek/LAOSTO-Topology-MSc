# Source code for "Topological properties of a two-dimensional electron gas at the LaAlO3 /SrTiO3 interface" MSc thesis

The code contains the implementation of the tight-binding model and calculations of topological invariants for 2DEG with specified classes for ToyModel (Ek + Zeeman + Rashba SOC + s-wave SC) and LAO/STO 2DEG. 

## Requirements
- C++17 compiler
- Eigen library (for linear algebra) with unsupported modules
- pybind11 library (for interfacing with Python)
- Python with numpy, scipy, and pfapack libraries (for sparse eigenvalue calculations and pfaffian calculations)
- CMake (for building the project)
