# ReducedFluids

Source code from the SIGGRAPH 2020 technical paper [Constraint Bubbles and Affine Regions: Reduced Fluid Models for Efficient Immersed Bubbles and Flexible Spatial Coarsening](https://cs.uwaterloo.ca/~rgoldade/reducedfluids/)

## Build

To build this project in Houdini (Linux):

1. Install [Eigen](http://eigen.tuxfamily.org/)

2. Make a folder called External the top level of the repository.

3. Clone the [Geometric Multigrd Pressure Solver](https://github.com/rgoldade/GeometricMultigridPressureSolver) into the External folder

4. Install Houdini 20.0 or higher.

5. Go to install folder (/opt/hfs.xx).

6. Type "source houdini_setup" to get the necessary environment variables.

7. Make a build folder the top level of the repository.

8. Run *cmake ..* in the build folder

9. Run *make* in the build folder.

10. Verify that it was added to Houdini by:
  - Launch Houdini.
  - Press "tab" in the Network Editor and select a "DOP Network" and place it in the editor.
  - Jump into the DOP Network, press "tab" again and verify that "HDK Affine Bubble Pressure Solver" is searchable.

To build this project in another OS, please refer to the [Houdini HDK](https://www.sidefx.com/docs/hdk/_h_d_k__intro__compiling.html).
