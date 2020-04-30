# ReducedFluids

Source code from the SIGGRAPH 2020 technical paper [Constraint Bubbles and Affine Regions: Reduced Fluid Models for Efficient Immersed Bubbles and Flexible Spatial Coarsening](https://cs.uwaterloo.ca/~rgoldade/reducedfluids/)

## Build

To build this project in Houdini (Linux):

1. Install [Eigen](http://eigen.tuxfamily.org/)

2. Clone the [Geometric Multigrd Pressure Solver](https://github.com/rgoldade/GeometricMultigridPressureSolver) into the External folder

3. Install Houdini 18.0 or higher.

4. Go to install folder (/opt/hfs.xx).

5. Type "source houdini_setup" to get the necessary environment variables.

6. Make a build folder the top level of the repository.

7. Run make in the build folder.

8. Verify that it was added to Houdini by:
  - Launch Houdini.
  - Press "tab" in the Network Editor and select a "DOP Network" and place it in the editor.
  - Jump into the DOP Network, press "tab" again and verify that "HDK Affine Bubble Pressure Solver" is searchable.
