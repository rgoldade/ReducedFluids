##
## Affine dense bubble solver
##

set(library_name HDK_AffineBubblePressureSolver)

# Add a library and its source files.
add_library( ${library_name} SHARED
    ../External/GeometricMultigridPressureSolver/Source/HDK_GeometricMultigridOperators.h
    ../External/GeometricMultigridPressureSolver/Source/HDK_GeometricMultigridOperators.cpp
    ../External/GeometricMultigridPressureSolver/Source/HDK_GeometricMultigridPoissonSolver.h
    ../External/GeometricMultigridPressureSolver/Source/HDK_GeometricMultigridPoissonSolver.cpp
    ../External/GeometricMultigridPressureSolver/Source/HDK_Utilities.h
    ../External/GeometricMultigridPressureSolver/Source/HDK_Utilities.cpp
    HDK_AffineBubblePressureSolver.cpp
    HDK_AffineBubblePressureSolver.h
    HDK_ReducedFluidUtilities.cpp
    HDK_ReducedFluidUtilities.h
)

# Link against the Houdini libraries, and add required include directories and
# compile definitions.
target_link_libraries( ${library_name} Houdini )

# Include ${CMAKE_CURRENT_BINARY_DIR} for the generated header.
target_include_directories( ${library_name} PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
)

# Sets several common target properties, such as the library's output directory.
houdini_configure_target( ${library_name} )

##
## Affine free surface solver
##

set(library_name HDK_AffineFreeSurfacePressureSolver)

# Add a library and its source files.
add_library( ${library_name} SHARED
	../External/GeometricMultigridPressureSolver/Source/HDK_Utilities.h
	../External/GeometricMultigridPressureSolver/Source/HDK_Utilities.cpp
    HDK_AffineFreeSurfacePressureSolver.cpp
    HDK_AffineFreeSurfacePressureSolver.h
    HDK_ReducedFluidUtilities.cpp
    HDK_ReducedFluidUtilities.h
)

# Link against the Houdini libraries, and add required include directories and
# compile definitions.
target_link_libraries( ${library_name} Houdini )

# Include ${CMAKE_CURRENT_BINARY_DIR} for the generated header.
target_include_directories( ${library_name} PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
)

# Sets several common target properties, such as the library's output directory.
houdini_configure_target( ${library_name} )

##
## Constraint bubble solver
##

set(library_name HDK_ConstraintBubblePressureSolver)

# Add a library and its source files.
add_library( ${library_name} SHARED
    ../External/GeometricMultigridPressureSolver/Source/HDK_GeometricMultigridOperators.h
    ../External/GeometricMultigridPressureSolver/Source/HDK_GeometricMultigridOperators.cpp
    ../External/GeometricMultigridPressureSolver/Source/HDK_GeometricMultigridPoissonSolver.h
    ../External/GeometricMultigridPressureSolver/Source/HDK_GeometricMultigridPoissonSolver.cpp
    ../External/GeometricMultigridPressureSolver/Source/HDK_Utilities.h
    ../External/GeometricMultigridPressureSolver/Source/HDK_Utilities.cpp
    HDK_ConstraintBubblePressureSolver.cpp
    HDK_ConstraintBubblePressureSolver.h
    HDK_ReducedFluidUtilities.cpp
    HDK_ReducedFluidUtilities.h
)

# Link against the Houdini libraries, and add required include directories and
# compile definitions.
target_link_libraries( ${library_name} Houdini )

# Include ${CMAKE_CURRENT_BINARY_DIR} for the generated header.
target_include_directories( ${library_name} PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
)

# Sets several common target properties, such as the library's output directory.
houdini_configure_target( ${library_name} )

##
## Two phase solver
##

set(library_name HDK_TwoPhasePressureSolver)

# Add a library and its source files.
add_library( ${library_name} SHARED
    ../External/GeometricMultigridPressureSolver/Source/HDK_Utilities.h
    ../External/GeometricMultigridPressureSolver/Source/HDK_Utilities.cpp
    HDK_ReducedFluidUtilities.cpp
    HDK_ReducedFluidUtilities.h
    HDK_TwoPhasePressureSolver.cpp
    HDK_TwoPhasePressureSolver.h
)

# Link against the Houdini libraries, and add required include directories and
# compile definitions.
target_link_libraries( ${library_name} Houdini )

# Include ${CMAKE_CURRENT_BINARY_DIR} for the generated header.
target_include_directories( ${library_name} PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
)

# Sets several common target properties, such as the library's output directory.
houdini_configure_target( ${library_name} )