# In this CMake file, we will find all required packages


# Find the Boost package - needed for unittests
find_package(Boost REQUIRED COMPONENTS program_options)

# Find Eigen3
find_package(Eigen3 3.3.4 REQUIRED)

# Find hf
find_package(hf 3.0.0 REQUIRED)

# Find the libint integral wrapper (also includes support for integral transformations)
find_package(libwint 3.0.0 REQUIRED)

# Find bmqc for bitset manipulations
find_package(bmqc 1.0.1 REQUIRED)

# Find numopt
find_package(numopt 1.2.0 REQUIRED)

# Find cpputil
find_package(cpputil 1.3.0 REQUIRED)
