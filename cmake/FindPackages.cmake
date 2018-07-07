# In this CMake file, we will find all required packages


# Find the Boost package and some required components
find_package(Boost REQUIRED COMPONENTS system thread program_options)

# Find Eigen3
find_package(Eigen3 3.3.4 REQUIRED)

# Find hf
find_package(hf 3.1.0 REQUIRED)

# Find the libint integral wrapper (also includes support for integral transformations)
find_package(libwint 3.1.0 REQUIRED)

# Find bmqc for bitset manipulations
find_package(bmqc 1.2.2 REQUIRED)

# Find numopt
find_package(numopt 1.4.0 REQUIRED)

# Find cpputil
find_package(cpputil 1.4.0 REQUIRED)
