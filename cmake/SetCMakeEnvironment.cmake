# In this CMake file, all CMake variables will be set


# Parse the project name into uppercase and lowercase.
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPERCASE)  # Uppercase is needed in version.hpp.in
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWERCASE)



# The name of the library should be equal to the project name
if(NOT LIBRARY_NAME)
    set(LIBRARY_NAME ${PROJECT_NAME})
endif()

# We want to make a static library
set(LIBRARY_TYPE STATIC)
set(EXPORT_TYPE ARCHIVE)



# Find the source folder
set(PROJECT_SOURCE_FOLDER ${CMAKE_SOURCE_DIR}/src)

# Find the source files
set(PROJECT_SOURCE_FILES
        ${PROJECT_SOURCE_FOLDER}/BaseCI.cpp
        ${PROJECT_SOURCE_FOLDER}/DOCI.cpp
        ${PROJECT_SOURCE_FOLDER}/FCI.cpp)

# Find the header folder
set(PROJECT_INCLUDE_FOLDER ${CMAKE_SOURCE_DIR}/include)

# Find the header files (not including version.hpp.in)
set(PROJECT_INCLUDE_FILES
        ${PROJECT_INCLUDE_FOLDER}/BaseCI.hpp
        ${PROJECT_INCLUDE_FOLDER}/ci.hpp
        ${PROJECT_INCLUDE_FOLDER}/DOCI.hpp
        ${PROJECT_INCLUDE_FOLDER}/FCI.hpp
        ${PROJECT_INCLUDE_FOLDER}/ONVExpansion.hpp
        ${PROJECT_INCLUDE_FOLDER}/ONVExpansionComponent.hpp
        ${PROJECT_INCLUDE_FOLDER}/version.hpp)

# Find the tests folder
set(PROJECT_TESTS_FOLDER ${CMAKE_SOURCE_DIR}/tests)

# Find the source files for the tests
set(PROJECT_TEST_SOURCE_FILES
        ${PROJECT_TESTS_FOLDER}/DOCI_Davidson_test.cpp
        ${PROJECT_TESTS_FOLDER}/DOCI_dense_test.cpp
        ${PROJECT_TESTS_FOLDER}/DOCI_OO_test.cpp
        ${PROJECT_TESTS_FOLDER}/DOCI_RDM_test.cpp
        ${PROJECT_TESTS_FOLDER}/DOCI_sparse_test.cpp
        ${PROJECT_TESTS_FOLDER}/DOCI_test.cpp
        ${PROJECT_TESTS_FOLDER}/FCI_dense_test.cpp
        ${PROJECT_TESTS_FOLDER}/FCI_Davidson_test.cpp
        ${PROJECT_TESTS_FOLDER}/FCI_RDM_test.cpp
        ${PROJECT_TESTS_FOLDER}/FCI_sparse_test.cpp
        ${PROJECT_TESTS_FOLDER}/FCI_test.cpp
        ${PROJECT_TESTS_FOLDER}/ONVExpansion_dynamic_bitset_test.cpp
        ${PROJECT_TESTS_FOLDER}/ONVExpansion_unsigned_test.cpp
        ${PROJECT_TESTS_FOLDER}/ONVExpansionComponent_dynamic_bitset_test.cpp
        ${PROJECT_TESTS_FOLDER}/ONVExpansionComponent_unsigned_test.cpp)

# Find the executable folder
set(PROJECT_EXECUTABLES_FOLDER ${CMAKE_SOURCE_DIR}/exe)

# Find the source files for the executables
set(PROJECT_EXECUTABLE_SOURCE_FILES
        ${PROJECT_EXECUTABLES_FOLDER}/main.cpp)


# Give the user the option to specify an installation prefix. If not given as -DINSTALLATION_PREFIX, defaults to /usr/local.
if(NOT INSTALLATION_PREFIX)
    set(INSTALLATION_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()
set(PROJECT_INSTALL_DIR ${INSTALLATION_PREFIX}/${PROJECT_NAME_LOWERCASE})
set(INCLUDE_INSTALL_DIR ${PROJECT_INSTALL_DIR}/include)
set(CMAKE_INSTALL_DIR ${PROJECT_INSTALL_DIR}/cmake)
set(LIBRARY_INSTALL_DIR ${PROJECT_INSTALL_DIR}/lib)


# Include the function that configures the executables
include(${CMAKE_SOURCE_DIR}/cmake/ConfigureExecutable.cmake)
