# In this CMake file, we will include the headers and link to the necessary libraries


# Include this project's headers
target_include_directories(${LIBRARY_NAME} PRIVATE ${PROJECT_INCLUDE_FOLDER})

# Include the boost headers (dynamic bitset)
target_include_directories(${LIBRARY_NAME} PUBLIC ${Boost_INCLUDE_DIRS})

# Include Eigen
target_link_libraries(${LIBRARY_NAME} PUBLIC Eigen3::Eigen)

# Include hf
target_include_directories(${LIBRARY_NAME} PUBLIC ${hf_INCLUDE_DIRS})
target_link_libraries(${LIBRARY_NAME} PUBLIC hf)

# Include bmqc
target_include_directories(${LIBRARY_NAME} PUBLIC ${bmqc_INCLUDE_DIRS})
target_link_libraries(${LIBRARY_NAME} PUBLIC bmqc)

# Include Spectra
target_include_directories(${LIBRARY_NAME} PUBLIC ${spectra_INCLUDE_DIRS})

# Include numopt
target_include_directories(${LIBRARY_NAME} PUBLIC ${numopt_INCLUDE_DIRS})
target_link_libraries(${LIBRARY_NAME} PUBLIC numopt)

# Include cpputil
target_include_directories(${LIBRARY_NAME} PRIVATE ${cpputil_INCLUDE_DIRS})
target_link_libraries(${LIBRARY_NAME} PRIVATE cpputil)
