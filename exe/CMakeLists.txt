# exe-level CMakeLists


# For every *.cpp-file (every file corresponding to one of the executables) ...
foreach(EXECUTABLE_SOURCE ${PROJECT_EXECUTABLE_SOURCE_FILES})
    # ... extract the filename without extension (NAME_WE) as a name for our executable ...
    get_filename_component(EXECUTABLE_NAME ${EXECUTABLE_SOURCE} NAME_WE)

    # ... add an executable based on the test source ...
    add_executable(${EXECUTABLE_NAME} ${EXECUTABLE_SOURCE})

    # ... configure (include headers and link libraries) the test ...
    configure_executable(${EXECUTABLE_NAME})
endforeach()
