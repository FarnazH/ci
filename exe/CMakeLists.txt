# exe-level CMakeLists.txt

add_executable(doci.x main.cpp)

configure_executable(doci.x)

add_executable(oo_doci_rhf.x oo_doci_main_rhf.cpp)
configure_executable(oo_doci_rhf.x)

add_executable(oo_doci_fno.x oo_doci_main_fno.cpp)
configure_executable(oo_doci_fno.x)
