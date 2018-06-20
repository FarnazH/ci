# exe-level CMakeLists.txt

add_executable(doci.x main.cpp)

configure_executable(doci.x)

add_executable(oo_doci_rhf.x oo_doci_rhf_main.cpp)
configure_executable(oo_doci_rhf.x)

add_executable(oo_doci_rhf_10-3_50000.x oo_doci_rhf_10-3_50000_main.cpp)
configure_executable(oo_doci_rhf_10-3_50000.x)

add_executable(oo_doci_fno.x oo_doci_fno_main.cpp)
configure_executable(oo_doci_fno.x)


add_executable(doci_rhf.x doci_rhf_main.cpp)
configure_executable(doci_rhf.x)

add_executable(doci_rhf_10-3_50000.x doci_rhf_10-3_50000_main.cpp)
configure_executable(doci_rhf_10-3_50000.x)

add_executable(doci_fno.x doci_fno_main.cpp)
configure_executable(doci_fno.x)
