# exe-level CMakeLists.txt

add_executable(doci.x main.cpp)
configure_executable(doci.x)


add_executable(fci_rhf.x fci_rhf_main.cpp)
configure_executable(fci_rhf.x)


add_executable(fci_lowdin.x fci_lowdin_main.cpp)
configure_executable(fci_lowdin.x)
