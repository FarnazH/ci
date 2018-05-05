
find_package( Boost COMPONENTS system thread program_options REQUIRED )

add_executable( doci.x
  main.cpp )

configure_executable( doci.x )


