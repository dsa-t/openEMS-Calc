cmake_minimum_required(VERSION 3.15)
project(emswrap)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find Boost libraries
find_package(Boost REQUIRED COMPONENTS program_options NO_MODULE)
# include_directories(${Boost_INCLUDE_DIRS})


# Find hdf5
find_package(HDF5 COMPONENTS C HL REQUIRED)

INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIRS})
link_directories(${HDF5_LIBRARY_DIRS})
# hdf5 compat
ADD_DEFINITIONS( -DH5_USE_16_API )


# Use HighFive
set(HIGHFIVE_USE_BOOST OFF)
set(HIGHFIVE_UNIT_TESTS OFF)
set(HIGHFIVE_BUILD_DOCS OFF)
set(HIGHFIVE_EXAMPLES OFF)
set(HighFive_FIND_QUIETLY ON)
add_subdirectory(highfive)

INCLUDE_DIRECTORIES(highfive) 


# Find CSXCAD
find_library(CSXCAD_LIBRARIES 
  NAMES CSXCAD
  HINTS ${CSXCAD_ROOT_DIR}/lib${LIB_SUFFIX}
  NO_CMAKE_FIND_ROOT_PATH
)
message(STATUS "CSXCAD_LIBRARIES: ${CSXCAD_LIBRARIES}" )
#TODO test if CSXCADs was found
find_path(CSXCAD_INCLUDE_DIR 
    NAMES ContinuousStructure.h 
    HINTS ${CSXCAD_ROOT_DIR}/include
    PATH_SUFFIXES "CSXCAD" ${CSXCAD_INCLUDE_DIR}
    NO_CMAKE_FIND_ROOT_PATH
)
message(STATUS "CSXCAD_INCLUDE_DIR: ${CSXCAD_INCLUDE_DIR}" )
INCLUDE_DIRECTORIES( ${CSXCAD_INCLUDE_DIR} ) 


# Find openEMS
find_library(openEMS_LIBRARIES
    NAMES openEMS
)
message(STATUS "openEMS_LIBRARIES: ${openEMS_LIBRARIES}" )
find_path(openEMS_INCLUDE_DIR 
    NAMES openems.h 
    HINTS ${openEMS_ROOT_DIR}/include
    PATH_SUFFIXES "openEMS" ${openEMS_INCLUDE_DIR}
)
message(STATUS "openEMS_INCLUDE_DIR: ${openEMS_INCLUDE_DIR}" )
INCLUDE_DIRECTORIES( ${openEMS_INCLUDE_DIR} )

# Find nf2ff
find_library(nf2ff_LIBRARIES
    NAMES nf2ff
)
message(STATUS "nf2ff_LIBRARIES: ${nf2ff_LIBRARIES}" )


# Define the source files for the project (all .cpp files in the src directory)
file(GLOB SOURCES "src/*.cpp")

# Add the include directory for header files
include_directories(include)

# Create the executable target
add_library(emswrap STATIC ${SOURCES})

# Link against the required libraries
target_link_libraries(emswrap
    ${CSXCAD_LIBRARIES}
    ${openEMS_LIBRARIES}
    ${nf2ff_LIBRARIES}
    HighFive
)

target_include_directories(emswrap PUBLIC 
    include
    ${CSXCAD_INCLUDE_DIR}
    ${openEMS_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
)

if(MSVC)
    # Disallow implicit linking for Boost
    target_compile_definitions( emswrap PUBLIC BOOST_ALL_NO_LIB )
endif()

install(TARGETS emswrap)