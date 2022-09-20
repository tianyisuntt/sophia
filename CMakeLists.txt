project(sophia)
cmake_minimum_required(VERSION 2.8)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -DNDEBUG")

#	
#	Specify paths to include files of PHAT and Gudhi (and decomment) to activate corresponding options:
#	
#set(PHAT /usr/local/include/)
#set(GUDHI /usr/local/include/)

if(PHAT AND EXISTS "${PHAT}/phat/representations/bit_tree_pivot_column.h")
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_REQUIRED_FLAGS ${OpenMP_CXX_FLAGS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    include_directories(${PHAT})
    add_definitions(-DPHAT)
  else()
    message("OpenMP not found, PHAT Option is disabled.")
  endif()
else()
  message("No or wrong path indicated, PHAT Option is disabled.")
endif()

if(GUDHI AND EXISTS "${GUDHI}/gudhi/Persistent_cohomology.h")
#	If the path for the Boost libraries is not standart, specify it manually (and decomment):
# set(BOOST_LIBRARYDIR /usr/local/lib/boost_stage64)
#	Also possible for the include files:
# set(BOOST_INCLUDE_DIRS /usr/local/include)
  find_package(Boost 1.6 REQUIRED COMPONENTS system filesystem chrono timer REQUIRED)
  if(Boost_FOUND)
    set(Boost_CXX_FLAGS "-DBOOST_SYSTEM_NO_DEPRECATED")
    set(CMAKE_REQUIRED_FLAGS ${Boost_CXX_FLAGS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Boost_CXX_FLAGS}")
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
    LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
    include_directories(${GUDHI})
    add_definitions(-DGUDHI)
  else()
    message("Boost Library not found, Gudhi Option is disabled.")
  endif()
else()
  message("No or wrong path indicated, Gudhi Option is disabled.")
endif()

include_directories(include)
aux_source_directory(./src SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})

