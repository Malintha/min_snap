cmake_minimum_required(VERSION 2.8.3)
project(qptest)

set(LDFLAGS= "-L. -lqpOASES.so")

set (CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
# set(CMAKE_CXX_FLAGS "-g -O0 -DSF_VISIBILITY -fvisibility=hidden -fno-strict-aliasing")

include_directories(
	/usr/local/include)

add_executable(qptest
  # multi_quad_traj.cpp
  example.cpp
  # main.cpp
  # solver.h
  # solver.cpp
  # yaml_test.h

)

target_link_libraries(qptest
  qpOASES
  yaml
  # pthread
)
