cmake_minimum_required(VERSION 2.8.3)
project(qptest)

set (CMAKE_CXX_STANDARD 11)

include_directories(
	/usr/local/include)

add_executable(qptest
  minsnap.cpp
)

target_link_libraries(qptest
  qpOASES
)
