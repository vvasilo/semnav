# MIT License (modified)

# Copyright (c) 2020 The Trustees of the University of Pennsylvania
# Authors:
# Vasileios Vasilopoulos <vvasilo@seas.upenn.edu>

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this **file** (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

cmake_minimum_required(VERSION 2.8.3)
project(semnav)

## Compile as C++11, supported in ROS Kinetic and newer
add_compile_options(-g -O3 -ffast-math -Wall)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  geometry_msgs
  nav_msgs
  sensor_msgs
  roscpp
  rosconsole
  rostime
  roscpp_serialization
  std_msgs
  message_filters
  message_generation
  object_pose_interface_msgs
  tf
)

# Find other libraries
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.45.0)
find_package(CGAL REQUIRED COMPONENTS Core)

###################################
## catkin specific configuration ##
###################################
catkin_package(
)

###########
## Build ##
###########
set(MAIN_DIR ${CMAKE_SOURCE_DIR})

## Specify additional locations of header files
include_directories(
  ${MAIN_DIR}/include
  ${catkin_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
)

## Declare C++ libraries
add_library(reactive_planner
  ${MAIN_DIR}/src/libraries/polygeom_lib.cpp
  ${MAIN_DIR}/src/libraries/reactive_planner_lib.cpp
)

target_link_libraries(reactive_planner
  ${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
  CGAL::CGAL
  CGAL::CGAL_Core
)

# Add executables
add_executable(navigation 
  ${MAIN_DIR}/src/navigation.cpp
)

add_executable(navigation_semantic 
  ${MAIN_DIR}/src/navigation_semantic.cpp
)

add_executable(navigation_humans 
  ${MAIN_DIR}/src/navigation_humans.cpp
)

add_executable(human_following 
  ${MAIN_DIR}/src/human_following.cpp
)

add_executable(human_following_signal 
  ${MAIN_DIR}/src/human_following_signal.cpp
)

add_executable(human_following_fallen
  ${MAIN_DIR}/src/human_following_fallen.cpp
)

add_executable(map_debug 
  ${MAIN_DIR}/src/map_debug.cpp
)

add_executable(fake_lidar_publisher 
  ${MAIN_DIR}/src/fake_lidar_publisher.cpp
)

add_executable(fake_odometry_publisher 
  ${MAIN_DIR}/src/fake_odometry_publisher.cpp
)

# Link
target_link_libraries(navigation
  ${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
  reactive_planner
)

target_link_libraries(navigation_semantic
  ${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
  reactive_planner
)

target_link_libraries(navigation_humans
  ${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
  reactive_planner
)

target_link_libraries(human_following
  ${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
  reactive_planner
)

target_link_libraries(human_following_signal
  ${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
  reactive_planner
)

target_link_libraries(human_following_fallen
  ${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
  reactive_planner
)

target_link_libraries(map_debug
  ${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
  reactive_planner
)

target_link_libraries(fake_lidar_publisher
  ${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
  reactive_planner
)

target_link_libraries(fake_odometry_publisher
  ${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
  reactive_planner
)

## Mark executables and/or libraries for installation
install(TARGETS reactive_planner
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
