cmake_minimum_required(VERSION 2.8.3)
project(object_pose_interface_msgs)

find_package(catkin REQUIRED 
  COMPONENTS
  roscpp
  message_generation
  std_msgs
  geometry_msgs
)

add_message_files(
  FILES
  Classification.msg
  KeypointDetection3D.msg
  KeypointDetections3D.msg
  ObjectClass.msg
  SemanticMapObject.msg
  SemanticMapObjectArray.msg
)

generate_messages(
  DEPENDENCIES
  std_msgs
  geometry_msgs
)

catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES object_pose_interface_msgs
  CATKIN_DEPENDS message_runtime std_msgs geometry_msgs
#  DEPENDS system_lib
)
