# Reactive Navigation with Semantic Feedback Using ROS

This package can be used for doubly reactive navigation with semantic feedback, using C++ and ROS. 

It has been tested with Ubuntu 18.04 and ROS Melodic, on three different robots: Turtlebot, Ghost Robotics Minitaur&trade; and Ghost Robotics Spirit&trade;.

Maintainer: Vasileios Vasilopoulos <vvasilo@seas.upenn.edu>

![](examples/human_following.gif)

## Relevant publications and packages
The code included in this package has been used in the following papers:
* V. Vasilopoulos, G. Pavlakos, S. L. Bowman, J. D. Caporale, K. Daniilidis, G. J. Pappas, D. E. Koditschek, "Reactive Semantic Planning in Unexplored Semantic Environments Using Deep Perceptual Feedback" (*under review*).
* V. Vasilopoulos, G. Pavlakos, K. Schmeckpeper, K. Daniilidis, D. E. Koditschek, "Reactive Navigation in Partially Familiar Planar Environments Using Semantic Perceptual Feedback", [arXiv](https://arxiv.org/abs/2002.08946): 2002.08946.

The doubly-reactive operations in the model space are based on the papers:
* Arslan, O., and Koditschek, D. E., "Exact Robot Navigation using Power Diagrams", *IEEE International Conference on Robotics and Automation* (ICRA '16), 2016, pp. 1-8.
* Arslan, O., and Koditschek, D. E., "Sensor-based Reactive Navigation in Unknown Convex Sphere Worlds", *The 12th International Workshop on the Algorithmic Foundations of Robotics* (WAFR), 2016.

## Hardware Setup
The package assumes that the robot possesses:
1. a LIDAR sensor, for estimating distance to unknown obstacles.
1. a way of generating a semantic map of its surroundings with familiar obstacles (see details in Semantic SLAM interfaces below).
1. a way of generating its own odometry estimate.

These three inputs are given as topics in the `navigation_*` launch files (see below).

## Prerequisites
* For our experiments, we use the [ZED Mini](https://github.com/stereolabs/zed-ros-wrapper) stereo camera. A resolution of 720HD@60Hz works well with an NVIDIA TX2 or an NVIDIA Xavier (make sure to run `sudo nvpmodel -m 0` to enable maximum performance first). 
* For reading a Hokuyo LIDAR sensor, we use the [urg_node](http://wiki.ros.org/urg_node) package.
* For LIDAR downsampling, we use the (forked and modified) [laser_scan_sparsifier](https://github.com/vvasilo/scan_tools/tree/indigo/laser_scan_sparsifier) package, included in [scan_tools](https://github.com/vvasilo/scan_tools). This package depends on [csm](https://github.com/AndreaCensi/csm) which must be installed first.
* We use the [robot_localization](http://wiki.ros.org/robot_localization) package for fusing odometry inputs from multiple sources.
* We use [Boost Geometry](https://www.boost.org/doc/libs/1_70_0/libs/geometry/doc/html/index.html) for basic operations with planar polygons, which must be already installed in your system. 
* For more advanced computational geometry operations, we use the [CGAL](https://www.cgal.org/index.html) library. See [here](https://www.cgal.org/download.html) for installation instructions.
* We implement the ear clipping triangulation method in C++ using the [earcut.hpp](https://github.com/mapbox/earcut.hpp) package, included [here](include). For the Python implementation, we use the [tripy](https://github.com/linuxlewis/tripy) package.
* Except for the ROS Python packages (already included with ROS), the following Python packages are also needed: `shapely`, `scipy` and `numpy`.
* For properly using the visualization functionalities in [visualization.py](src/libraries/visualization.py), we need the Python modules `matplotlib` and `imagemagick`.
* For benchmark experiments with Vicon, the [motion_capture_system](https://github.com/KumarRobotics/motion_capture_system) is used.

You can install all the prerequisites, by first independently installing the [ZED SDK](https://www.stereolabs.com/developers/) on your machine, and then running the following commands:
```bash
sudo apt-get install ros-melodic-urg-node ros-melodic-robot-localization python-shapely python-scipy python-numpy libcgal-dev
cd ~/catkin_ws/src
git clone https://github.com/stereolabs/zed-ros-wrapper.git
git clone https://github.com/AndreaCensi/csm.git
git clone https://github.com/vvasilo/scan_tools.git
git clone https://github.com/KumarRobotics/motion_capture_system.git
catkin build csm
catkin build
pip install tripy
```

## Installation
Once all the prerequisites above are satisfied, install with
```bash
cd ~/catkin_ws/src
git clone https://github.com/vvasilo/semnav.git
cp -r semnav/extras/object_pose_interface_msgs .
catkin build
```

##  Semantic SLAM interfaces
This package needs an external Semantic SLAM engine, *not included here*. However, *any* such engine can be used. The only restriction is associated with the type of messages used, i.e., the semantic map has to be given in a specific way.

In our implementation, these messages are included in a separate package called `object_pose_interface_msgs`. We include pointers to the necessary message formats in the [extras](extras) folder. 

We provide the semantic map in the form of a `SemanticMapObjectArray` message. Each `SemanticMapObject` in the array has a `classification` and `pose` element, as well as a number of 3D `keypoints`.

Using [semslam_polygon_publisher.py](src/tracking/semslam_polygon_publisher.py), we project those `keypoints` on the horizontal plane of motion, and republish the semantic map object with a CCW-oriented `polygon2d` element (i.e., the projection of this 3D object on the 2D plane). To do so, we use a pre-defined object mesh, given in the form of a .mat file. The `mesh_location` for all objects is defined in the associated `tracking_*` launch file (see below), and we include examples for different objects [here](extras/meshes).

**Note**: If the user knows the 2D polygon directly, the above procedure is not necessary - only the `polygon2d` element of each `SemanticMapObject` is used for navigation.

## Types of files and libraries
* The main reactive planning library is [reactive_planner_lib.cpp](src/libraries/reactive_planner_lib.cpp) (in Python: [reactive_planner_lib.py](src/libraries/reactive_planner_lib.py)), which uses functions from [polygeom_lib.cpp](src/libraries/polygeom_lib.cpp) (in Python: [polygeom_lib.py](src/libraries/polygeom_lib.py)). This file includes the functionality for the diffeomorphism construction using either the ear clipping algorithm (see the functions `diffeoTreeTriangulation` and `polygonDiffeoTriangulation`), or convex decomposition (see the functions `diffeoTreeConvex` and `polygonDiffeoConvex`). In the C++ implementation, [reactive_planner_lib.cpp](src/libraries/reactive_planner_lib.cpp) and [polygeom_lib.cpp](src/libraries/polygeom_lib.cpp) are built together into a single library.
* We include several navigation nodes, depending on each geometric or semantic task:
    1. [navigation.cpp](src/navigation.cpp) is the basic navigation node. The user has to specify a geometric target that the robot needs to reach.
    1. [navigation_humans.cpp](src/navigation_humans.cpp) also uses moving humans as obstacles. It is assumed that the robot detects humans with the standard 24-keypoint interface (ROS message `KeypointDetections3D` included [here](extras/object_pose_interface_msgs)).
    1. [navigation_semantic.cpp](src/navigation_semantic.cpp) lets the robot navigate to a predefined geometric target, until it sees a desired object. Then, it tries to reach an obstacle free side in front or behind that object.
    1. [human_following.cpp](src/human_following.cpp) is the basic human following node. The robot navigates to a predefined geometric target, until it sees a human. Then, it follows the closest human, or navigates to the last position it saw a human.
    1. [human_following_fallen.cpp](src/human_following_fallen.cpp) lets the robot navigate to a predefined geometric target, until it sees a human falling down - then it proceeds to approach him.
    1. [human_following_signal.cpp](src/human_following_signal.cpp) lets the robot navigate to a predefined geometric target, until it sees a human. Then it follows the closest human (or navigates to the last position it saw one), until a stop gesture is given (left or right hand raised) - after that, it navigates back to its starting position.
* The [visualization.py](src/libraries/visualization.py) script visualizes properties of the diffeomorphism construction, using the Python implementations ([reactive_planner_lib.py](src/libraries/reactive_planner_lib.py) and [polygeom_lib.py](src/libraries/polygeom_lib.py)).

**Note**: The Python libraries ([reactive_planner_lib.py](src/libraries/reactive_planner_lib.py) and [polygeom_lib.py](src/libraries/polygeom_lib.py)) are also used in a separate MATLAB simulation (see [semnav_matlab](https://github.com/vvasilo/semnav_matlab)).

## Usage
To use the code on a real robot, you need to launch one of each type of launch files below:
* The files with name `bringup_*` launch the sensors for each corresponding robot. For example, the file [bringup_turtlebot.launch](launch/bringup_turtlebot.launch) launches:
    1. the stereo camera launch file (`zed_no_tf.launch`).
    1. the Vicon launch file (if present).
    1. the Kobuki node to bring up Turtlebot's control.
    1. the [urg_node](http://wiki.ros.org/urg_node) node for the Hokuyo LIDAR sensor.
    1. the [laser_scan_sparsifier](https://github.com/vvasilo/scan_tools/tree/indigo/laser_scan_sparsifier) node for downsampling the LIDAR data.
* The files with name `tracking_*` launch the tracking files needed for semantic navigation. For example, the file [tracking_turtlebot_semslam_onboard.launch](launch/tracking_turtlebot_semslam_onboard.launch) launches: 
    1. the corresponding semantic SLAM launch file from the semantic_slam package.
    1. the necessary tf transforms (e.g., between the camera and the robot and between the LIDAR and the robot) for this particular robot.
    1. the [semslam_polygon_publisher.py](src/tracking/semslam_polygon_publisher.py) node that subscribes to the output of the semantic SLAM and publishes 2D polygons on the plane.
* The files with name `navigation_*` launch the reactive controller. For example, the file [navigation_turtlebot_onboard.launch](launch/navigation_turtlebot_onboard.launch) launches the main navigation [node](src/navigation.cpp) for Turtlebot, which subscribes to:
    1. the local odometry node (in this case provided directly by the ZED stereo camera).
    1. the LIDAR data, after downsampling.
    1. the 2D polygons from [semslam_polygon_publisher.py](src/tracking/semslam_polygon_publisher.py).
    1. necessary tf updates to correct local odometry as new updates from the semantic SLAM pipeline become available.

We also include a [debugging launch file](launch/navigation_debug.launch), that communicates with fake [LIDAR](src/fake_lidar_publisher.cpp), [odometry](src/fake_odometry_publisher.cpp) and [semantic map](src/map_debug.cpp) publishers.

![](examples/spirit_outdoor.gif)