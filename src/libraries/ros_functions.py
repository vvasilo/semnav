"""
MIT License (modified)

Copyright (c) 2020 The Trustees of the University of Pennsylvania
Authors:
Vasileios Vasilopoulos <vvasilo@seas.upenn.edu>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this **file** (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

# General ROS and Python imports
import roslib, rospy, struct, math, tf, numpy, os, sys
import scipy.io as sio
import message_filters

# Import ROS messages
from nav_msgs.msg import Odometry
from nav_msgs.msg import Path as NavPath
from std_msgs.msg import Float32, UInt32
from sensor_msgs.msg import LaserScan
from shapely.geometry import Polygon as ShapelyPolygon
from geometry_msgs.msg import Twist, Point32, PoseStamped, PointStamped
from geometry_msgs.msg import PolygonStamped as GeometryMsgsPolygon


def publish_freespace(Publisher,LF,Header,FrameId):
	"""
	Function that publishes a Polygon message showing the current local freespace
	
	Input:
		1) Publisher: PolygonStamped ROS publisher
		2) LF: Local freespace polygon
		3) Header: ROS header to be appended to Polygon message
		4) FrameId: String defining the frame_id of the Polygon message
	"""
	# Publish freespace polygon
	polygon_msg = GeometryMsgsPolygon()
	polygon_msg.header = Header
	polygon_msg.header.frame_id = FrameId
	if LF.any() and ShapelyPolygon(LF).is_valid:
		numvertex = LF.shape[0]
		for i in range(0, numvertex):
			polygon_msg.polygon.points.append(Point32(x=LF[i][0], y=LF[i][1], z=0.))
	else:
		polygon_msg.polygon.points.append(Point32(x=0., y=0., z=0.))
	Publisher.publish(polygon_msg)
	return


def publish_path(Publisher,Path,Header,FrameId):
	"""
	Function that publishes a Path message showing the current path
	
	Input:
		1) Publisher: Path ROS publisher
		2) Path: 2D path as an array of points
		3) Header: ROS header to be appended to Path message
		4) FrameId: String defining the frame_id of the Path message
	"""
	# Publish path
	path_msg = NavPath()
	path_msg.header = Header
	path_msg.header.frame_id = FrameId
	numpoints = Path.shape[0]
	for i in range(0, numpoints):
		pose = PoseStamped()
		pose.header = Header
		pose.pose.position.x = Path[i][0]
		pose.pose.position.y = Path[i][1]
		pose.pose.position.z = 0.
		path_msg.poses.append(pose)
	Publisher.publish(path_msg)
	return


def publish_point(Publisher,Point,Header,FrameId):
	"""
	Function that publishes a Point message showing a position
	Input:
		1) Publisher: PointStamped ROS publisher
		2) Point: Array with the desired point
		3) Header: ROS header to be appended to Point message
		4) FrameId: String defining the frame_id of the Point message
	"""
	point_msg = PointStamped()
	point_msg.header = Header
	point_msg.header.frame_id = FrameId
	point_msg.point.x = Point[0]
	point_msg.point.y = Point[1]
	point_msg.point.z = 0.
	Publisher.publish(point_msg)
	return


def publish_twist(Publisher,LinearCmd,AngularCmd,HeightCmd):
	"""
	Function that publishes a Twist message
	
	Input:
		1) Publisher: Twist ROS publisher
		2) LinearCmd: Linear command
		3) AngularCmd: Angular command
		4) HeightCmd: Robot height command (only applicable for Minitaur)
	"""
    # Publish twist
	twist = Twist()
	twist.linear.x = LinearCmd
	twist.linear.z = HeightCmd
	twist.angular.z = AngularCmd
	Publisher.publish(twist)
	ros_sleep(0.01)
	return


def publish_behavior_id(Publisher,BehaviorIdCmd):
	"""
	Function that publishes a behavior ID command (only applicable for Minitaur)
	
	Input:
		1) Publisher: UInt32 ROS publisher
		2) BehaviorIdCmd: Behavior ID command (selects behavior)
	"""
	Publisher.publish(BehaviorIdCmd)
	ros_sleep(0.01)
	return


def publish_behavior_mode(Publisher,BehaviorModeCmd):
	"""
	Function that publishes a behavior mode command (only applicable for Minitaur)
	
	Input:
		1) Publisher: UInt32 ROS publisher
		2) BehaviorModeCmd: Behavior mode command (starts or stops behavior)
	"""
	Publisher.publish(BehaviorModeCmd)
	ros_sleep(0.01)
	return


def ros_sleep(SleepTime):
	"""
	Function that sleeps for a specified amount of time

	Input:
		1) SleepTime: Time to sleep in sec
	"""
	rospy.sleep(SleepTime)
	return


def get_current_time():
	"""
	Function that returns current ROS time
	"""
	return rospy.get_time()