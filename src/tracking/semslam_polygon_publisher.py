#!/usr/bin/env python

"""
MIT License (modified)

Copyright (c) 2019 The Trustees of the University of Pennsylvania
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

# Python imports
import numpy as np
import copy
import shapely.geometry
import scipy.io as sio
import time
from shapely.geometry import Polygon
from shapely.ops import cascaded_union

# ROS imports
import rospy
import tf
import tf2_ros
import std_msgs.msg
from std_msgs.msg import Float32, UInt8, Float32MultiArray
from geometry_msgs.msg import PoseWithCovarianceStamped, TransformStamped, PointStamped, PolygonStamped, Point32
from nav_msgs.msg import Odometry
from object_pose_interface_msgs.msg import SemanticMapObject, SemanticMapObjectArray
from visualization_msgs.msg import MarkerArray
import message_filters




# Global variables
global pub_semantic, mesh_location


def object_subscriber(data):
	"""
	Function that transforms detected objects to 2D polygons
	"""
	# Global variables
	global pub_semantic, mesh_location

	# Check all incoming objects
	for obj_id in range(len(data.objects)):
		# Find current time
		before_time = time.clock()

		# Load mesh faces
		mat_file = sio.loadmat(mesh_location + data.objects[obj_id].classification.type.name + '.mat')
		mesh_faces = mat_file['faces']
		keypoints_raw = mat_file['keypoints_raw']
		keypoints_raw = keypoints_raw-np.array([keypoints_raw.mean(axis=1)]).transpose()

		# Find rotation and translation
		translation = np.array([[data.objects[obj_id].pose.pose.position.x], [data.objects[obj_id].pose.pose.position.y], [data.objects[obj_id].pose.pose.position.z]])
		rotation = tf.transformations.quaternion_matrix([data.objects[obj_id].pose.pose.orientation.x, data.objects[obj_id].pose.pose.orientation.y, data.objects[obj_id].pose.pose.orientation.z, data.objects[obj_id].pose.pose.orientation.w])

		# Find array of points in 3D
		keypoints = np.empty((0,3))
		keypoints_msg = PolygonStamped()
		keypoints_msg.header = data.header
		for point_id in range(np.size(keypoints_raw,1)):
			keypoint_transformed = np.array(translation+np.matmul(rotation[0:3,0:3],np.array([keypoints_raw[:,point_id]]).transpose())).transpose()
			keypoint_to_append = Point32()
			keypoint_to_append.x = keypoint_transformed[0,0]
			keypoint_to_append.y = keypoint_transformed[0,1]
			keypoint_to_append.z = keypoint_transformed[0,2]
			keypoints_msg.polygon.points.append(keypoint_to_append)
			keypoints = np.append(keypoints,keypoint_transformed, axis=0)
		data.objects[obj_id].keypoints = keypoints_msg
		
		# Span the triangular faces of the mesh to create the 2D projection
		polygons_for_union = []
		for face_id in range(mesh_faces.shape[0]):
			polygons_for_union.append(Polygon(np.array([keypoints[mesh_faces[face_id][0]-1][0:2],keypoints[mesh_faces[face_id][1]-1][0:2],keypoints[mesh_faces[face_id][2]-1][0:2]])))
		
		projection_polygon = cascaded_union(polygons_for_union)

		# Simplify and orient the final polygon CCW
		projection_polygon = projection_polygon.simplify(0.05, preserve_topology = True)
		projection_polygon = shapely.geometry.polygon.orient(projection_polygon, 1.0)

		# Return the final polygon
		polygon_out = PolygonStamped()
		polygon_out.header = data.header
		for point_id in range(len(projection_polygon.exterior.coords.xy[0])):
			point_to_append = Point32()
			point_to_append.x = projection_polygon.exterior.coords.xy[0][point_id]
			point_to_append.y = projection_polygon.exterior.coords.xy[1][point_id]
			point_to_append.z = 0.0
			polygon_out.polygon.points.append(point_to_append)
		data.objects[obj_id].polygon2d = polygon_out
		
	
	pub_semantic.publish(data)

	return
	

def init():
	# Global variables
	global pub_semantic, mesh_location

	# Initialize node
	rospy.init_node('semslam_polygon_publisher', anonymous = True)

	# Find parameters
	mesh_location = rospy.get_param('~mesh_location')
	object_poses_topic = rospy.get_param('~object_poses_topic')
	semantic_map_topic = rospy.get_param('~semantic_map_topic')

	# Initialize ROS publishers
	pub_semantic = rospy.Publisher(semantic_map_topic, SemanticMapObjectArray, queue_size=1)

	# Define object pose subscriber
	rospy.Subscriber(object_poses_topic, SemanticMapObjectArray, object_subscriber)

	# Keep running
	rospy.spin()


if __name__ == '__main__':
	try:
		init()
	except rospy.ROSInterruptException: pass
