// MIT License (modified)

// Copyright (c) 2020 The Trustees of the University of Pennsylvania
// Authors:
// Vasileios Vasilopoulos <vvasilo@seas.upenn.edu>

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this **file** (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <reactive_planner_lib.h>

class MapDebugNode {
	public:
		// Constructor
		MapDebugNode(ros::NodeHandle* nodehandle) : nh_(*nodehandle) {
			// Find parameters
			nh_.getParam("pub_semantic_topic", pub_semantic_topic_);

			// Initialize publishers
			pub_semantic_map_ = nh_.advertise<object_pose_interface_msgs::SemanticMapObjectArray>(pub_semantic_topic_, 1, true);

			// Spin
			while (ros::ok()) {
				// std::vector<std::vector<double>> polygon1 = {{2.518, 0.5048}, {1.83, 0.2963}, {2.0430, -0.2348}, {2.406, -0.8039}, {2.655, -0.0533}, {2.518, 0.5048}};
				// std::vector<std::vector<double>> polygon1 = {{0.0, 0.0}, {5.0, 0.0}, {5.0, 5.0}, {0.0, 5.0}, {0.0, 4.0}, {4.0, 4.0}, {4.0, 1.0}, {0.0, 1.0}, {0.0, 0.0}};
				// std::vector<std::vector<double>> polygon1 = {{2.8, 1.4}, {6.2, 1.4}, {6.2, 10.0}, {2.8, 10.0}, {2.8, 8.6}, {4.8, 8.6}, {4.8, 7.0}, {2.8, 7.0}, {2.8, 4.4}, {4.8, 4.4}, {4.8, 2.8}, {2.8, 2.8}, {2.8, 1.4}};
				// std::vector<std::vector<double>> polygon1 = {{7.0, 9.5}, {6.0, 7.6}, {4.0, 7.2}, {5.4, 5.6}, {5.1, 3.5}, {7.0, 4.4}, {8.9, 3.5}, {8.6, 5.6}, {10.0, 7.2}, {8.0, 7.6}, {7.0, 9.5}};
				// std::vector<std::vector<double>> polygon1 = {{7.0, 7.0}, {7.0, 6.0}, {8.0, 6.0}, {8.0, 1.0}, {7.0, 1.0}, {7.0, 0.0}, {10.0, 0.0}, {10.0, 1.0}, {9.0, 1.0}, {9.0, 6.0}, {10.0, 6.0}, {10.0, 7.0}, {7.0, 7.0}};
				std::vector<std::vector<double>> polygon1 = {{0.0, 0.0}, {0.5, 0.0}, {0.5, -1.0}, {1.5, -1.0}, {1.5, 1.0}, {-1.0, 1.0}, {-1.0, -3.0}, {3.0, -3.0}, {3.0, 1.0}, {2.0, 1.0}, {2.0, -2.0}, {0.0, -2.0}, {0.0, 0.0}};
				object_pose_interface_msgs::SemanticMapObjectArray polygon_list_msg;
				object_pose_interface_msgs::SemanticMapObject polygon1_msg = populate_polygon_msg(polygon1);
				polygon_list_msg.objects.push_back(polygon1_msg);
				publish_semantic_map(polygon_list_msg);
			}
		}

		void publish_semantic_map(object_pose_interface_msgs::SemanticMapObjectArray semantic_map) {
			pub_semantic_map_.publish(semantic_map);
			return;
		}

		object_pose_interface_msgs::SemanticMapObject populate_polygon_msg(std::vector<std::vector<double>> polygon_in) {
			object_pose_interface_msgs::SemanticMapObject polygon_out;
			for (size_t i = 0; i < polygon_in.size(); i++) {
				geometry_msgs::Point32 point_new;
				point_new.x = polygon_in[i][0];
				point_new.y = polygon_in[i][1];
				polygon_out.polygon2d.polygon.points.push_back(point_new);
			}
			return polygon_out;
		}
	
	private:
		// Nodehandle
		ros::NodeHandle nh_;

		// Parameters
		std::string pub_semantic_topic_;
		ros::Publisher pub_semantic_map_;
};

int main(int argc, char** argv) {
	// ROS setups
	ros::init(argc, argv, "map_debug");

	// ROS nodehandle
	ros::NodeHandle nh("~");

	// Start map_debug node
	MapDebugNode mapDebug(&nh);

	return 0;
}
