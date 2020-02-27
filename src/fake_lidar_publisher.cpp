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

class FakeLidarPublisherNode {
	public:
		// Constructor
		FakeLidarPublisherNode(ros::NodeHandle* nodehandle) : nh_(*nodehandle) {
			// Find parameters
			nh_.getParam("pub_lidar_topic", pub_lidar_topic_);

			// Initialize publishers
			pub_lidar_ = nh_.advertise<sensor_msgs::LaserScan>(pub_lidar_topic_, 50);

			// LIDAR data
			unsigned int num_readings = 100;
			double laser_frequency = 30;
			double ranges[num_readings];
			double intensities[num_readings];
			ros::Rate r(laser_frequency);

			// Spin
			while (ros::ok()) {
				for (unsigned int i = 0; i < num_readings; ++i) {
					ranges[i] = 20.0;
					intensities[i] = 100.0;
				}
				sensor_msgs::LaserScan scan;
				scan.header.stamp = ros::Time::now();
				scan.header.frame_id = "laser";
				scan.angle_min = -2.35;
				scan.angle_max = 2.35;
				scan.angle_increment = 4.7/(num_readings-1);
				scan.time_increment = (1/laser_frequency) / (num_readings);
				scan.range_min = 20.0;
				scan.range_max = 20.0;
				scan.ranges.resize(num_readings);
				scan.intensities.resize(num_readings);
				for (unsigned int i = 0; i < num_readings; ++i) {
					scan.ranges[i] = ranges[i];
					scan.intensities[i] = intensities[i];
				}
				publish_lidar(scan);
				r.sleep();
			}
		}

		void publish_lidar(sensor_msgs::LaserScan lidar) {
			pub_lidar_.publish(lidar);
			return;
		}
	
	private:
		// Nodehandle
		ros::NodeHandle nh_;

		// Parameters
		std::string pub_lidar_topic_;
		ros::Publisher pub_lidar_;
};

int main(int argc, char** argv) {
	// ROS setups
	ros::init(argc, argv, "fake_lidar_publisher");

	// ROS nodehandle
	ros::NodeHandle nh("~");

	// Start fake LIDAR publisher node
	FakeLidarPublisherNode fakeLidarPublisher(&nh);

	return 0;
}