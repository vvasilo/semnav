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

class FakeOdometryPublisherNode {
	public:
		// Constructor
		FakeOdometryPublisherNode(ros::NodeHandle* nodehandle) : nh_(*nodehandle) {
			// Find parameters
			nh_.getParam("pub_odom_topic", pub_odom_topic_);

			// Initialize publishers
			pub_odom_ = nh_.advertise<nav_msgs::Odometry>(pub_odom_topic_, 50);

			// Odometry data
			double odom_frequency = 30;
			ros::Rate r(odom_frequency);

			// Spin
			while (ros::ok()) {
				nav_msgs::Odometry odom_msg;
				odom_msg.header.stamp = ros::Time::now();
				odom_msg.child_frame_id = "robot";
				odom_msg.pose.pose.position.x = -1.5;
				odom_msg.pose.pose.position.y = 0.0;
				odom_msg.pose.pose.orientation.w = 1.0;
				publish_odom(odom_msg);
				r.sleep();
			}
		}

		void publish_odom(nav_msgs::Odometry odom_data) {
			pub_odom_.publish(odom_data);
			return;
		}
	
	private:
		// Nodehandle
		ros::NodeHandle nh_;

		// Parameters
		std::string pub_odom_topic_;
		ros::Publisher pub_odom_;
};

int main(int argc, char** argv) {
	// ROS setups
	ros::init(argc, argv, "fake_odometry_publisher");

	// ROS nodehandle
	ros::NodeHandle nh("~");

	// Start fake odometry publisher node
	FakeOdometryPublisherNode fakeOdometryPublisher(&nh);

	return 0;
}