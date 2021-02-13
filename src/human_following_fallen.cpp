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

#define BEHAVIOR_SIT 0
#define BEHAVIOR_STAND 1
#define BEHAVIOR_WALK 2

#define MODE_STAND 0
#define MODE_START 1

// Define properties for dilations
const int points_per_circle = 5;
bg::strategy::buffer::join_miter join_strategy_input;
bg::strategy::buffer::end_flat end_strategy_input;
bg::strategy::buffer::point_circle point_strategy_input;
bg::strategy::buffer::side_straight side_strategy_input;

class HumanFollowingFallenNode {
	public:
		// Constructor
		HumanFollowingFallenNode(ros::NodeHandle* nodehandle) : nh_(*nodehandle) {
			// Find parameters
			nh_.getParam("pub_twist_topic", pub_twist_topic_);
			nh_.getParam("pub_behaviorID_topic", pub_behaviorID_topic_);
			nh_.getParam("pub_behaviorMode_topic", pub_behaviorMode_topic_);
			nh_.getParam("sub_laser_topic", sub_laser_topic_);
			nh_.getParam("sub_robot_topic", sub_robot_topic_);
			nh_.getParam("sub_semantic_topic", sub_semantic_topic_);
			nh_.getParam("sub_human_topic", sub_human_topic_);
			nh_.getParam("world_frame_id", world_frame_id_);
			nh_.getParam("odom_frame_id", odom_frame_id_);
			nh_.getParam("camera_optical_frame_id", camera_optical_frame_id_);
			nh_.getParam("laser_frame_id", laser_frame_id_);

			nh_.getParam("RobotRadius", RobotRadius_);
			nh_.getParam("ObstacleDilation", ObstacleDilation_);
			nh_.getParam("WalkHeight", WalkHeight_);

			nh_.getParam("AllowableRange", AllowableRange_);
			nh_.getParam("CutoffRange", CutoffRange_);

			nh_.getParam("ForwardLinCmdLimit", ForwardLinCmdLimit_);
			nh_.getParam("BackwardLinCmdLimit", BackwardLinCmdLimit_);
			nh_.getParam("AngCmdLimit", AngCmdLimit_);
			nh_.getParam("RFunctionExponent", RFunctionExponent_);
			nh_.getParam("Epsilon", Epsilon_);
			nh_.getParam("VarEpsilon", VarEpsilon_);
			nh_.getParam("Mu1", Mu1_);
			nh_.getParam("Mu2", Mu2_);
			nh_.getParam("SemanticMapUpdateRate", DiffeoTreeUpdateRate_);
			DiffeoParams_ = DiffeoParamsClass(RFunctionExponent_, Epsilon_, VarEpsilon_, Mu1_, Mu2_, {{-100.0, -100.0}, {100.0, -100.0}, {100.0, 100.0}, {-100.0, 100.0}, {-100.0, -100.0}});

			nh_.getParam("LinearGain", LinearGain_);
			nh_.getParam("AngularGain", AngularGain_);

			nh_.getParam("Goal_x", Goal_x_);
			nh_.getParam("Goal_y", Goal_y_);
			Goal_.set<0>(Goal_x_);
			Goal_.set<1>(Goal_y_);
			nh_.getParam("Tolerance", Tolerance_);

			nh_.getParam("LowpassCutOff", LowpassCutoff_);
			nh_.getParam("LowpassSampling", LowpassSampling_);
			nh_.getParam("LowpassOrder", LowpassOrder_);
			nh_.getParam("LowpassSamples", LowpassSamples_);

			nh_.getParam("DebugFlag", DebugFlag_);

			// Initialize publishers
			pub_behaviorID_ = nh_.advertise<std_msgs::UInt32>(pub_behaviorID_topic_, 1, true);
			pub_behaviorMode_ = nh_.advertise<std_msgs::UInt32>(pub_behaviorMode_topic_, 1, true);
			pub_twist_ = nh_.advertise<geometry_msgs::Twist>(pub_twist_topic_, 1, true);

			// Register callbacks
			message_filters::Subscriber<sensor_msgs::LaserScan> sub_laser(nh_, sub_laser_topic_, 1);
			message_filters::Subscriber<nav_msgs::Odometry> sub_robot(nh_, sub_robot_topic_, 1);
			typedef message_filters::sync_policies::ApproximateTime<sensor_msgs::LaserScan,nav_msgs::Odometry> SyncPolicy;
			typedef message_filters::Synchronizer<SyncPolicy> Sync;
			boost::shared_ptr<Sync> sync;
			sync.reset(new Sync(SyncPolicy(100000), sub_laser, sub_robot));
			sync->registerCallback(boost::bind(&HumanFollowingFallenNode::control_callback, this, _1, _2));

			ros::Subscriber sub_semantic = nh_.subscribe(sub_semantic_topic_, 1, &HumanFollowingFallenNode::diffeo_tree_update, this);

			ros::Subscriber sub_human = nh_.subscribe(sub_human_topic_, 1, &HumanFollowingFallenNode::human_update, this);

			// Publish zero commands
			publish_behavior_id(BEHAVIOR_STAND);
			ros::Duration(5.0).sleep();
			publish_behavior_id(BEHAVIOR_WALK);
			publish_twist(0.0, 0.0);

			// Spin
			ros::MultiThreadedSpinner spinner(4);
			spinner.spin();
		}

		void publish_twist(double LinearCmd, double AngularCmd) {
			geometry_msgs::Twist commandTwist;
			commandTwist.linear.x = LinearCmd;
			commandTwist.angular.z = AngularCmd;
			pub_twist_.publish(commandTwist);
			return;
		}

		void publish_behavior_id(uint16_t BehaviorIdCmd) {
			std_msgs::UInt32 commandId;
			commandId.data = BehaviorIdCmd;
			pub_behaviorID_.publish(commandId);
			return;
		}

		void publish_behavior_mode(uint16_t BehaviorModeCmd) {
			std_msgs::UInt32 commandMode;
			commandMode.data = BehaviorModeCmd;
			pub_behaviorMode_.publish(commandMode);
			return;
		}

		void diffeoTrees_cout(std::vector<std::vector<PolygonClass>> diffeoTreeArray) {
			// Mapper
			std::ofstream svg("/home/kodlab/tree.svg");
			bg::svg_mapper<point> mapper(svg, 1000, 1000);
			std::vector<polygon> polygon_vector, polygon_tilde_vector;
			std::vector<point> point_vector;

			// Print polygon information
			std::cout << "Number of polygons: " << diffeoTreeArray.size() << std::endl;
			for (size_t i = 0; i < diffeoTreeArray.size(); i++) {
				std::cout << "Now printing tree for polygon " << i << std::endl;
				for (size_t j = 0; j < diffeoTreeArray[i].size(); j++) {
					std::cout << "Polygon " << j << " index: " << diffeoTreeArray[i][j].get_index() << std::endl;
					std::cout << "Polygon " << j << " depth: " << diffeoTreeArray[i][j].get_depth() << std::endl;
					std::cout << "Polygon " << j << " predecessor: " << diffeoTreeArray[i][j].get_predecessor() << std::endl;
					std::cout << "Polygon " << j << " radius: " << diffeoTreeArray[i][j].get_radius() << std::endl;
					std::cout << "Polygon " << j << " center: " << bg::dsv(diffeoTreeArray[i][j].get_center()) << std::endl;
					std::vector<point> polygon_vertices = diffeoTreeArray[i][j].get_vertices();
					polygon_vertices.push_back(polygon_vertices[0]);
					std::vector<point> polygon_vertices_tilde = diffeoTreeArray[i][j].get_vertices_tilde();
					polygon_vertices_tilde.push_back(polygon_vertices_tilde[0]);
					std::vector<point> augmented_polygon_vertices = diffeoTreeArray[i][j].get_augmented_vertices();
					augmented_polygon_vertices.push_back(augmented_polygon_vertices[0]);
					polygon polygon_vertices_polygon = BoostPointToBoostPoly(polygon_vertices);
					polygon polygon_vertices_tilde_polygon = BoostPointToBoostPoly(polygon_vertices_tilde);
					polygon augmented_polygon_vertices_polygon = BoostPointToBoostPoly(augmented_polygon_vertices);
					polygon_vector.push_back(polygon_vertices_polygon);
					polygon_tilde_vector.push_back(polygon_vertices_tilde_polygon);
					point_vector.push_back(diffeoTreeArray[i][j].get_center());
					std::cout << "Polygon " << j << " vertices: " << bg::dsv(polygon_vertices_polygon) << std::endl;
					std::cout << "Polygon " << j << " augmented vertices: " << bg::dsv(augmented_polygon_vertices_polygon) << std::endl;
					std::cout << "Polygon " << j << " size of r_t " << diffeoTreeArray[i][j].get_r_t().size() << std::endl;
					std::cout << "Polygon " << j << " size of r_n " << diffeoTreeArray[i][j].get_r_n().size() << std::endl;
					std::cout << "Polygon " << j << " collar: " << bg::dsv(polygon_vertices_tilde_polygon) << std::endl;
					std::cout << "Polygon " << j << " size of r_tilde_t " << diffeoTreeArray[i][j].get_r_tilde_t().size() << std::endl;
					std::cout << "Polygon " << j << " size of r_tilde_n " << diffeoTreeArray[i][j].get_r_tilde_n().size() << std::endl;
					std::cout << "Polygon " << j << " polygons are valid: " << bg::is_valid(polygon_vertices_polygon) << " " << bg::is_valid(polygon_vertices_tilde_polygon) << std::endl;
					std::cout << " " << std::endl;
				}
			}

			// Plot
			std::cout << polygon_vector.size() << std::endl;
			for (size_t i = 0; i < polygon_vector.size(); i++) {
				mapper.add(polygon_vector[i]);
				mapper.add(polygon_tilde_vector[i]);
				mapper.add(point_vector[i]);
			}
			for (size_t i = 0; i < polygon_vector.size(); i++) {
				mapper.map(polygon_vector[i], "fill-opacity:0.3;fill:rgb(51,51,153);stroke:rgb(51,51,153);stroke-width:5", 5);
				mapper.map(polygon_tilde_vector[i], "fill-opacity:0.3;fill:rgb(255,0,0);stroke:rgb(255,0,0);stroke-width:3", 3);
				mapper.map(point_vector[i], "fill-opacity:0.3;fill:rgb(153,204,0);stroke:rgb(153,204,0);stroke-width:5", 5);
			}
		}

		void human_update(const object_pose_interface_msgs::KeypointDetections3D::ConstPtr& human_data) {
			/**
			 * Function that updates the target destination based on the position of the human
			 * 
			 * Input:
			 * 	1) human_data: A KeypointDetections3D object
			 */

			if (human_data->detections.size() == 0) {
				return;
			} else {
				// Get previous start signal
				bool previous_start_signal;
				{
					std::lock_guard<std::mutex> lock(mutex_);
					previous_start_signal = StartSignal_;
				}

				// Find human position in the frame of the camera by taking the mean of all keypoints
				std::vector<double> mean_x(human_data->detections.size(), 0.0);
				std::vector<double> mean_y(human_data->detections.size(), 0.0);
				std::vector<double> mean_z(human_data->detections.size(), 0.0);
				std::vector<double> distance_from_camera(human_data->detections.size(), 0.0);
				for (size_t i = 0; i < human_data->detections.size(); i++) {
					for (size_t j = 0; j < human_data->detections[i].x.size(); j++) {
						mean_x[i] = mean_x[i] + human_data->detections[i].x[j];
						mean_y[i] = mean_y[i] + human_data->detections[i].y[j];
						mean_z[i] = mean_z[i] + human_data->detections[i].z[j];
					}
					mean_x[i] = mean_x[i]/human_data->detections[i].x.size();
					mean_y[i] = mean_y[i]/human_data->detections[i].y.size();
					mean_z[i] = mean_z[i]/human_data->detections[i].z.size();

					distance_from_camera[i] = sqrt(pow(mean_x[i],2) + pow(mean_y[i],2) + pow(mean_z[i],2));
				}

				// Find the human that is closest to the camera
				size_t min_dist_element = std::distance(distance_from_camera.begin(), std::min_element(distance_from_camera.begin(), distance_from_camera.end()));

				// Check if human is fallen to initiate start signal
				std::vector<double> global_x(human_data->detections[min_dist_element].x.size(), 0.0);
				std::vector<double> global_y(human_data->detections[min_dist_element].y.size(), 0.0);
				std::vector<double> global_z(human_data->detections[min_dist_element].z.size(), 0.0);
				listener_.waitForTransform(world_frame_id_, camera_optical_frame_id_, human_data->header.stamp, ros::Duration(2.0));
				for (size_t j = 0; j < human_data->detections[min_dist_element].x.size(); j++) {
					geometry_msgs::PointStamped pointCamera, pointMap;
					pointCamera.header.stamp = human_data->header.stamp;
					pointCamera.header.frame_id = camera_optical_frame_id_;
					pointCamera.point.x = human_data->detections[min_dist_element].x[j];
					pointCamera.point.y = human_data->detections[min_dist_element].y[j];
					pointCamera.point.z = human_data->detections[min_dist_element].z[j];
					try {
						listener_.waitForTransform(world_frame_id_, ros::Time(0), camera_optical_frame_id_, human_data->header.stamp, world_frame_id_, ros::Duration(1.0));
						listener_.transformPoint(world_frame_id_, ros::Time(0), pointCamera, world_frame_id_, pointMap);
					} catch (tf::TransformException &ex) {
						ROS_ERROR("%s",ex.what());
						return;
					}
					global_x[j] = pointMap.point.x;
					global_y[j] = pointMap.point.y;
					global_z[j] = pointMap.point.z;
				}
				size_t max_z_element = std::distance(global_z.begin(), std::max_element(global_z.begin(), global_z.end()));
				bool start_signal = false;
				if ((global_z[max_z_element] < 0.8) || (previous_start_signal)) {
					start_signal = true;
				}

				// Update
				{
					std::lock_guard<std::mutex> lock(mutex_);
					Goal_camera_x_ = mean_x[min_dist_element];
					Goal_camera_y_ = mean_y[min_dist_element];
					Goal_camera_z_ = mean_z[min_dist_element];
					Goal_camera_stamp_ = human_data->header.stamp;
					StartSignal_ = start_signal;
					HumanSeen_ = true;
				}
				
				return;
			}
		}

		void diffeo_tree_update(const object_pose_interface_msgs::SemanticMapObjectArray::ConstPtr& semantic_map_data) {
			/**
			 * Function that updates the semantic map polygons to be used by the control callback
			 * 
			 * Input:
			 * 	1) semantic_map_data: A SemanticMapObjectArray object
			 */

			// Check if update is needed
			// std::cout << DiffeoTreeUpdateRate_ << std::endl;
			// std::cout << ros::Time::now().toSec() - DiffeoTreeUpdateTime_ << std::endl;
			if (ros::Time::now().toSec() - DiffeoTreeUpdateTime_ < (1.0/DiffeoTreeUpdateRate_)) {
				return;
			} else {
				// Count time
				double start_time = ros::Time::now().toSec();

				// Initialize polygon lists
				std::vector<polygon> polygon_list;
				std::vector<polygon> polygon_list_merged;

				// Span the incoming message to add all the polygons
				for (size_t i = 0; i < semantic_map_data->objects.size(); i++) {
					// Extract points of the polygon
					std::vector<point> polygon_in_points;
					for (size_t j = 0; j < semantic_map_data->objects[i].polygon2d.polygon.points.size(); j++) {
						polygon_in_points.push_back(point(semantic_map_data->objects[i].polygon2d.polygon.points[j].x, semantic_map_data->objects[i].polygon2d.polygon.points[j].y));
					}
					polygon polygon_in = BoostPointToBoostPoly(polygon_in_points);

					// Dilate the polygon by the robot radius and append it to the polygon list
					multi_polygon output;
					bg::strategy::buffer::distance_symmetric<double> distance_strategy(ObstacleDilation_);
					bg::buffer(polygon_in, output, distance_strategy, side_strategy_input, join_strategy_input, end_strategy_input, point_strategy_input);
					polygon_list.push_back(output.front());
				}

				// Span all the found polygons to check for intersections between the known polygons and keep only the merged polygons
				multi_polygon output_union;
				if (polygon_list.size() >= 1) {
					output_union.push_back(polygon_list.back());
					polygon_list.pop_back();
					while (!polygon_list.empty()) {
						polygon next_polygon = polygon_list.back();
						polygon_list.pop_back();
						multi_polygon temp_result;
						bg::union_(output_union, next_polygon, temp_result);
						output_union = temp_result;
					}
				}
				for (size_t i = 0; i < output_union.size(); i++) {
					// polygon ch_component;
					// bg::convex_hull(output_union[i], ch_component);
					polygon simplified_component;
					bg::simplify(output_union[i], simplified_component, 0.2);
					polygon_list_merged.push_back(simplified_component);
				}

				// Find diffeomorphism trees for all merged polygons
				std::vector<std::vector<PolygonClass>> localDiffeoTreeArray;
				for (size_t i = 0; i < polygon_list_merged.size(); i++) {
					std::cout << bg::dsv(polygon_list_merged[i]) << std::endl;
					std::vector<PolygonClass> tree;
					diffeoTreeConvex(BoostPointToStd(BoostPolyToBoostPoint(polygon_list_merged[i])), DiffeoParams_, &tree);
					localDiffeoTreeArray.push_back(tree);
				}

				// Update
				{
					std::lock_guard<std::mutex> lock(mutex_);
					DiffeoTreeArray_.clear();
					PolygonList_.clear();
					DiffeoTreeArray_.assign(localDiffeoTreeArray.begin(), localDiffeoTreeArray.end());
					PolygonList_.assign(polygon_list_merged.begin(), polygon_list_merged.end());
				}
				
				if (DebugFlag_) {
					diffeoTrees_cout(DiffeoTreeArray_);
				}
				ROS_WARN_STREAM("[Navigation] Updated trees in " << ros::Time::now().toSec()-start_time << " seconds.");

				// Update time
				DiffeoTreeUpdateTime_ = ros::Time::now().toSec();
			}
			return;
		}

		void control_callback(const sensor_msgs::LaserScan::ConstPtr& lidar_data, const nav_msgs::Odometry::ConstPtr& robot_data) {
			/**
			 * Callback function that implements the main part of the reactive controller
			 * 
			 * Input:
			 * 	1) lidar_data: Data received from the LIDAR sensor
			 * 	2) robot_data: Data received from the robot odometry topic
			 */

			// Make a local copy
			std::vector<polygon> localPolygonList;
			std::vector<std::vector<PolygonClass>> localDiffeoTreeArray;
			double localGoalCameraX, localGoalCameraY, localGoalCameraZ;
			ros::Time localGoalCameraStamp;
			bool localStartSignal, human_seen;
			{
				std::lock_guard<std::mutex> lock(mutex_);
				localPolygonList.assign(PolygonList_.begin(), PolygonList_.end());
				localDiffeoTreeArray.assign(DiffeoTreeArray_.begin(), DiffeoTreeArray_.end());
				localGoalCameraX = Goal_camera_x_;
				localGoalCameraY = Goal_camera_y_;
				localGoalCameraZ = Goal_camera_z_;
				localGoalCameraStamp = Goal_camera_stamp_;
				localStartSignal = StartSignal_;
				human_seen = HumanSeen_;
			}

			// Compute before time
			double before_time = ros::Time::now().toSec();

			// Assuming the incoming odometry message is in the odom frame, transform to map frame
			geometry_msgs::PoseStamped odomPose, mapPose;
			odomPose.header.stamp = ros::Time(0);
			odomPose.header.frame_id = odom_frame_id_;
			odomPose.pose = robot_data->pose.pose;
			try {
				listener_.waitForTransform(world_frame_id_, odom_frame_id_, ros::Time(0), ros::Duration(1.0));
				listener_.transformPose(world_frame_id_, odomPose, mapPose);
			} catch (tf::TransformException &ex) {
				ROS_ERROR("%s",ex.what());
				return;
			}

			// Get robot position and orientation
			tf::Quaternion rotation = tf::Quaternion(mapPose.pose.orientation.x, 
													 mapPose.pose.orientation.y,
													 mapPose.pose.orientation.z,
													 mapPose.pose.orientation.w);
			tf::Matrix3x3 m(rotation);
			double roll, pitch, yaw;
			m.getRPY(roll, pitch, yaw);
			double x_robot_position = mapPose.pose.position.x;
			double y_robot_position = mapPose.pose.position.y;

			// Register robot state - Compensate for LIDAR to camera transform
			double RobotPositionX = x_robot_position;
			double RobotPositionY = y_robot_position;
			RobotPitch_ = pitch;
			RobotPosition_.set<0>(RobotPositionX-RobotRadius_*cos(yaw));
			RobotPosition_.set<1>(RobotPositionY-RobotRadius_*sin(yaw));
			RobotOrientation_ = yaw;

			// Transform camera coordinates to global coordinates
			geometry_msgs::PointStamped localGoalCamera, localGoalMap;
			localGoalCamera.header.stamp = localGoalCameraStamp;
			localGoalCamera.header.frame_id = camera_optical_frame_id_;
			localGoalCamera.point.x = localGoalCameraX;
			localGoalCamera.point.y = localGoalCameraY;
			localGoalCamera.point.z = localGoalCameraZ;
			try {
				listener_.waitForTransform(world_frame_id_, ros::Time(0), camera_optical_frame_id_, localGoalCameraStamp, world_frame_id_, ros::Duration(1.0));
				listener_.transformPoint(world_frame_id_, ros::Time(0), localGoalCamera, world_frame_id_, localGoalMap);
				Goal_x_ = localGoalMap.point.x - 0.1*cos(atan2(localGoalMap.point.y-RobotPositionY,localGoalMap.point.x-RobotPositionX));
				Goal_y_ = localGoalMap.point.y - 0.1*sin(atan2(localGoalMap.point.y-RobotPositionY,localGoalMap.point.x-RobotPositionX));
			} catch (tf::TransformException &ex) {
				ROS_ERROR("%s",ex.what());
			}
			if (localStartSignal && human_seen) {
				Goal_.set<0>(Goal_x_);
				Goal_.set<1>(Goal_y_);
			} else {
				Goal_.set<0>(0.0);
				Goal_.set<1>(0.0);
			}

			// Construct LIDAR object
			LIDARClass LIDAR;
			constructLIDAR2D(lidar_data, CutoffRange_, AllowableRange_, RobotPitch_, &LIDAR);

			// Complete LIDAR readings
			completeLIDAR2D(&LIDAR);
			// ROS_INFO_STREAM("[Navigation] Constructed LIDAR with " << LIDAR.RangeMeasurements.size() << " rays and " << LIDAR.Angle.size() << " angles.");

			// Set the LIDAR rays that hit known obstacles to the LIDAR range
			for (size_t i = 0; i < localPolygonList.size(); i++) {
				compensateObstacleLIDAR2D(RobotPosition_, RobotOrientation_, localPolygonList[i], &LIDAR);
			}
			// ROS_INFO_STREAM("[Navigation] Compensated for known obstacles.");

			// Find list of polygon objects in the model layer based on the known obstacles
			std::vector<polygon> KnownObstaclesModel;
			for (size_t i = 0; i < localDiffeoTreeArray.size(); i++) {
				std::vector<double> theta = linspace(-M_PI, M_PI, 15);
				std::vector<std::vector<double>> model_polygon_coords;
				point root_center = localDiffeoTreeArray[i].back().get_center();
				double root_radius = localDiffeoTreeArray[i].back().get_radius();
				for (size_t j = 0; j < theta.size(); j++) {
					model_polygon_coords.push_back({root_center.get<0>()+root_radius*cos(theta[j]), root_center.get<1>()+root_radius*sin(theta[j])});
				}
				KnownObstaclesModel.push_back(BoostPointToBoostPoly(StdToBoostPoint(model_polygon_coords)));
			}
			// ROS_INFO_STREAM("[Navigation] Constructed model space.");

			// Find the diffeomorphism and its jacobian at the robot position, along with the necessary second derivatives
			std::vector<double> RobotPositionTransformed = {RobotPosition_.get<0>(), RobotPosition_.get<1>()};
			std::vector<std::vector<double>> RobotPositionTransformedD = {{1.0, 0.0}, {0.0, 1.0}};
			std::vector<double> RobotPositionTransformedDD = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			for (size_t i = 0; i < localDiffeoTreeArray.size(); i++) {
				OutputStructVector TempTransformation = polygonDiffeoConvex(RobotPositionTransformed, localDiffeoTreeArray[i], DiffeoParams_);

				std::vector<double> TempPositionTransformed = TempTransformation.Value;
				std::vector<std::vector<double>> TempPositionTransformedD = TempTransformation.Jacobian;
				std::vector<double> TempPositionTransformedDD = TempTransformation.JacobianD;

				double res1 = TempPositionTransformedD[0][0]*RobotPositionTransformedDD[0] + TempPositionTransformedD[0][1]*RobotPositionTransformedDD[4] + RobotPositionTransformedD[0][0]*(TempPositionTransformedDD[0]*RobotPositionTransformedD[0][0] + TempPositionTransformedDD[1]*RobotPositionTransformedD[1][0]) + RobotPositionTransformedD[1][0]*(TempPositionTransformedDD[2]*RobotPositionTransformedD[0][0] + TempPositionTransformedDD[3]*RobotPositionTransformedD[1][0]);
				double res2 = TempPositionTransformedD[0][0]*RobotPositionTransformedDD[1] + TempPositionTransformedD[0][1]*RobotPositionTransformedDD[5] + RobotPositionTransformedD[0][0]*(TempPositionTransformedDD[0]*RobotPositionTransformedD[0][1] + TempPositionTransformedDD[1]*RobotPositionTransformedD[1][1]) + RobotPositionTransformedD[1][0]*(TempPositionTransformedDD[2]*RobotPositionTransformedD[0][1] + TempPositionTransformedDD[3]*RobotPositionTransformedD[1][1]);
				double res3 = TempPositionTransformedD[0][0]*RobotPositionTransformedDD[2] + TempPositionTransformedD[0][1]*RobotPositionTransformedDD[6] + RobotPositionTransformedD[0][1]*(TempPositionTransformedDD[0]*RobotPositionTransformedD[0][0] + TempPositionTransformedDD[1]*RobotPositionTransformedD[1][0]) + RobotPositionTransformedD[1][1]*(TempPositionTransformedDD[2]*RobotPositionTransformedD[0][0] + TempPositionTransformedDD[3]*RobotPositionTransformedD[1][0]);
				double res4 = TempPositionTransformedD[0][0]*RobotPositionTransformedDD[3] + TempPositionTransformedD[0][1]*RobotPositionTransformedDD[7] + RobotPositionTransformedD[0][1]*(TempPositionTransformedDD[0]*RobotPositionTransformedD[0][1] + TempPositionTransformedDD[1]*RobotPositionTransformedD[1][1]) + RobotPositionTransformedD[1][1]*(TempPositionTransformedDD[2]*RobotPositionTransformedD[0][1] + TempPositionTransformedDD[3]*RobotPositionTransformedD[1][1]);
				double res5 = TempPositionTransformedD[1][0]*RobotPositionTransformedDD[0] + TempPositionTransformedD[1][1]*RobotPositionTransformedDD[4] + RobotPositionTransformedD[0][0]*(TempPositionTransformedDD[4]*RobotPositionTransformedD[0][0] + TempPositionTransformedDD[5]*RobotPositionTransformedD[1][0]) + RobotPositionTransformedD[1][0]*(TempPositionTransformedDD[6]*RobotPositionTransformedD[0][0] + TempPositionTransformedDD[7]*RobotPositionTransformedD[1][0]);
				double res6 = TempPositionTransformedD[1][0]*RobotPositionTransformedDD[1] + TempPositionTransformedD[1][1]*RobotPositionTransformedDD[5] + RobotPositionTransformedD[0][0]*(TempPositionTransformedDD[4]*RobotPositionTransformedD[0][1] + TempPositionTransformedDD[5]*RobotPositionTransformedD[1][1]) + RobotPositionTransformedD[1][0]*(TempPositionTransformedDD[6]*RobotPositionTransformedD[0][1] + TempPositionTransformedDD[7]*RobotPositionTransformedD[1][1]);
				double res7 = TempPositionTransformedD[1][0]*RobotPositionTransformedDD[2] + TempPositionTransformedD[1][1]*RobotPositionTransformedDD[6] + RobotPositionTransformedD[0][1]*(TempPositionTransformedDD[4]*RobotPositionTransformedD[0][0] + TempPositionTransformedDD[5]*RobotPositionTransformedD[1][0]) + RobotPositionTransformedD[1][1]*(TempPositionTransformedDD[6]*RobotPositionTransformedD[0][0] + TempPositionTransformedDD[7]*RobotPositionTransformedD[1][0]);
				double res8 = TempPositionTransformedD[1][0]*RobotPositionTransformedDD[3] + TempPositionTransformedD[1][1]*RobotPositionTransformedDD[7] + RobotPositionTransformedD[0][1]*(TempPositionTransformedDD[4]*RobotPositionTransformedD[0][1] + TempPositionTransformedDD[5]*RobotPositionTransformedD[1][1]) + RobotPositionTransformedD[1][1]*(TempPositionTransformedDD[6]*RobotPositionTransformedD[0][1] + TempPositionTransformedDD[7]*RobotPositionTransformedD[1][1]);

				RobotPositionTransformedDD[0] = res1;
				RobotPositionTransformedDD[1] = res2;
				RobotPositionTransformedDD[2] = res3;
				RobotPositionTransformedDD[3] = res4;
				RobotPositionTransformedDD[4] = res5;
				RobotPositionTransformedDD[5] = res6;
				RobotPositionTransformedDD[6] = res7;
				RobotPositionTransformedDD[7] = res8;

				RobotPositionTransformedD = MatrixMatrixMultiplication(TempPositionTransformedD, RobotPositionTransformedD);

				RobotPositionTransformed = TempPositionTransformed;
			}

			// Make a point for the transformed robot position
			point RobotPositionTransformedPoint = point(RobotPositionTransformed[0],RobotPositionTransformed[1]);
			// ROS_INFO_STREAM("[Navigation] Found diffeomorphism.");

			// Find alpha1, alpha2, beta1, beta2
			double alpha1 = -(RobotPositionTransformedD[1][0]*cos(RobotOrientation_) + RobotPositionTransformedD[1][1]*sin(RobotOrientation_));
			double beta1 = RobotPositionTransformedDD[0]*pow(cos(RobotOrientation_),2) + (RobotPositionTransformedDD[1]+RobotPositionTransformedDD[2])*sin(RobotOrientation_)*cos(RobotOrientation_) + RobotPositionTransformedDD[3]*pow(sin(RobotOrientation_),2);
			double alpha2 = RobotPositionTransformedD[0][0]*cos(RobotOrientation_) + RobotPositionTransformedD[0][1]*sin(RobotOrientation_);
			double beta2 = RobotPositionTransformedDD[4]*pow(cos(RobotOrientation_),2) + (RobotPositionTransformedDD[5]+RobotPositionTransformedDD[6])*sin(RobotOrientation_)*cos(RobotOrientation_) + RobotPositionTransformedDD[7]*pow(sin(RobotOrientation_),2);

			// Find transformed orientation
			double RobotOrientationTransformed = atan2(RobotPositionTransformedD[1][0]*cos(RobotOrientation_)+RobotPositionTransformedD[1][1]*sin(RobotOrientation_), RobotPositionTransformedD[0][0]*cos(RobotOrientation_)+RobotPositionTransformedD[0][1]*sin(RobotOrientation_));

			// Read LIDAR data in the model space to account for the known obstacles
			LIDARClass LIDARmodel_known;
			LIDARmodel_known.RangeMeasurements = LIDAR.RangeMeasurements;
			LIDARmodel_known.Angle = LIDAR.Angle;
			LIDARmodel_known.Range = LIDAR.Range;
			LIDARmodel_known.Infinity = LIDAR.Infinity;
			LIDARmodel_known.MinAngle = LIDAR.MinAngle;
			LIDARmodel_known.MaxAngle = LIDAR.MaxAngle;
			LIDARmodel_known.Resolution = LIDAR.Resolution;
			LIDARmodel_known.NumSample = LIDAR.NumSample;
			readLIDAR2D(point(RobotPositionTransformed[0], RobotPositionTransformed[1]), RobotOrientationTransformed, KnownObstaclesModel, LIDAR.Range, LIDAR.MinAngle, LIDAR.MaxAngle, LIDAR.NumSample, &LIDARmodel_known);
			// ROS_INFO_STREAM("[Navigation] Constructed known model space LIDAR with " << LIDARmodel_known.RangeMeasurements.size() << " rays and " << LIDARmodel_known.Angle.size() << " angles.");

			// Translate LIDAR data from the unknown obstacles to the transformed robot state
			LIDARClass LIDARmodel_unknown;
			LIDARmodel_unknown.RangeMeasurements = LIDAR.RangeMeasurements;
			LIDARmodel_unknown.Angle = LIDAR.Angle;
			LIDARmodel_unknown.Range = LIDAR.Range;
			LIDARmodel_unknown.Infinity = LIDAR.Infinity;
			LIDARmodel_unknown.MinAngle = LIDAR.MinAngle;
			LIDARmodel_unknown.MaxAngle = LIDAR.MaxAngle;
			LIDARmodel_unknown.Resolution = LIDAR.Resolution;
			LIDARmodel_unknown.NumSample = LIDAR.NumSample;
			translateLIDAR2D(RobotPosition_, RobotOrientation_, point(RobotPositionTransformed[0], RobotPositionTransformed[1]), RobotOrientationTransformed, ObstacleDilation_, &LIDARmodel_unknown);
			// ROS_INFO_STREAM("[Navigation] Constructed unknown model space LIDAR with " << LIDARmodel_unknown.RangeMeasurements.size() << " rays and " << LIDARmodel_unknown.Angle.size() << " angles.");

			// Build final model LIDAR object
			std::vector<double> newRangeMeasurements(LIDAR.RangeMeasurements.size(), 0.0);
			for (size_t i = 0; i < LIDAR.RangeMeasurements.size(); i++) {
				newRangeMeasurements[i] = std::min(LIDARmodel_known.RangeMeasurements[i], LIDARmodel_unknown.RangeMeasurements[i]);
			}
			LIDARClass LIDARmodel(newRangeMeasurements, LIDAR.Range-bg::distance(RobotPositionTransformedPoint, RobotPosition_), LIDAR.Infinity, LIDAR.MinAngle, LIDAR.MaxAngle, LIDAR.Resolution);
			// ROS_INFO_STREAM("[Navigation] Constructed model space LIDAR with " << LIDARmodel.RangeMeasurements.size() << " rays and " << LIDARmodel.Angle.size() << " angles.");

			// Find local freespace; the robot radius can be zero because we have already dilated the obstacles
			polygon LF_model = localfreespaceLIDAR2D(RobotPositionTransformedPoint, RobotOrientationTransformed, 0.0, &LIDARmodel);
			ROS_INFO_STREAM("[Navigation] Constructed model freespace.");

			// Find projected goal
			point LGL_model = localgoal_linearLIDAR2D(RobotPositionTransformedPoint, RobotOrientationTransformed, LF_model, Goal_);
			ROS_INFO_STREAM("[Navigation] Computed linear local goal.");
			point LGA1_model = localgoalLIDAR2D(LF_model, Goal_);
			ROS_INFO_STREAM("[Navigation] Computed angular local goal 1.");
			point LGA2_model = localgoal_angularLIDAR2D(RobotPositionTransformedPoint, RobotOrientationTransformed, LF_model, Goal_);
			ROS_INFO_STREAM("[Navigation] Computer angular local goal 2.");
			point LGA_model(LGA1_model.get<0>(), LGA1_model.get<1>()); // avoid division by zero
			ROS_INFO_STREAM("[Navigation] Computed model space projections.");

			// Compute the basis for the virtual control inputs
			double tV = (LGL_model.get<0>()-RobotPositionTransformed[0])*cos(RobotOrientationTransformed) + (LGL_model.get<1>()-RobotPositionTransformed[1])*sin(RobotOrientationTransformed);
			double tW1 = (LGA_model.get<0>()-RobotPositionTransformed[0])*cos(RobotOrientationTransformed) + (LGA_model.get<1>()-RobotPositionTransformed[1])*sin(RobotOrientationTransformed);
			double tW2 = -(LGA_model.get<0>()-RobotPositionTransformed[0])*sin(RobotOrientationTransformed) + (LGA_model.get<1>()-RobotPositionTransformed[1])*cos(RobotOrientationTransformed);

			// Compute the basis for transforming to actual control inputs
			double e_norm = sqrt(pow((RobotPositionTransformedD[0][0]*cos(RobotOrientation_)+RobotPositionTransformedD[0][1]*sin(RobotOrientation_)),2) + pow((RobotPositionTransformedD[1][0]*cos(RobotOrientation_)+RobotPositionTransformedD[1][1]*sin(RobotOrientation_)),2));
			double dksi_dpsi = MatrixDeterminant(RobotPositionTransformedD)/pow(e_norm,2);
			double DksiCosSin = (alpha1*beta1 + alpha2*beta2)/pow(e_norm,2);
			ROS_WARN_STREAM("[Navigation] Diffeomorphism jacobian determinant is " << MatrixDeterminant(RobotPositionTransformedD) << ", robot position is " << bg::dsv(RobotPosition_) << " and robot positon in model space is " << bg::dsv(RobotPositionTransformedPoint));

			// Compute commands by accounting for limits
			double LinearCtrlGain, AngularCtrlGain;
			std::vector<double> vector_to_check_1 = {LinearGain_, ForwardLinCmdLimit_*e_norm/abs(tV), 0.5*AngCmdLimit_*dksi_dpsi*e_norm/(abs(tV*DksiCosSin))};
			LinearCtrlGain = *std::min_element(vector_to_check_1.begin(), vector_to_check_1.end());

			std::vector<double> vector_to_check_2 = {AngularGain_, 0.5*AngCmdLimit_*dksi_dpsi/(abs(atan2(tW2,tW1)))};
			AngularCtrlGain = *std::min_element(vector_to_check_2.begin(), vector_to_check_2.end());

			// Compute virtual and actual inputs
			double dV_virtual = LinearCtrlGain*tV;
			double LinearCmd = dV_virtual/e_norm;
			double dW_virtual = AngularCtrlGain*atan2(tW2,tW1);
			double AngularCmd = (dW_virtual-LinearCmd*DksiCosSin)/dksi_dpsi;

			// Stop if the distance from the goal is less than delta
			if (bg::distance(RobotPosition_, Goal_) < Tolerance_) {
				LinearCmd = 0.0;
				AngularCmd = 0.0;
				publish_behavior_id(BEHAVIOR_STAND);
				ros::Duration(5.0).sleep();
				publish_behavior_id(BEHAVIOR_SIT);

				ROS_WARN_STREAM("[Navigation] Successfully navigated to goal and stopped...");
			}

			// Apply limits
			if (LinearCmd > ForwardLinCmdLimit_) LinearCmd = ForwardLinCmdLimit_;
			if (LinearCmd < BackwardLinCmdLimit_) LinearCmd = BackwardLinCmdLimit_;
			if (AngularCmd < -AngCmdLimit_) AngularCmd = -AngCmdLimit_;
			if (AngularCmd > AngCmdLimit_) AngularCmd = AngCmdLimit_;

			// Publish twist
			publish_twist(LinearCmd, AngularCmd);

			// Print debug information
			if (DebugFlag_) {
				std::cout << "Local free space in model space: " << bg::dsv(LF_model) << std::endl;
				std::cout << "Local linear goal in model space: " << bg::dsv(LGL_model) << std::endl;
				std::cout << "Local angular goal 1 in model space: " << bg::dsv(LGA1_model) << std::endl;
				std::cout << "Local angular goal 2 in model space: " << bg::dsv(LGA2_model) << std::endl;
				std::cout << "Local angular goal in model space: " << bg::dsv(LGA_model) << std::endl;
				std::cout << "Robot position in model space: " << bg::dsv(RobotPositionTransformedPoint) << std::endl;
			}

			// Print time
			ROS_WARN_STREAM("[Navigation] Linear command is " << LinearCmd << " and angular command is " << AngularCmd);
			ROS_WARN_STREAM("[Navigation] Command update for " << int(localDiffeoTreeArray.size()) << " polygons in " << ros::Time::now().toSec()-before_time << " seconds.");

			return;
		}
	
	private:
		// Nodehandle
		ros::NodeHandle nh_;

		// Parameters
		std::string pub_twist_topic_;
		std::string pub_behaviorID_topic_;
		std::string pub_behaviorMode_topic_;
		std::string sub_laser_topic_;
		std::string sub_robot_topic_;
		std::string sub_human_topic_;
		std::string sub_semantic_topic_;
		std::string world_frame_id_;
		std::string odom_frame_id_;
		std::string camera_optical_frame_id_;
		std::string laser_frame_id_;

		std::string target_object_;
		double target_object_length_;
		double target_object_width_;

		ros::Publisher pub_behaviorID_;
		ros::Publisher pub_behaviorMode_;
		ros::Publisher pub_twist_;

		ros::Publisher pub_semantic_map_;

		double RobotRadius_;
		double ObstacleDilation_;
		double WalkHeight_;
		double AllowableRange_;
		double CutoffRange_;

		double ForwardLinCmdLimit_;
		double BackwardLinCmdLimit_;
		double AngCmdLimit_;

		double RFunctionExponent_;
		double Epsilon_;
		double VarEpsilon_;
		double Mu1_;
		double Mu2_;
		double DiffeoTreeUpdateRate_;
		DiffeoParamsClass DiffeoParams_;

		double LinearGain_;
		double AngularGain_;

		double Goal_x_ = 0.0;
		double Goal_y_ = 0.0;
		point Goal_;
		double Goal_camera_x_ = 0.0;
		double Goal_camera_y_ = 0.0;
		double Goal_camera_z_ = 0.0;
		ros::Time Goal_camera_stamp_;
		double Tolerance_;

		double LowpassCutoff_;
		double LowpassSampling_;
		double LowpassOrder_;
		double LowpassSamples_;

		point RobotPosition_ = point(0.0, 0.0);
		double RobotOrientation_ = 0.0;
		double RobotPitch_ = 0.0;

		std::vector<polygon> PolygonList_;
		std::vector<std::vector<PolygonClass>> DiffeoTreeArray_;

		double DiffeoTreeUpdateTime_ = 0.0;

		bool DebugFlag_ = false;
		bool StartSignal_ = false;
		bool HumanSeen_ = false;

		tf::TransformListener listener_;

		std::mutex mutex_;
};

int main(int argc, char** argv) {
	// ROS setups
	ros::init(argc, argv, "human_following_fallen");

	// ROS nodehandle
	ros::NodeHandle nh("~");

	// Start navigation node
	HumanFollowingFallenNode humanFollowingFallenNode(&nh);

	return 0;
}
