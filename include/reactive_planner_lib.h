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

// Boost imports
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/polygon.hpp>

// ROS imports
#include <ros/ros.h>
#include <tf/tf.h>
#include <tf/transform_datatypes.h>
#include <tf/transform_listener.h>
#include <tf2/LinearMath/Quaternion.h>
#include <tf2/LinearMath/Matrix3x3.h>
#include <message_filters/subscriber.h>
#include <message_filters/synchronizer.h>
#include <message_filters/sync_policies/approximate_time.h>
#include <sensor_msgs/LaserScan.h>
#include <nav_msgs/Odometry.h>
#include <std_msgs/UInt32.h>
#include <geometry_msgs/Twist.h>
#include <geometry_msgs/PolygonStamped.h>
#include <geometry_msgs/PointStamped.h>
#include <geometry_msgs/Point32.h>
#include <geometry_msgs/PoseStamped.h>
#include <object_pose_interface_msgs/SemanticMapObjectArray.h>
#include <object_pose_interface_msgs/KeypointDetections3D.h>

// Local imports
#include <polygeom_lib.h>

// Other imports
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <mutex>

// Define namespaces
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// Define various geometries
using point = bg::model::point<double, 2, bg::cs::cartesian>;
using polygon = bg::model::polygon<point, false, true>;
using line = bg::model::linestring<point>;
using multi_point = bg::model::multi_point<point>;
using multi_polygon = bg::model::multi_polygon<polygon>;


// RULES: 1) Use RobotPosition as a point
//        2) Use RobotPosition separately from RobotOrientation
//        3) Pass LIDAR as a pointer


std::vector<double> linspace(double a, double b, uint16_t n);


template<typename T>
class SmartVector : public std::vector<T> {
    public:
        // act like operator []
        T operator()(size_t _Pos) {
            return (*this)[_Pos];
        }

        // act like MATLAB operator ()
        SmartVector<T> operator()(std::vector<size_t>& positions) {
            SmartVector<T> sub;
            sub.resize(positions.size());
            size_t sub_i = 0;
            for(std::vector<size_t>::iterator pit = positions.begin(); pit != positions.end() ; pit++, sub_i++) {
                sub[sub_i] = (*this)[*pit];
            }
            return sub;
        }
};


struct OutputStructVector {
    /** 
     * Struct that includes the value, jacobian and derivatives of the jacobian for a vector-valued function at a specific point
     * 
     * Properties:
     *  1) Value: Value of the function
     *  2) Jacobian: Jacobian of the function
     *  3) JacobianD: Derivatives of the jacobian in the order 11_x, 11_y, 12_x, 12_y, 21_x, 21_y, 22_x, 22_y
     */
    std::vector<double> Value;
    std::vector<std::vector<double>> Jacobian;
    std::vector<double> JacobianD;
};


struct OutputStructScalar {
    /** 
     * Struct that includes the value, jacobian and hessian for a scalar-valued function at a specific point
     * 
     * Properties:
     *  1) Value: Value of the function
     *  2) Gradient: Gradient of the function
     *  3) Hessian: Hessian of the function
     */
    double Value;
    std::vector<double> Gradient;
    std::vector<std::vector<double>> Hessian;
};


class LIDARClass {
    /** 
     * Class that describes a LIDAR object and is updated as new measurements are received
     * 
     * Properties:
     *  1) RangeMeasurements: Range measurements received
     *  2) Range: Range of the sensor
     *  3) Infinity: Range to be considered as infinity
     *  4) MinAngle: Minimum angle of the sensor
     *  5) MaxAngle: Maximum angle of the sensor
     *  6) Resolution: Sensor angular resolution
     */
    public:
        std::vector<double> RangeMeasurements; // Array of ranges
        double Range; // double
        double Infinity; // double
        double MinAngle; // double
        double MaxAngle; // double
        double Resolution; // double
        uint16_t NumSample; // integer
        std::vector<double> Angle; // Array of angles

        LIDARClass(std::vector<double> RangeMeasurements_in, double Range_in, double Infinity_in, double MinAngle_in, double MaxAngle_in, double Resolution_in) {
            RangeMeasurements = RangeMeasurements_in;
            Range = Range_in;
            Infinity = Infinity_in;
            MinAngle = MinAngle_in;
            MaxAngle = MaxAngle_in;
            NumSample =  RangeMeasurements_in.size();
            Angle = linspace(MinAngle, MaxAngle, NumSample);
	        Resolution = Angle[1]-Angle[0];
        }

        LIDARClass() {}
};


class DiffeoParamsClass {
    /** 
     * Class that describes the diffeomorphism parameters
     * 
     * Properties:
     *  1) p: R-function exponent
     *  2) epsilon: Distance for the switches
     *  3) varepsilon: Distance allowed to dilate the polygon
     *  4) mu_1: Switch exponent mu1
     *  5) mu_2: Switch exponent mu2
     *  6) workspace: Polygonal boundary of the workspace
     */
    public:
        DiffeoParamsClass() {}
        
        DiffeoParamsClass(double p_in, double epsilon_in, double varepsilon_in, double mu_1_in, double mu_2_in, std::vector<std::vector<double>> workspace_in) {
            this->p = p_in;
            this->epsilon = epsilon_in;
            this->varepsilon = varepsilon_in;
            this->mu_1 = mu_1_in;
            this->mu_2 = mu_2_in;
            this->workspace = workspace_in;
        }

        double get_p() const {return this->p;}
        double get_epsilon() const {return this->epsilon;}
        double get_varepsilon() const {return this->varepsilon;}
        double get_mu_1() const {return this->mu_1;}
        double get_mu_2() const {return this->mu_2;}
        std::vector<std::vector<double>> get_workspace() const {return this->workspace;}

        void set_p(double p_in) {this->p = p_in;}
        void set_epsilon(double epsilon_in) {this->epsilon = epsilon_in;}
        void set_varepsilon(double varepsilon_in) {this->varepsilon = varepsilon_in;}
        void set_mu_1(double mu_1_in) {this->mu_1 = mu_1_in;}
        void set_mu_2(double mu_2_in) {this->mu_2 = mu_2_in;}
        void set_workspace(std::vector<std::vector<double>> workspace_in) {this->workspace = workspace_in;}
    
    private:
        double p;
        double epsilon;
        double varepsilon;
        double mu_1;
        double mu_2;
        std::vector<std::vector<double>> workspace;
};


void completeLIDAR2D(LIDARClass* LIDAR);


void constructLIDAR2D(const sensor_msgs::LaserScan::ConstPtr& DataLIDAR, double CutoffRange, double AllowableRange, double Pitch, LIDARClass* LIDAR);


std::vector<point> obstaclePointsLIDAR2D(point RobotPosition, double RobotOrientation, LIDARClass* LIDAR);


void compensateObstacleLIDAR2D(point RobotPosition, double RobotOrientation, polygon Obstacle, LIDARClass* LIDAR);


void readLIDAR2D(point RobotPosition, double RobotOrientation, std::vector<polygon> Obstacles, double Range, double MinAngle, double MaxAngle, uint16_t NumSample, LIDARClass* virtualLIDAR);


void translateLIDAR2D(point RobotPosition, double RobotOrientation, point RobotPositionTransformed, double RobotOrientationTransformed, double RobotRadius, LIDARClass* LIDAR);


std::vector<bool> localminLIDAR2D(LIDARClass* LIDAR);


polygon localworkspaceLIDAR2D(point RobotPosition, double RobotOrientation, double RobotRadius, LIDARClass* LIDAR);


polygon localfreespaceLIDAR2D(point RobotPosition, double RobotOrientation, double RobotRadius, LIDARClass* LIDAR);


line localfreespace_linearLIDAR2D(point RobotPosition, double RobotOrientation, polygon LF);


line localfreespace_angularLIDAR2D(point RobotPosition, double RobotOrientation, polygon LF, point Goal);


point localgoalLIDAR2D(polygon LF, point Goal);


point localgoal_linearLIDAR2D(point RobotPosition, double RobotOrientation, polygon LF, point Goal);


point localgoal_angularLIDAR2D(point RobotPosition, double RobotOrientation, polygon LF, point Goal);


void diffeoTreeTriangulation(std::vector<std::vector<double>> PolygonVertices, DiffeoParamsClass DiffeoParams, std::vector<TriangleClass> *tree);


void diffeoTreeConvex(std::vector<std::vector<double>> PolygonVertices, DiffeoParamsClass DiffeoParams, std::vector<PolygonClass> *tree);


OutputStructVector polygonDiffeoTriangulation(std::vector<double> Position, std::vector<TriangleClass> DiffeoTree, DiffeoParamsClass DiffeoParams);


OutputStructVector polygonDiffeoConvex(std::vector<double> Position, std::vector<PolygonClass> DiffeoTree, DiffeoParamsClass DiffeoParams);


OutputStructVector triangleDiffeo(std::vector<double> Position, TriangleClass Triangle, DiffeoParamsClass DiffeoParams);


OutputStructVector polygonDiffeo(std::vector<double> Position, PolygonClass PolygonUsed, DiffeoParamsClass DiffeoParams);


OutputStructScalar triangleSwitch(std::vector<double> Position, TriangleClass Triangle, DiffeoParamsClass DiffeoParams);


OutputStructScalar polygonSwitch(std::vector<double> Position, PolygonClass PolygonUsed, DiffeoParamsClass DiffeoParams);


OutputStructScalar triangleDeformingFactor(std::vector<double> Position, TriangleClass Triangle);


OutputStructScalar polygonDeformingFactor(std::vector<double> Position, PolygonClass PolygonUsed);


OutputStructScalar triangleBetaSwitch(std::vector<double> Position, TriangleClass Triangle, DiffeoParamsClass DiffeoParams);


OutputStructScalar polygonBetaSwitch(std::vector<double> Position, PolygonClass PolygonUsed, DiffeoParamsClass DiffeoParams);


OutputStructScalar triangleGammaSwitch(std::vector<double> Position, TriangleClass Triangle, DiffeoParamsClass DiffeoParams);


OutputStructScalar polygonGammaSwitch(std::vector<double> Position, PolygonClass PolygonUsed, DiffeoParamsClass DiffeoParams);


OutputStructScalar triangleOutsideImplicit(std::vector<double> Position, TriangleClass Triangle, DiffeoParamsClass DiffeoParams);


OutputStructScalar polygonOutsideImplicit(std::vector<double> Position, PolygonClass PolygonUsed, DiffeoParamsClass DiffeoParams);


OutputStructScalar triangleInsideImplicit(std::vector<double> Position, TriangleClass Triangle, DiffeoParamsClass DiffeoParams);


OutputStructScalar polygonInsideImplicit(std::vector<double> Position, PolygonClass PolygonUsed, DiffeoParamsClass DiffeoParams);