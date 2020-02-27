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


// Define properties for dilations
const int points_per_circle = 5;
bg::strategy::buffer::join_round join_strategy(points_per_circle);
bg::strategy::buffer::end_flat end_strategy;
bg::strategy::buffer::point_circle circle_strategy;
bg::strategy::buffer::side_straight side_strategy;


std::vector<double> linspace(double a, double b, uint16_t n) {
    std::vector<double> array;
    double step = (b-a) / (n-1);
    array.push_back(a);

    while(array.size() < n) {
        a += step;           // could recode to better handle rounding errors
	    array.push_back(a);
    }
    array.pop_back();
    array.push_back(b);
    return array;
}


void completeLIDAR2D(LIDARClass* LIDAR) {
    /**
     * Function that completes missing LIDAR data due to limited field of view
     * 
     * Input: 
     *  1) LIDAR: Incomplete LIDAR object
     */
    if ((LIDAR->MaxAngle-LIDAR->MinAngle) < 2*M_PI) {
        std::vector<double> tempR = LIDAR->RangeMeasurements;
        std::vector<double> tempAngle = LIDAR->Angle;
        double tempResolution = LIDAR->Resolution;

        // Updated LIDAR model
        LIDAR->MaxAngle = LIDAR->MinAngle + 2*M_PI;
        LIDAR->NumSample = uint16_t(round((2*M_PI/tempResolution)+1));
        LIDAR->Resolution = (LIDAR->MaxAngle-LIDAR->MinAngle)/(LIDAR->NumSample-1);
        LIDAR->Angle = linspace(LIDAR->MinAngle, LIDAR->MaxAngle, LIDAR->NumSample);

        // Completed Range Data
        std::vector<double> R(LIDAR->NumSample, LIDAR->Range);
        for (size_t i = 0 ; i < tempAngle.size() ; i++) {
            size_t index = size_t(floor((fmod(tempAngle[i]-LIDAR->MinAngle+LIDAR->Resolution/2,2*M_PI))/LIDAR->Resolution));
            R[index] = tempR[i];
        }
        LIDAR->RangeMeasurements.assign(R.begin(), R.end());
    }
    
    return;
}


void constructLIDAR2D(const sensor_msgs::LaserScan::ConstPtr& DataLIDAR, double CutoffRange, double AllowableRange, double Pitch, LIDARClass *LIDAR) {
    /** 
     * Function that constructs a LIDAR object to be used by the reactive planner
     * 
     * Input:
     *  1) DataLIDAR: Received LIDAR data
     *  2) CutoffRange: Cutoff range below which the range measurement is ignored
     *  3) AllowableRange: Maximum allowable LIDAR range
     *  4) Pitch: Robot pitch to be considered for range compensation (default is 0)
     *  5) LIDAR: Constructed LIDAR object
     * 
     */
    
    // LIDAR operations
    double Range = AllowableRange;
    double Infinity = AllowableRange+10.;
    double MinAngle = double(DataLIDAR->angle_min);
    double MaxAngle = double(DataLIDAR->angle_max);
    double Resolution = double(DataLIDAR->angle_increment);
    double NumSample =  DataLIDAR->ranges.size();
    std::vector<double> Angle = linspace(MinAngle, MaxAngle, NumSample);
	Resolution = Angle[1]-Angle[0];
    std::vector<double> RangeMeasurements(DataLIDAR->ranges.begin(), DataLIDAR->ranges.end());

    for (size_t i = 0; i < RangeMeasurements.size() ; i++) {
        // Project on the robot plane
        RangeMeasurements[i] = RangeMeasurements[i] * cos(Pitch);

        // Reject NaNs
        if (isnan(RangeMeasurements[i])) {
            RangeMeasurements[i] = Infinity;
        }

        // Apply cutoff range
        if (RangeMeasurements[i] <= CutoffRange) {
            RangeMeasurements[i] = Infinity;
        }
    }
    
    // Construct LIDAR object
    LIDAR->RangeMeasurements = RangeMeasurements;
    LIDAR->Range = Range;
    LIDAR->Infinity = Infinity;
    LIDAR->MinAngle = MinAngle;
    LIDAR->MaxAngle = MaxAngle;
    LIDAR->Resolution = Resolution;
    LIDAR->NumSample = NumSample;
    LIDAR->Angle = Angle;
    
    return;
}


std::vector<point> obstaclePointsLIDAR2D(point RobotPosition, double RobotOrientation, LIDARClass* LIDAR) {
    /** 
     * Function that returns the coordinates of observed obstacle points from the LIDAR
     * 
     * Input:
     *  1) RobotPosition: Robot position
     *  2) RobotOrientation: Robot orientation
     *  3) LIDAR: Current LIDAR object
     * 
     * Output:
     *  1) PointsAll: Coordinates of the observed obstacle points from the LIDAR
     */

    // Find observed LIDAR points coordinates
    std::vector<point> PointsAll(LIDAR->RangeMeasurements.size());
    for (size_t i = 0 ; i < LIDAR->RangeMeasurements.size() ; i++) {
        point new_point(RobotPosition.get<0>()+LIDAR->RangeMeasurements[i]*cos(LIDAR->Angle[i]+RobotOrientation), RobotPosition.get<1>()+LIDAR->RangeMeasurements[i]*sin(LIDAR->Angle[i]+RobotOrientation));
        PointsAll[i] = new_point;
    }
    
    return PointsAll;
}


void compensateObstacleLIDAR2D(point RobotPosition, double RobotOrientation, polygon Obstacle, LIDARClass* LIDAR) {
    /** 
     * Function that checks if the LIDAR hits a specific obstacle whose polygon is known
     * 
     * Input:
     *  1) RobotPosition: Robot position
     *  2) RobotOrientation: Robot orientation
     *  3) Obstacle: shapely.geometry.Polygon obstacle defining the polygonal obstacle
     *  4) LIDAR: Current LIDAR object
     */
    
    // Find the indices that correspond to a LIDAR range less than the maximum range
    for (size_t i = 0 ; i < LIDAR->RangeMeasurements.size() ; i++) {
        if (bg::within(point(RobotPosition.get<0>()+LIDAR->RangeMeasurements[i]*cos(LIDAR->Angle[i]+RobotOrientation), RobotPosition.get<1>()+LIDAR->RangeMeasurements[i]*sin(LIDAR->Angle[i]+RobotOrientation)), Obstacle)) {
            LIDAR->RangeMeasurements[i] = LIDAR->Range;
        }
    }

    return;
}


void readLIDAR2D(point RobotPosition, double RobotOrientation, std::vector<polygon> Obstacles, double Range, double MinAngle, double MaxAngle, uint16_t NumSample, LIDARClass *virtualLIDAR) {
    /** 
     * Function that generates a virtual LIDAR object, based on the position of the robot and the surrounding obstacles
     * 
     * Input:
     *  1) RobotPosition: Robot position
     *  2) RobotOrientation: Robot orientation
     *  3) Obstacles: shapely.geometry.Polygon obstacle array defining the polygonal obstacles
     *  4) Range: Range of the LIDAR object to be generated
     *  5) MinAngle: Minimum angle of the LIDAR object to be generated
     *  6) MaxAngle: Maximum angle of the LIDAR object to be generated
     *  7) NumSample: Number of samples used in the process
     *  8) virtualLIDAR: Complete virtual LIDAR object
     * 
     */

    // Create a dummy origin
    point origin(0.0, 0.0);

    // Initialize LIDAR object
    double Resolution = (MaxAngle-MinAngle)/(NumSample-1);
    std::vector<double> RangeMeasurements(NumSample, Range);
    double Infinity = Range + 20.0;
    virtualLIDAR->RangeMeasurements = RangeMeasurements;
    virtualLIDAR->Range = Range;
    virtualLIDAR->Infinity = Infinity;
    virtualLIDAR->MinAngle = MinAngle;
    virtualLIDAR->MaxAngle = MaxAngle;
    virtualLIDAR->Resolution = Resolution;

    // Translation in the physical space to bring the points to the origin
    bg::strategy::transform::translate_transformer<double, 2, 2> translate(-RobotPosition.get<0>(), -RobotPosition.get<1>());

    // Rotation matrix from the global frame to the local sensor frame (rotate by -RobotOrientation) - Boost rotates CW
    bg::strategy::transform::rotate_transformer<bg::radian, double, 2, 2> rotate(RobotOrientation);

    // Determine distance to the workspace boundary and obstacles
    for (size_t co = 0 ; co < Obstacles.size() ; co++) {
        // Obstacle in the local sensor frame
        std::vector<point> obstacle_points = Obstacles[co].outer();

        // Bring the obstacle to the sensor frame
        for (size_t i = 0; i < obstacle_points.size(); i++) {
            point output_1, output_2;
            bg::transform(obstacle_points[i], output_1, translate);
            bg::transform(output_1, output_2, rotate);
            obstacle_points[i].set<0>(output_2.get<0>());
            obstacle_points[i].set<1>(output_2.get<1>());
        }

        // Compute distance to every obstacle edge
        for (size_t cv = 0; cv < obstacle_points.size(); cv++) {
            size_t cn = (cv+1)%(obstacle_points.size());

            point vc = obstacle_points[cv];
            point vn = obstacle_points[cn];

            // Compute the distance to the origin
            double dist = std::min(bg::distance(vc, origin), bg::distance(vn, origin));
            double w = (vn.get<0>()*(vn.get<0>()-vc.get<0>()) + vn.get<1>()*(vn.get<1>()-vc.get<1>()))/pow(bg::distance(vn,vc),2);
            if (w>=0.0 && w<=1.0) {
                point vx(w*vc.get<0>() + (1-w)*vn.get<0>(), w*vc.get<1>() + (1-w)*vn.get<1>());
                dist = std::min(dist, bg::distance(vx, origin));
            }

            double ac = atan2(vc.get<1>(), vc.get<0>());
            double an = atan2(vn.get<1>(), vn.get<0>());

            bool flagDist = (dist <= virtualLIDAR->Range);
            bool flagAngle = (std::min(ac,an) <= virtualLIDAR->MaxAngle) && (std::max(ac,an) >= virtualLIDAR->MinAngle);
            // Compute LIDAR range if the obstacle segment is in the sensing region
            if (flagDist && flagAngle) {
                // Closest LIDAR ray index
                size_t I = size_t(round((std::max(std::min(ac,an),virtualLIDAR->MinAngle)-virtualLIDAR->MinAngle)/virtualLIDAR->Resolution));
                I = (I%virtualLIDAR->NumSample);

                // Compute the intersection of the LIDAR ray with the sensor footprint
                point vtemp(cos(virtualLIDAR->Angle[I]), sin(virtualLIDAR->Angle[I]));
                point vRtemp(-sin(virtualLIDAR->Angle[I]), cos(virtualLIDAR->Angle[I]));
                double w = -bg::dot_product(vn,vRtemp)/bg::dot_product(point(vc.get<0>()-vn.get<0>(), vc.get<1>()-vn.get<1>()),vRtemp);
                if (w>=0.0 && w<=1.0) {
                    point xtemp(w*vc.get<0>()+(1-w)*vn.get<0>(), w*vc.get<1>()+(1-w)*vn.get<1>());
                    if (bg::dot_product(xtemp,vtemp) >= 0.0) {
                        virtualLIDAR->RangeMeasurements[I] = std::min(virtualLIDAR->RangeMeasurements[I], bg::distance(xtemp,origin));
                    }
                }

                // Compute the intersection of adjacent LIDAR rays
                size_t J = ((I+1)%virtualLIDAR->NumSample);
                bool flagValid = true;
                while (flagValid && (J != I)) {
                    point vtemp(cos(virtualLIDAR->Angle[J]), sin(virtualLIDAR->Angle[J]));
                    point vRtemp(-sin(virtualLIDAR->Angle[J]), cos(virtualLIDAR->Angle[J]));
                    double w = -bg::dot_product(vn,vRtemp)/bg::dot_product(point(vc.get<0>()-vn.get<0>(), vc.get<1>()-vn.get<1>()),vRtemp);
                    if (w>=0.0 && w<=1.0) {
                        point xtemp(w*vc.get<0>()+(1-w)*vn.get<0>(), w*vc.get<1>()+(1-w)*vn.get<1>());
                        if (bg::dot_product(xtemp,vtemp) >= 0.0) {
                            virtualLIDAR->RangeMeasurements[J] = std::min(virtualLIDAR->RangeMeasurements[J], bg::distance(xtemp,origin));
                            J = ((J+1)%virtualLIDAR->NumSample);
                        } else {
                            flagValid = false;
                        }
                    } else {
                        flagValid = false;
                    }
                }

                J = (I-1)%virtualLIDAR->NumSample;
                flagValid = true;
                while (flagValid && (J != I)) {
                    point vtemp(cos(virtualLIDAR->Angle[J]), sin(virtualLIDAR->Angle[J]));
                    point vRtemp(-sin(virtualLIDAR->Angle[J]), cos(virtualLIDAR->Angle[J]));
                    double w = -bg::dot_product(vn,vRtemp)/bg::dot_product(point(vc.get<0>()-vn.get<0>(), vc.get<1>()-vn.get<1>()),vRtemp);
                    if (w>=0.0 && w<=1.0) {
                        point xtemp(w*vc.get<0>()+(1-w)*vn.get<0>(), w*vc.get<1>()+(1-w)*vn.get<1>());
                        if (bg::dot_product(xtemp,vtemp) >= 0.0) {
                            virtualLIDAR->RangeMeasurements[J] = std::min(virtualLIDAR->RangeMeasurements[J], bg::distance(xtemp,origin));
                            J = ((J-1)%virtualLIDAR->NumSample);
                        } else {
                            flagValid = false;
                        }
                    } else {
                        flagValid = false;
                    }
                }
            }
        }
    }
    
    return;
}


void translateLIDAR2D(point RobotPosition, double RobotOrientation, point RobotPositionTransformed, double RobotOrientationTransformed, double RobotRadius, LIDARClass* LIDAR) {
    /** 
     * Rebase LIDAR readings from RobotState to RobotStateTransformed
     * 
     * Input:
     *  1) RobotPosition: Robot position
     *  2) RobotOrientation: Robot orientation
     *  3) RobotPositionTransformed: Transformed robot position
     *  4) RobotOrientationTransformed: Transformed robot orientation
     *  5) RobotRadius: Robot radius
     *  6) LIDAR: Current LIDAR object
     * 
     */

    // Create a dummy origin
    point origin(0.0, 0.0);

    // Account for the robot radius
    std::transform(LIDAR->RangeMeasurements.begin(), LIDAR->RangeMeasurements.end(), LIDAR->RangeMeasurements.begin(), std::bind2nd(std::minus<double>(), RobotRadius));
    
    // Find obstacle points
    std::vector<point> obstacle_points = obstaclePointsLIDAR2D(RobotPosition, RobotOrientation, LIDAR);

    // Translation in the physical space to bring the points to the origin
    bg::strategy::transform::translate_transformer<double, 2, 2> translate(-RobotPosition.get<0>(), -RobotPosition.get<1>());

    // Rotation matrix from the global frame to the local sensor frame (rotate by -RobotOrientation) - Boost rotates CW
    bg::strategy::transform::rotate_transformer<bg::radian, double, 2, 2> rotate(RobotOrientation);

    // Rotation matrix from the local sensor frame to the global frame in the transformed space - Boost rotates CW
    bg::strategy::transform::rotate_transformer<bg::radian, double, 2, 2> rotate_model(-RobotOrientationTransformed);

    // Apply rotations and translations to find the equivalent distance
    std::vector<double> R(LIDAR->RangeMeasurements.size(), LIDAR->Range);
    for (size_t i = 0 ; i < obstacle_points.size() ; i++) {
        point output_1, output_2, output_3;
        bg::transform(obstacle_points[i], output_1, translate);
        bg::transform(output_1, output_2, rotate);
        bg::transform(output_2, output_3, rotate_model);
        R[i] = bg::distance(output_3, origin);
    }

    // Update LIDAR object
    LIDAR->RangeMeasurements = R;

    return;
}


std::vector<bool> localminLIDAR2D(LIDARClass* LIDAR) {
    /**
     * Function that finds the indices of local minima of the LIDAR range data
     * 
     * Input:
     *  1) LIDAR: Current LIDAR object
     * 
     * Output:
     *  1) Imin: Indices of local minima of the LIDAR range data
     */

    // Vectors used for comparison
    std::vector<double> Rp;
    std::vector<double> Rn;

    // Compute the indices of strictly local minima of the LIDAR range data
    if ((LIDAR->MaxAngle-LIDAR->MinAngle)<2*M_PI) {
        // Assume that the region outside the angular range of LIDAR is free
        Rp.push_back(LIDAR->Range);
        Rp.insert(Rp.begin()+1, LIDAR->RangeMeasurements.begin(), LIDAR->RangeMeasurements.end()-1);
        Rn.insert(Rn.begin(), LIDAR->RangeMeasurements.begin()+1, LIDAR->RangeMeasurements.end());
        Rn.push_back(LIDAR->Range);
    } else {
        Rp.push_back(LIDAR->RangeMeasurements[LIDAR->RangeMeasurements.size()-2]);
        Rp.insert(Rp.begin()+1, LIDAR->RangeMeasurements.begin(), LIDAR->RangeMeasurements.end()-1);
        Rn.insert(Rn.begin(), LIDAR->RangeMeasurements.begin()+1, LIDAR->RangeMeasurements.end());
        Rn.push_back(LIDAR->RangeMeasurements[1]);
    }

    // Logical tests
    std::vector<bool> Imin(LIDAR->RangeMeasurements.size(), false);
    for (size_t i = 0; i < LIDAR->RangeMeasurements.size(); i++) {
        if (((LIDAR->RangeMeasurements[i]<=Rp[i]) && (LIDAR->RangeMeasurements[i]<Rn[i])) || ((LIDAR->RangeMeasurements[i]<Rp[i]) && (LIDAR->RangeMeasurements[i]<=Rn[i]))) {
            Imin[i] = true;
        }
    }
    return Imin;
}


polygon localworkspaceLIDAR2D(point RobotPosition, double RobotOrientation, double RobotRadius, LIDARClass* LIDAR) {
    /**
     * Function that returns the local workspace
     * 
     * Input:
     *  1) RobotPosition: Robot position
     *  2) RobotOrientation: Robot orientation
     *  3) RobotRadius: Robot radius
     *  4) LIDAR: Current LIDAR object, assumed to be in the (-pi, pi) range
     * 
     * Output:
     *  1) LW: Local workspace polygon
     */

    // Initialize vector of ranges
    double epsilon = 0.000000001;
    std::vector<double> R(LIDAR->RangeMeasurements.size(), LIDAR->Range);
    polygon LW;
    if (*std::min_element(LIDAR->RangeMeasurements.begin(), LIDAR->RangeMeasurements.end()) < epsilon) {
        return LW;
    } else {
        // Modify range data defining the local workspace and the local footprint
        std::vector<point> ListOfPointsFootprint(R.size());
        for (size_t i = 0; i < R.size(); i++) {
            R[i] = 0.5*(LIDAR->RangeMeasurements[i]+RobotRadius);
            ListOfPointsFootprint[i] = point(R[i]*cos(LIDAR->Angle[i]+RobotOrientation), R[i]*sin(LIDAR->Angle[i]+RobotOrientation));
        }

        // Initialize the local workspace with the minimum square that respects the LIDAR's sensing range
        std::vector<point> ListOfPoints = {point(-10,-10), point(10,-10), point(10,10), point(-10,10), point(-10,-10)};
        LW = BoostPointToBoostPoly(ListOfPoints);
        std::vector<bool> Imin = localminLIDAR2D(LIDAR);

        for (size_t k = 0; k < Imin.size(); k++) {
            if (Imin[k]) {
                if (bg::area(LW) < 0.01) {
                    return LW;
                } else {
                    // Local minimum parameters
                    double Ak = LIDAR->Angle[k]; // Angle
                    double Rk = R[k]; // Range

                    // Separating hyperplane parameters
                    point n(-cos(Ak+RobotOrientation), -sin(Ak+RobotOrientation)); // Hyperplane normal
                    point m(-Rk*n.get<0>(), -Rk*n.get<1>()); // A point on the separating hyperplane

                    // Update the local workspace by taking its intersection with the associated halfplane
                    LW = cvxpolyxhplane(LW, m, n);
                }
            } else {
                continue;
            }
        }

        // Local workspace footprint
        polygon LocalFootprint = BoostPointToBoostPoly(ListOfPointsFootprint);

        // Update local workspace
        if (bg::is_valid(LW) && bg::is_valid(LocalFootprint)) {
            std::deque<polygon> intersection_output;
            bg::intersection(LW, LocalFootprint, intersection_output);
            polygon LW_out;
            if (intersection_output.empty()) {
                return LW_out;
            } else {
                std::vector<std::vector<double>> transform = {{1, 0, RobotPosition.get<0>()}, {0, 1, RobotPosition.get<1>()}, {0, 0, 1}};
                bg::strategy::transform::matrix_transformer<double, 2, 2> transform_strategy(transform[0][0], transform[0][1], transform[0][2],
                                                                                             transform[1][0], transform[1][1], transform[1][2],
                                                                                             transform[2][0], transform[2][1], transform[2][2]);
                bg::transform(intersection_output[0], LW_out, transform_strategy);

                // Make local workspace convex
                polygon LW_final;
                bg::convex_hull(LW_out, LW_final);

                return LW_final;
            }            
        }
    }
    return LW;
}


polygon localfreespaceLIDAR2D(point RobotPosition, double RobotOrientation, double RobotRadius, LIDARClass* LIDAR) {
    /**
     * Function that returns the local workspace
     * 
     * Input:
     *  1) RobotPosition: Robot position
     *  2) RobotOrientation: Robot orientation
     *  3) RobotRadius: Robot radius
     *  4) LIDAR: Current LIDAR object - assumed to be in the (-pi, pi) range
     * 
     * Output:
     *  1) LF: Local freespace polygon
     */

    // Initialize vector of ranges
    double epsilon = 0.000000001;
    std::vector<double> R = LIDAR->RangeMeasurements;
    polygon LF;
    if (*std::min_element(LIDAR->RangeMeasurements.begin(), LIDAR->RangeMeasurements.end()) < epsilon) {
        return LF;
    } else {
        // Modify range data defining the local workspace and the local footprint
        std::vector<point> ListOfPointsFootprint(R.size());
        for (size_t i = 0; i < R.size(); i++) {
            R[i] = 0.5*(LIDAR->RangeMeasurements[i]-RobotRadius);
            ListOfPointsFootprint[i] = point(R[i]*cos(LIDAR->Angle[i]+RobotOrientation), R[i]*sin(LIDAR->Angle[i]+RobotOrientation));
        }

        // Initialize the local workspace with the minimum square that respects the LIDAR's sensing range
        std::vector<point> ListOfPoints = {point(-10.0,-10.0), point(10.0,-10.0), point(10.0,10.0), point(-10.0,10.0), point(-10.0,-10.0)};
        LF = BoostPointToBoostPoly(ListOfPoints);
        std::vector<bool> Imin = localminLIDAR2D(LIDAR);

        for (size_t k = 0; k < Imin.size(); k++) {
            if (Imin[k]) {
                if (bg::area(LF) < 0.01) {
                    return LF;
                } else {
                    // Local minimum parameters
                    double Ak = LIDAR->Angle[k]; // Angle
                    double Rk = R[k]; // Range

                    // Separating hyperplane parameters
                    point n(-cos(Ak+RobotOrientation), -sin(Ak+RobotOrientation)); // Hyperplane normal
                    point m(-Rk*n.get<0>(), -Rk*n.get<1>()); // A point on the separating hyperplane

                    // Update the local workspace by taking its intersection with the associated halfplane
                    LF = cvxpolyxhplane(LF, m, n);
                }
            } else {
                continue;
            }
        }

        // Local workspace footprint
        polygon LocalFootprint = BoostPointToBoostPoly(ListOfPointsFootprint);

        // Update local workspace
        if (bg::is_valid(LF) && bg::is_valid(LocalFootprint)) {
            std::deque<polygon> intersection_output;
            bg::intersection(LF, LocalFootprint, intersection_output);
            polygon LF_out;
            if (intersection_output.empty()) {
                return LF_out;
            } else {
                std::vector<std::vector<double>> transform = {{1, 0, RobotPosition.get<0>()}, {0, 1, RobotPosition.get<1>()}, {0, 0, 1}};
                bg::strategy::transform::matrix_transformer<double, 2, 2> transform_strategy(transform[0][0], transform[0][1], transform[0][2],
                                                                                             transform[1][0], transform[1][1], transform[1][2],
                                                                                             transform[2][0], transform[2][1], transform[2][2]);
                bg::transform(intersection_output[0], LF_out, transform_strategy);

                // Make local workspace convex
                polygon LF_final;
                bg::convex_hull(LF_out, LF_final);

                return LF_final;
            }            
        }
    }
    return LF;
}


line localfreespace_linearLIDAR2D(point RobotPosition, double RobotOrientation, polygon LF) {
    /** 
     * Function that returns the linear local freespace as the intersection of the local freespace with the current heading line
     * 
     * Input:
     *  1) RobotPosition: Robot position
     *  2) RobotOrientation: Robot orientation
     *  3) LF: Local freespace
     * 
     * Output:
     *  1) LFL: Linear freespace
     */

    // Initialize line
    line LFL;
    if (bg::area(LF) < 0.01) {
        LFL.push_back(RobotPosition);
        LFL.push_back(RobotPosition);
        return LFL;
    } else {
        LFL.push_back(RobotPosition);
        LFL.push_back(polyxray(LF, RobotPosition, point(cos(RobotOrientation), sin(RobotOrientation))));
        return LFL;
    }
}


line localfreespace_angularLIDAR2D(point RobotPosition, double RobotOrientation, polygon LF, point Goal) {
    /** 
     * Function that returns the angular local freespace as the intersection of the local freespace with the line connecting the robot to the goal
     * 
     * Input:
     *  1) RobotPosition: Robot position
     *  2) RobotOrientation: Robot orientation
     *  3) LF: Local freespace
     *  4) Goal: Designated goal
     * 
     * Output:
     *  1) LFA: Angular freespace
     */

    // Initialize line
    line LFA;
    if (bg::area(LF) < 0.01) {
        LFA.push_back(RobotPosition);
        LFA.push_back(RobotPosition);
        return LFA;
    } else {
        LFA.push_back(RobotPosition);
        LFA.push_back(Goal);
        std::deque<line> LFA_deque;
        line LFA_final;
        bg::intersection(LFA, LF, LFA_deque);
        if (LFA_deque.size() > 0) {
            LFA_final = LFA_deque[0];
        }
        return LFA_final;
    }
}


point localgoalLIDAR2D(polygon LF, point Goal) {
    /** 
     * Function that computes the local goal as the projection of the global goal on the local freespace
     * 
     * Input:
     *  1) LF: Local freespace
     *  2) Goal: Designated goal
     * 
     * Output:
     *  1) LG: Local goal
     */

    // Compute local goal - the closest point of local free space to the global goal
    if (bg::area(LF) < 0.01) {
        return Goal;
    } else {
        if (bg::within(Goal, LF)) {
            return Goal;
        } else {
            ProjectionResultStruct projection = polydist(LF, Goal);
            return projection.projected_point;
        }
    }
}


point localgoal_linearLIDAR2D(point RobotPosition, double RobotOrientation, polygon LF, point Goal) {
    /** 
     * Function that computes the linear local goal as the projection of the global goal on the linear local freespace
     * 
     * Input:
     *  1) RobotPosition: Robot position
     *  2) RobotOrientation: Robot orientation
     *  3) LF: Local freespace
     *  4) Goal: Global goal
     * 
     * Output:
     *  1) LG: Local linear goal
     */

    // Compute linear local free space
    line LFL = localfreespace_linearLIDAR2D(RobotPosition, RobotOrientation, LF);

    // Compute local goal for unicycle
    if (bg::is_empty(LFL)) {
        return RobotPosition;
    } else {
        ProjectionResultStruct projection = linedist(LFL, Goal);
        return projection.projected_point;
    }

}


point localgoal_angularLIDAR2D(point RobotPosition, double RobotOrientation, polygon LF, point Goal) {
    /** 
     * Function that computes the linear local goal as the projection of the global goal on the linear local freespace
     * 
     * Input:
     *  1) RobotPosition: Robot position
     *  2) RobotOrientation: Robot orientation
     *  3) LF: Local freespace
     *  4) Goal: Global goal
     * 
     * Output:
     *  1) LG: Local angular goal
     */

    // Compute angular local free space
    line LFA = localfreespace_angularLIDAR2D(RobotPosition, RobotOrientation, LF, Goal);

    // Compute local goal for unicycle
    if (bg::is_empty(LFA)) {
        return Goal;
    } else {
        ProjectionResultStruct projection = linedist(LFA, Goal);
        return projection.projected_point;
    }
}

void diffeoTreeTriangulation(std::vector<std::vector<double>> PolygonVertices, DiffeoParamsClass DiffeoParams, std::vector<TriangleClass> *tree) {
    /**
     * Function that calculates the triangulation tree of a polygon and augments it with properties used in semantic navigation (based on the ear clipping method)
     * 
     * Input:
     *  1) PolygonVertices: Vertex Coordinates of input polygon (start and end vertices must be the same)
     *  2) DiffeoParams: Options for the diffeomorphism construction
     *  3) tree: Stack of triangles with all the desired properties
     * 
     */

    // Create a dummy origin
    point origin(0.0, 0.0);

    // Unpack diffeomorphism parameters
    double varepsilon = DiffeoParams.get_varepsilon();
    std::vector<std::vector<double>> workspace = DiffeoParams.get_workspace();

    // Construct a polygon based on the input vertices - TODO: Check that the polygon is closed
    polygon PolygonIn = BoostPointToBoostPoly(StdToBoostPoint(PolygonVertices));

    // Construct a line and a polygon based on the workspace
    line workspaceLine = BoostPointToBoostLine(StdToBoostPoint(workspace));
    polygon workspacePolygon = BoostPointToBoostPoly(StdToBoostPoint(workspace));

    // Check if the polygon intersects the workspace boundary
    if (bg::intersects(PolygonIn, workspaceLine)) {
        // Compute the intersection with the workspace
        std::deque<polygon> polygon_to_use;
        bg::intersection(PolygonIn, workspacePolygon, polygon_to_use);

        // Find the vertices of the polygon
        std::vector<std::vector<double>> PolygonVertexList = BoostPointToStd(BoostPolyToBoostPoint(polygon_to_use[0]));

        // Compute the triangulation tree of the polygon with its dual (adjacency) graph
        polytriangulation(PolygonVertexList, workspace, true, tree);

        // Find the center and the adjacency edge to the boundary
        std::vector<point> last_triangle_vertices = tree->back().get_vertices();
        std::vector<double> dist_vector(3, 0.0);
        for (size_t i = 0; i < 3; i++) {
            dist_vector[i] = bg::distance(last_triangle_vertices[i], workspaceLine);
        }
        size_t min_dist_element = std::distance(dist_vector.begin(), std::min_element(dist_vector.begin(), dist_vector.end()));
        point median_point;
        switch(min_dist_element) {
            case 0:
                if(dist_vector[1] <= dist_vector[2]) {
                    tree->back().set_vertices({last_triangle_vertices[0], last_triangle_vertices[1], last_triangle_vertices[2]});
                    tree->back().set_adj_edge({last_triangle_vertices[0], last_triangle_vertices[1]});
                    median_point.set<0>(0.5*(last_triangle_vertices[0].get<0>()+last_triangle_vertices[1].get<0>()));
                    median_point.set<1>(0.5*(last_triangle_vertices[0].get<1>()+last_triangle_vertices[1].get<1>()));
                } else {
                    tree->back().set_vertices({last_triangle_vertices[2], last_triangle_vertices[0], last_triangle_vertices[1]});
                    tree->back().set_adj_edge({last_triangle_vertices[2], last_triangle_vertices[0]});
                    median_point.set<0>(0.5*(last_triangle_vertices[2].get<0>()+last_triangle_vertices[0].get<0>()));
                    median_point.set<1>(0.5*(last_triangle_vertices[2].get<1>()+last_triangle_vertices[0].get<1>()));
                }
                break;
            case 1:
                if(dist_vector[0] <= dist_vector[2]) {
                    tree->back().set_vertices({last_triangle_vertices[0], last_triangle_vertices[1], last_triangle_vertices[2]});
                    tree->back().set_adj_edge({last_triangle_vertices[0], last_triangle_vertices[1]});
                    median_point.set<0>(0.5*(last_triangle_vertices[0].get<0>()+last_triangle_vertices[1].get<0>()));
                    median_point.set<1>(0.5*(last_triangle_vertices[0].get<1>()+last_triangle_vertices[1].get<1>()));
                } else {
                    tree->back().set_vertices({last_triangle_vertices[1], last_triangle_vertices[2], last_triangle_vertices[0]});
                    tree->back().set_adj_edge({last_triangle_vertices[1], last_triangle_vertices[2]});
                    median_point.set<0>(0.5*(last_triangle_vertices[1].get<0>()+last_triangle_vertices[2].get<0>()));
                    median_point.set<1>(0.5*(last_triangle_vertices[1].get<1>()+last_triangle_vertices[2].get<1>()));
                }
                break;
            case 2:
                if(dist_vector[0] <= dist_vector[1]) {
                    tree->back().set_vertices({last_triangle_vertices[2], last_triangle_vertices[0], last_triangle_vertices[1]});
                    tree->back().set_adj_edge({last_triangle_vertices[2], last_triangle_vertices[0]});
                    median_point.set<0>(0.5*(last_triangle_vertices[2].get<0>()+last_triangle_vertices[0].get<0>()));
                    median_point.set<1>(0.5*(last_triangle_vertices[2].get<1>()+last_triangle_vertices[0].get<1>()));
                } else {
                    tree->back().set_vertices({last_triangle_vertices[1], last_triangle_vertices[2], last_triangle_vertices[0]});
                    tree->back().set_adj_edge({last_triangle_vertices[1], last_triangle_vertices[2]});
                    median_point.set<0>(0.5*(last_triangle_vertices[1].get<0>()+last_triangle_vertices[2].get<0>()));
                    median_point.set<1>(0.5*(last_triangle_vertices[1].get<1>()+last_triangle_vertices[2].get<1>()));
                }
                break;
        }
        last_triangle_vertices = tree->back().get_vertices();
        point median_ray(median_point.get<0>()-last_triangle_vertices[2].get<0>(), median_point.get<1>()-last_triangle_vertices[2].get<1>());
        median_ray.set<0>(median_ray.get<0>()/bg::distance(median_ray,origin));
        median_ray.set<1>(median_ray.get<1>()/bg::distance(median_ray,origin));
        point last_triangle_center(median_point.get<0>()+1.0*median_ray.get<0>(), median_point.get<1>()+1.0*median_ray.get<1>());
        tree->back().set_center(last_triangle_center);

        // Compute the tangent and normal vectors of the root triangle
        std::vector<point> r_t_vector = {point((last_triangle_vertices[1].get<0>()-last_triangle_vertices[0].get<0>())/bg::distance(last_triangle_vertices[0],last_triangle_vertices[1]), (last_triangle_vertices[1].get<1>()-last_triangle_vertices[0].get<1>())/bg::distance(last_triangle_vertices[0],last_triangle_vertices[1])), point((last_triangle_vertices[2].get<0>()-last_triangle_vertices[1].get<0>())/bg::distance(last_triangle_vertices[1],last_triangle_vertices[2]), (last_triangle_vertices[2].get<1>()-last_triangle_vertices[1].get<1>())/bg::distance(last_triangle_vertices[1],last_triangle_vertices[2])), point((last_triangle_vertices[0].get<0>()-last_triangle_vertices[2].get<0>())/bg::distance(last_triangle_vertices[0],last_triangle_vertices[2]), (last_triangle_vertices[0].get<1>()-last_triangle_vertices[2].get<1>())/bg::distance(last_triangle_vertices[0],last_triangle_vertices[2]))};
        std::vector<point> r_n_vector = {point(-r_t_vector[0].get<1>(), r_t_vector[0].get<0>()), point(-r_t_vector[1].get<1>(), r_t_vector[1].get<0>()), point(-r_t_vector[2].get<1>(), r_t_vector[2].get<0>())};
        tree->back().set_r_t(r_t_vector);
        tree->back().set_r_n(r_n_vector);

        // Find the remaining tangents and normals from vertices 0 and 1 to the center
        std::vector<point> r_center_t_vector = {point((last_triangle_center.get<0>()-last_triangle_vertices[0].get<0>())/bg::distance(last_triangle_vertices[0],last_triangle_center), (last_triangle_center.get<1>()-last_triangle_vertices[0].get<1>())/bg::distance(last_triangle_vertices[0],last_triangle_center)), point((last_triangle_vertices[1].get<0>()-last_triangle_center.get<0>())/bg::distance(last_triangle_vertices[1],last_triangle_center), (last_triangle_vertices[1].get<1>()-last_triangle_center.get<1>())/bg::distance(last_triangle_vertices[1],last_triangle_center))};
        std::vector<point> r_center_n_vector = {point(-r_center_t_vector[0].get<1>(), r_center_t_vector[0].get<0>()), point(-r_center_t_vector[1].get<1>(), r_center_t_vector[1].get<0>())};
        tree->back().set_r_center_t(r_center_t_vector);
        tree->back().set_r_center_n(r_center_n_vector);

        // Compute the dilated polygon and truncate it by the rays emanating from the center
        bg::strategy::buffer::distance_symmetric<double> distance_strategy(varepsilon);
        polygon last_triangle_polygon = BoostPointToBoostPoly({last_triangle_vertices[0], last_triangle_vertices[1], last_triangle_vertices[2], last_triangle_vertices[0]});
        multi_polygon last_triangle_multipolygon_dilated;
        bg::buffer(last_triangle_polygon, last_triangle_multipolygon_dilated, distance_strategy, side_strategy, join_strategy, end_strategy, circle_strategy);
        polygon last_triangle_polygon_dilated = last_triangle_multipolygon_dilated.front();
        polygon intersect_1 = cvxpolyxhplane(last_triangle_polygon_dilated, last_triangle_center, tree->back().get_r_center_n().front());
        polygon intersect_2 = cvxpolyxhplane(intersect_1, last_triangle_center, tree->back().get_r_center_n().back());

        // Compute the intersection with the workspace
        std::deque<polygon> output_1;
        std::deque<polygon> output_2;
        polygon final_polygon;
        bg::intersection(intersect_2, workspacePolygon, output_1);
        bg::union_(output_1[0], BoostPointToBoostPoly({last_triangle_center, last_triangle_vertices[1], last_triangle_vertices[2], last_triangle_vertices[0], last_triangle_center}), output_2);
        bg::simplify(output_2[0], final_polygon, 0.02);
        std::vector<point> last_triangle_vertices_tilde = BoostPolyToBoostPoint(final_polygon);
        last_triangle_vertices_tilde.pop_back();
        tree->back().set_vertices_tilde(last_triangle_vertices_tilde);

        // Find the tangent and normal vectors for the generated polygonal collar
        std::vector<point> r_tilde_t_vector, r_tilde_n_vector;
        for (size_t i = 0; i < last_triangle_vertices_tilde.size(); i++) {
            size_t j = (i+1)%(last_triangle_vertices_tilde.size());
            double dist_ij = bg::distance(last_triangle_vertices_tilde[i],last_triangle_vertices_tilde[j]);
            r_tilde_t_vector.push_back(point((last_triangle_vertices_tilde[j].get<0>()-last_triangle_vertices_tilde[i].get<0>())/dist_ij, (last_triangle_vertices_tilde[j].get<1>()-last_triangle_vertices_tilde[i].get<1>())/dist_ij));
            r_tilde_n_vector.push_back(point(-(last_triangle_vertices_tilde[j].get<1>()-last_triangle_vertices_tilde[i].get<1>())/dist_ij, (last_triangle_vertices_tilde[j].get<0>()-last_triangle_vertices_tilde[i].get<0>())/dist_ij));
        }
        tree->back().set_r_tilde_t(r_tilde_t_vector);
        tree->back().set_r_tilde_n(r_tilde_n_vector);

        // Add a dummy radius
        tree->back().set_radius(0.0);
    } else {
        // Compute the triangulation tree of the polygon with its dual (adjacency) graph
        polytriangulation(PolygonVertices, workspace, false, tree);

        // Find the center and radius of the root
        std::vector<point> last_triangle_vertices = tree->back().get_vertices();
        tree->back().set_center(point((last_triangle_vertices[0].get<0>()+last_triangle_vertices[1].get<0>()+last_triangle_vertices[2].get<0>())/3.0, (last_triangle_vertices[0].get<1>()+last_triangle_vertices[1].get<1>()+last_triangle_vertices[2].get<1>())/3.0));
        tree->back().set_radius(0.8*bg::distance(tree->back().get_center(), BoostPointToBoostLine({last_triangle_vertices[0], last_triangle_vertices[1], last_triangle_vertices[2], last_triangle_vertices[0]})));

        // Compute the tangent and normal vectors of the root triangle
        std::vector<point> r_t_vector = {point((last_triangle_vertices[1].get<0>()-last_triangle_vertices[0].get<0>())/bg::distance(last_triangle_vertices[0],last_triangle_vertices[1]), (last_triangle_vertices[1].get<1>()-last_triangle_vertices[0].get<1>())/bg::distance(last_triangle_vertices[0],last_triangle_vertices[1])), point((last_triangle_vertices[2].get<0>()-last_triangle_vertices[1].get<0>())/bg::distance(last_triangle_vertices[1],last_triangle_vertices[2]), (last_triangle_vertices[2].get<1>()-last_triangle_vertices[1].get<1>())/bg::distance(last_triangle_vertices[1],last_triangle_vertices[2])), point((last_triangle_vertices[0].get<0>()-last_triangle_vertices[2].get<0>())/bg::distance(last_triangle_vertices[0],last_triangle_vertices[2]), (last_triangle_vertices[0].get<1>()-last_triangle_vertices[2].get<1>())/bg::distance(last_triangle_vertices[0],last_triangle_vertices[2]))};
        std::vector<point> r_n_vector = {point(-r_t_vector[0].get<1>(), r_t_vector[0].get<0>()), point(-r_t_vector[1].get<1>(), r_t_vector[1].get<0>()), point(-r_t_vector[2].get<1>(), r_t_vector[2].get<0>())};
        tree->back().set_r_t(r_t_vector);
        tree->back().set_r_n(r_n_vector);

        // Find the polygonal color for the root by dilating the triangle by varepsilon
        bg::strategy::buffer::distance_symmetric<double> distance_strategy(varepsilon);
        polygon last_triangle_polygon = BoostPointToBoostPoly({last_triangle_vertices[0], last_triangle_vertices[1], last_triangle_vertices[2], last_triangle_vertices[0]});
        multi_polygon last_triangle_multipolygon_dilated;
        bg::buffer(last_triangle_polygon, last_triangle_multipolygon_dilated, distance_strategy, side_strategy, join_strategy, end_strategy, circle_strategy);
        polygon last_triangle_polygon_dilated = last_triangle_multipolygon_dilated.front();

        // Compute the intersection with the workspace
        std::deque<polygon> output_1;
        polygon final_polygon;
        bg::intersection(last_triangle_polygon_dilated, workspacePolygon, output_1);
        bg::simplify(output_1[0], final_polygon, 0.02);
        std::vector<point> last_triangle_vertices_tilde = BoostPolyToBoostPoint(final_polygon);
        last_triangle_vertices_tilde.pop_back();
        tree->back().set_vertices_tilde(last_triangle_vertices_tilde);

        // Find the tangent and normal vectors for the generated polygonal collar
        std::vector<point> r_tilde_t_vector, r_tilde_n_vector;
        for (size_t i = 0; i < last_triangle_vertices_tilde.size(); i++) {
            size_t j = (i+1)%(last_triangle_vertices_tilde.size());
            double dist_ij = bg::distance(last_triangle_vertices_tilde[i],last_triangle_vertices_tilde[j]);
            r_tilde_t_vector.push_back(point((last_triangle_vertices_tilde[j].get<0>()-last_triangle_vertices_tilde[i].get<0>())/dist_ij, (last_triangle_vertices_tilde[j].get<1>()-last_triangle_vertices_tilde[i].get<1>())/dist_ij));
            r_tilde_n_vector.push_back(point(-(last_triangle_vertices_tilde[j].get<1>()-last_triangle_vertices_tilde[i].get<1>())/dist_ij, (last_triangle_vertices_tilde[j].get<0>()-last_triangle_vertices_tilde[i].get<0>())/dist_ij));
        }
        tree->back().set_r_tilde_t(r_tilde_t_vector);
        tree->back().set_r_tilde_n(r_tilde_n_vector);
    }

    // Identify all the children properties
    for (size_t i = 0; i < tree->size()-1; i++) {
        // Compute the tangent and normal vectors of the child hyperplanes
        // r0 is always the shared edge between the parent and the child, r1 and r2 the rest in CCW order
        std::vector<point> triangle_vertices = (*tree)[i].get_vertices();
        std::vector<point> r_t_vector = {point((triangle_vertices[1].get<0>()-triangle_vertices[0].get<0>())/bg::distance(triangle_vertices[0],triangle_vertices[1]), (triangle_vertices[1].get<1>()-triangle_vertices[0].get<1>())/bg::distance(triangle_vertices[0],triangle_vertices[1])), point((triangle_vertices[2].get<0>()-triangle_vertices[1].get<0>())/bg::distance(triangle_vertices[1],triangle_vertices[2]), (triangle_vertices[2].get<1>()-triangle_vertices[1].get<1>())/bg::distance(triangle_vertices[1],triangle_vertices[2])), point((triangle_vertices[0].get<0>()-triangle_vertices[2].get<0>())/bg::distance(triangle_vertices[0],triangle_vertices[2]), (triangle_vertices[0].get<1>()-triangle_vertices[2].get<1>())/bg::distance(triangle_vertices[0],triangle_vertices[2]))};
        std::vector<point> r_n_vector = {point(-r_t_vector[0].get<1>(), r_t_vector[0].get<0>()), point(-r_t_vector[1].get<1>(), r_t_vector[1].get<0>()), point(-r_t_vector[2].get<1>(), r_t_vector[2].get<0>())};
        (*tree)[i].set_r_t(r_t_vector);
        (*tree)[i].set_r_n(r_n_vector);

        // Find the median from the 3rd point to the shared edge and from that compute the center for the purging transformation
        std::vector<point> triangle_adj_edge = (*tree)[i].get_adj_edge();
        point median_point;
        median_point.set<0>(0.5*(triangle_adj_edge[0].get<0>()+triangle_adj_edge[1].get<0>()));
        median_point.set<1>(0.5*(triangle_adj_edge[0].get<1>()+triangle_adj_edge[1].get<1>()));
        point median_ray(median_point.get<0>()-triangle_vertices[2].get<0>(), median_point.get<1>()-triangle_vertices[2].get<1>());
        median_ray.set<0>(median_ray.get<0>()/bg::distance(median_ray,origin));
        median_ray.set<1>(median_ray.get<1>()/bg::distance(median_ray,origin));
        std::vector<point> predecessor_triangle_vertices = (*tree)[(*tree)[i].get_predecessor()].get_vertices();
        point intersection_point = polyxray(BoostPointToBoostPoly({predecessor_triangle_vertices[0], predecessor_triangle_vertices[1], predecessor_triangle_vertices[2], predecessor_triangle_vertices[0]}), median_point, median_ray);
        point triangle_center(0.2*median_point.get<0>()+0.8*intersection_point.get<0>(), 0.2*median_point.get<1>()+0.8*intersection_point.get<1>());
        (*tree)[i].set_center(triangle_center);

        // Find the remaining tangents and normals from vertices 0 and 1 to the center
        std::vector<point> r_center_t_vector = {point((triangle_center.get<0>()-triangle_vertices[0].get<0>())/bg::distance(triangle_vertices[0],triangle_center), (triangle_center.get<1>()-triangle_vertices[0].get<1>())/bg::distance(triangle_vertices[0],triangle_center)), point((triangle_vertices[1].get<0>()-triangle_center.get<0>())/bg::distance(triangle_vertices[1],triangle_center), (triangle_vertices[1].get<1>()-triangle_center.get<1>())/bg::distance(triangle_vertices[1],triangle_center))};
        std::vector<point> r_center_n_vector = {point(-r_center_t_vector[0].get<1>(), r_center_t_vector[0].get<0>()), point(-r_center_t_vector[1].get<1>(), r_center_t_vector[1].get<0>())};
        (*tree)[i].set_r_center_t(r_center_t_vector);
        (*tree)[i].set_r_center_n(r_center_n_vector);

        // Compute the dilated polygon and truncate it by the rays emanating from the center
        bg::strategy::buffer::distance_symmetric<double> distance_strategy(varepsilon);
        polygon triangle_polygon = BoostPointToBoostPoly({triangle_vertices[0], triangle_vertices[1], triangle_vertices[2], triangle_vertices[0]});
        multi_polygon triangle_multipolygon_dilated;
        bg::buffer(triangle_polygon, triangle_multipolygon_dilated, distance_strategy, side_strategy, join_strategy, end_strategy, circle_strategy);
        polygon triangle_polygon_dilated = triangle_multipolygon_dilated.front();
        polygon intersect_1 = cvxpolyxhplane(triangle_polygon_dilated, triangle_center, (*tree)[i].get_r_center_n().front());
        polygon intersect_2 = cvxpolyxhplane(intersect_1, triangle_center, (*tree)[i].get_r_center_n().back());
        polygon candidate_polygon = intersect_2;

        // Check for collisions with all the triangles that will succeed i in the diffeomorphism construction except for its parent
        for (size_t j = i+1; j < tree->size(); j++) {
            if (j == (*tree)[i].get_predecessor()) {
                continue;
            } else {
                std::vector<point> triangle_to_test_vertices = (*tree)[j].get_vertices();
                polygon triangle_to_test = BoostPointToBoostPoly({triangle_to_test_vertices[0], triangle_to_test_vertices[1], triangle_to_test_vertices[2], triangle_to_test_vertices[0]});
                std::deque<polygon> difference_output;
                bg::difference(candidate_polygon, triangle_to_test, difference_output);

                // If the difference operation created a multipolygon, keep only the polygon that contains the barycenter of the extended triangle
                point point_to_consider((triangle_vertices[0].get<0>()+triangle_vertices[1].get<0>()+triangle_center.get<0>())/3.0, (triangle_vertices[0].get<1>()+triangle_vertices[1].get<1>()+triangle_center.get<1>())/3.0);
                BOOST_FOREACH(polygon const& difference_component, difference_output) {
                    if (bg::within(point_to_consider, difference_component)) {
                        candidate_polygon = difference_component;
                        break;
                    }
                }
            }
        }

        // Extract final vertices
        polygon candidate_polygon_simplified = candidate_polygon;
        bg::simplify(candidate_polygon, candidate_polygon_simplified, 0.02);
        std::vector<std::vector<double>> candidate_polygon_vertices = BoostPointToStd(BoostPolyToBoostPoint(candidate_polygon_simplified));
        candidate_polygon_vertices.pop_back();

        // Decompose the polygon into its convex pieces and find the piece that includes the barycenter of the extended triangle
        CGAL_Polygon_2 cgal_polygon;
        CGAL_Polygon_list partition_polys;
        CGAL_Traits partition_traits;
	    CGAL::set_pretty_mode(std::cout);
        std::vector<point> final_polygon_vertices;
        for (size_t k = 0; k < candidate_polygon_vertices.size(); k++) {
            cgal_polygon.push_back(CGAL_Point_2(candidate_polygon_vertices[k][0], candidate_polygon_vertices[k][1]));
        }
        if (cgal_polygon.is_simple()) {
            // std::cout << cgal_polygon << std::endl;
            CGAL::greene_approx_convex_partition_2(cgal_polygon.vertices_begin(), 
                                                cgal_polygon.vertices_end(), 
                                                std::back_inserter(partition_polys),
                                                partition_traits);
            assert(CGAL::convex_partition_is_valid_2(cgal_polygon.vertices_begin(),
                                                    cgal_polygon.vertices_end(),
                                                    partition_polys.begin(),
                                                    partition_polys.end(),
                                                    partition_traits));
            CGAL_Point_2 cgal_point_to_consider((triangle_vertices[0].get<0>()+triangle_vertices[1].get<0>()+triangle_center.get<0>())/3.0, (triangle_vertices[0].get<1>()+triangle_vertices[1].get<1>()+triangle_center.get<1>())/3.0);
            for (size_t k = 0; k < partition_polys.size(); k++) {
                auto it = std::next(partition_polys.begin(), k);
                CGAL_Polygon_2 polygon_to_consider = *it;
                final_polygon_vertices = {};
                for (size_t l = 0; l < polygon_to_consider.size(); l++) {
                    final_polygon_vertices.push_back(point(polygon_to_consider.vertex(l).x(), polygon_to_consider.vertex(l).y()));
                }
                final_polygon_vertices.push_back(final_polygon_vertices[0]);
                if (bg::within(point((triangle_vertices[0].get<0>()+triangle_vertices[1].get<0>()+triangle_center.get<0>())/3.0, (triangle_vertices[0].get<1>()+triangle_vertices[1].get<1>()+triangle_center.get<1>())/3.0), BoostPointToBoostPoly(final_polygon_vertices))) {
                    break;
                }
            }
        } else {
            final_polygon_vertices = BoostPolyToBoostPoint(intersect_2);
        }

        // Generate the outer polygonal collar
        polygon final_polygon = BoostPointToBoostPoly(final_polygon_vertices);
        std::deque<polygon> output;
        bg::intersection(final_polygon, workspacePolygon, output);
        std::vector<point> output_vertices = BoostPolyToBoostPoint(output[0]);
        output_vertices.pop_back();
        (*tree)[i].set_vertices_tilde(output_vertices);

        // Find the tangent and normal vectors for the generated polygonal collar
        std::vector<point> r_tilde_t_vector, r_tilde_n_vector;
        for (size_t k = 0; k < output_vertices.size(); k++) {
            size_t l = (k+1)%(output_vertices.size());
            double dist_kl = bg::distance(output_vertices[k],output_vertices[l]);
            r_tilde_t_vector.push_back(point((output_vertices[l].get<0>()-output_vertices[k].get<0>())/dist_kl, (output_vertices[l].get<1>()-output_vertices[k].get<1>())/dist_kl));
            r_tilde_n_vector.push_back(point(-(output_vertices[l].get<1>()-output_vertices[k].get<1>())/dist_kl, (output_vertices[l].get<0>()-output_vertices[k].get<0>())/dist_kl));
        }
        (*tree)[i].set_r_tilde_t(r_tilde_t_vector);
        (*tree)[i].set_r_tilde_n(r_tilde_n_vector);
    }

    return;
}


void diffeoTreeConvex(std::vector<std::vector<double>> PolygonVertices, DiffeoParamsClass DiffeoParams, std::vector<PolygonClass> *tree) {
    /**
     * Function that calculates the convex decomposition tree of a polygon and augments it with properties used in semantic navigation
     * 
     * Input:
     *  1) PolygonVertices: Vertex Coordinates of input polygon (start and end vertices must be the same)
     *  2) DiffeoParams: Options for the diffeomorphism construction
     *  3) tree: Stack of triangles with all the desired properties
     * 
     */

    // Create a dummy origin
    point origin(0.0, 0.0);

    // Unpack diffeomorphism parameters
    double varepsilon = DiffeoParams.get_varepsilon();
    std::vector<std::vector<double>> workspace = DiffeoParams.get_workspace();

    // Construct a polygon based on the input vertices
    polygon PolygonIn = BoostPointToBoostPoly(StdToBoostPoint(PolygonVertices));

    // Construct a line and a polygon based on the workspace
    line workspaceLine = BoostPointToBoostLine(StdToBoostPoint(workspace));
    polygon workspacePolygon = BoostPointToBoostPoly(StdToBoostPoint(workspace));

    // Check if the polygon intersects the workspace boundary
    if (bg::intersects(PolygonIn, workspaceLine)) {
        // Compute the intersection with the workspace
        std::deque<polygon> polygon_to_use;
        bg::intersection(PolygonIn, workspacePolygon, polygon_to_use);

        // Find the vertices of the polygon
        std::vector<std::vector<double>> PolygonVertexList = BoostPointToStd(BoostPolyToBoostPoint(polygon_to_use[0]));

        // Compute the convex decomposition tree of the polygon with its dual (adjacency) graph
        polyconvexdecomposition(PolygonVertexList, workspace, true, tree);

        // Find the adjacency edge to the boundary
        std::vector<point> last_polygon_vertices = tree->back().get_vertices();
        std::vector<double> dist_vector(last_polygon_vertices.size(), 0.0);
        for (size_t i = 0; i < last_polygon_vertices.size(); i++) {
            dist_vector[i] = bg::distance(last_polygon_vertices[i], workspaceLine);
        }
        size_t min_dist_element = std::distance(dist_vector.begin(), std::min_element(dist_vector.begin(), dist_vector.end()));
        std::vector<point> new_root_vertices;
        if (dist_vector[(min_dist_element+1)%last_polygon_vertices.size()] >= dist_vector[(min_dist_element-1)%last_polygon_vertices.size()]) {
            for (size_t j = 0; j < last_polygon_vertices.size(); j++) {
                new_root_vertices.push_back(point(last_polygon_vertices[(min_dist_element-1+j)%last_polygon_vertices.size()].get<0>(), last_polygon_vertices[(min_dist_element-1+j)%last_polygon_vertices.size()].get<1>()));
            }
        } else {
            for (size_t j = 0; j < last_polygon_vertices.size(); j++) {
                new_root_vertices.push_back(point(last_polygon_vertices[(min_dist_element+j)%last_polygon_vertices.size()].get<0>(), last_polygon_vertices[(min_dist_element+j)%last_polygon_vertices.size()].get<1>()));
            }
        }
        tree->back().set_vertices(new_root_vertices);
        last_polygon_vertices = tree->back().get_vertices();
        std::vector<point> last_polygon_adj_edge = {last_polygon_vertices[0], last_polygon_vertices[1]};
        tree->back().set_adj_edge(last_polygon_adj_edge);

        // Find the center of transformation
        point median_point(0.5*(last_polygon_vertices[1].get<0>()+last_polygon_vertices[0].get<0>()), 0.5*(last_polygon_vertices[1].get<1>()+last_polygon_vertices[0].get<1>()));
        point median_ray(median_point.get<0>()-last_polygon_vertices[2].get<0>(), median_point.get<1>()-last_polygon_vertices[2].get<1>());
        median_ray.set<0>(median_ray.get<0>()/bg::distance(median_ray,origin));
        median_ray.set<1>(median_ray.get<1>()/bg::distance(median_ray,origin));
        point last_polygon_center(median_point.get<0>()+0.1*median_ray.get<0>(), median_point.get<1>()+0.1*median_ray.get<1>());
        tree->back().set_center(last_polygon_center);

        // Find the tangent and normal vectors for the generated polygon
        std::vector<point> r_t_vector, r_n_vector;
        for (size_t i = 0; i < last_polygon_vertices.size(); i++) {
            size_t j = (i+1)%(last_polygon_vertices.size());
            double dist_ij = bg::distance(last_polygon_vertices[i],last_polygon_vertices[j]);
            r_t_vector.push_back(point((last_polygon_vertices[j].get<0>()-last_polygon_vertices[i].get<0>())/dist_ij, (last_polygon_vertices[j].get<1>()-last_polygon_vertices[i].get<1>())/dist_ij));
            r_n_vector.push_back(point(-(last_polygon_vertices[j].get<1>()-last_polygon_vertices[i].get<1>())/dist_ij, (last_polygon_vertices[j].get<0>()-last_polygon_vertices[i].get<0>())/dist_ij));
        }
        tree->back().set_r_t(r_t_vector);
        tree->back().set_r_n(r_n_vector);

        // Find the remaining tangents and normals from vertices 0 and 1 to the center
        std::vector<point> r_center_t_vector = {point((last_polygon_center.get<0>()-last_polygon_vertices[0].get<0>())/bg::distance(last_polygon_vertices[0],last_polygon_center), (last_polygon_center.get<1>()-last_polygon_vertices[0].get<1>())/bg::distance(last_polygon_vertices[0],last_polygon_center)), point((last_polygon_vertices[1].get<0>()-last_polygon_center.get<0>())/bg::distance(last_polygon_vertices[1],last_polygon_center), (last_polygon_vertices[1].get<1>()-last_polygon_center.get<1>())/bg::distance(last_polygon_vertices[1],last_polygon_center))};
        std::vector<point> r_center_n_vector = {point(-r_center_t_vector[0].get<1>(), r_center_t_vector[0].get<0>()), point(-r_center_t_vector[1].get<1>(), r_center_t_vector[1].get<0>())};
        tree->back().set_r_center_t(r_center_t_vector);
        tree->back().set_r_center_n(r_center_n_vector);

        // Compute the dilated polygon and truncate it by the rays emanating from the center
        bg::strategy::buffer::distance_symmetric<double> distance_strategy(varepsilon);
        std::vector<point> polygon_to_consider = last_polygon_vertices;
        polygon_to_consider.push_back(last_polygon_vertices[0]);
        polygon last_polygon = BoostPointToBoostPoly(polygon_to_consider);
        multi_polygon last_polygon_multipolygon_dilated;
        bg::buffer(last_polygon, last_polygon_multipolygon_dilated, distance_strategy, side_strategy, join_strategy, end_strategy, circle_strategy);
        polygon last_triangle_polygon_dilated = last_polygon_multipolygon_dilated.front();
        polygon intersect_1 = cvxpolyxhplane(last_triangle_polygon_dilated, last_polygon_center, tree->back().get_r_center_n().front());
        polygon intersect_2 = cvxpolyxhplane(intersect_1, last_polygon_center, tree->back().get_r_center_n().back());

        // Compute the intersection with the workspace
        std::deque<polygon> output_1, output_2;
        polygon final_polygon;
        bg::intersection(intersect_2, workspacePolygon, output_1);
        bg::union_(output_1[0], BoostPointToBoostPoly({last_polygon_center, last_polygon_vertices[1], last_polygon_vertices[2], last_polygon_vertices[0], last_polygon_center}), output_2);
        bg::simplify(output_2[0], final_polygon, 0.02);
        std::vector<point> last_polygon_vertices_tilde = BoostPolyToBoostPoint(final_polygon);
        last_polygon_vertices_tilde.pop_back();
        tree->back().set_vertices_tilde(last_polygon_vertices_tilde);

        // Find the tangent and normal vectors for the generated polygonal collar
        std::vector<point> r_tilde_t_vector, r_tilde_n_vector;
        for (size_t i = 0; i < last_polygon_vertices_tilde.size(); i++) {
            size_t j = (i+1)%(last_polygon_vertices_tilde.size());
            double dist_ij = bg::distance(last_polygon_vertices_tilde[i],last_polygon_vertices_tilde[j]);
            r_tilde_t_vector.push_back(point((last_polygon_vertices_tilde[j].get<0>()-last_polygon_vertices_tilde[i].get<0>())/dist_ij, (last_polygon_vertices_tilde[j].get<1>()-last_polygon_vertices_tilde[i].get<1>())/dist_ij));
            r_tilde_n_vector.push_back(point(-(last_polygon_vertices_tilde[j].get<1>()-last_polygon_vertices_tilde[i].get<1>())/dist_ij, (last_polygon_vertices_tilde[j].get<0>()-last_polygon_vertices_tilde[i].get<0>())/dist_ij));
        }
        tree->back().set_r_tilde_t(r_tilde_t_vector);
        tree->back().set_r_tilde_n(r_tilde_n_vector);

        // Finally, compute the augmented inner polygon that includes the center of deformation and update
        last_polygon_vertices.insert(last_polygon_vertices.begin()+1, last_polygon_center);
        r_t_vector.insert(r_t_vector.begin(), r_center_t_vector[1]);
        r_t_vector.insert(r_t_vector.begin(), r_center_t_vector[0]);
        r_n_vector.insert(r_n_vector.begin(), r_center_n_vector[1]);
        r_n_vector.insert(r_n_vector.begin(), r_center_n_vector[0]);
        tree->back().set_augmented_vertices(last_polygon_vertices);
        tree->back().set_r_t(r_t_vector);
        tree->back().set_r_n(r_n_vector);

        // Add a dummy radius
        tree->back().set_radius(0.0);
    } else {
        // Compute the convex decomposition tree of the polygon with its dual (adjacency) graph
        polyconvexdecomposition(PolygonVertices, workspace, false, tree);

        // Find the center of the root
        std::vector<point> last_polygon_vertices = tree->back().get_vertices();
        tree->back().set_augmented_vertices(last_polygon_vertices);
        double sum_x = 0.0, sum_y = 0.0;
        for (size_t i = 0; i < last_polygon_vertices.size(); i++) {
            sum_x = sum_x + last_polygon_vertices[i].get<0>();
            sum_y = sum_y + last_polygon_vertices[i].get<1>();
        }
	double polygon_vertex_size = (double)last_polygon_vertices.size();
        tree->back().set_center(point(sum_x/polygon_vertex_size, sum_y/polygon_vertex_size));

        // Find the radius of the root
        std::vector<point> polygon_to_consider = last_polygon_vertices;
        polygon_to_consider.push_back(last_polygon_vertices[0]);
        tree->back().set_radius(0.8*bg::distance(tree->back().get_center(), BoostPointToBoostLine(polygon_to_consider)));

        // Compute the tangent and normal vectors of the root polygon
        std::vector<point> r_t_vector, r_n_vector;
        for (size_t i = 0; i < last_polygon_vertices.size(); i++) {
            size_t j = (i+1)%(last_polygon_vertices.size());
            double dist_ij = bg::distance(last_polygon_vertices[i],last_polygon_vertices[j]);
            r_t_vector.push_back(point((last_polygon_vertices[j].get<0>()-last_polygon_vertices[i].get<0>())/dist_ij, (last_polygon_vertices[j].get<1>()-last_polygon_vertices[i].get<1>())/dist_ij));
            r_n_vector.push_back(point(-(last_polygon_vertices[j].get<1>()-last_polygon_vertices[i].get<1>())/dist_ij, (last_polygon_vertices[j].get<0>()-last_polygon_vertices[i].get<0>())/dist_ij));
        }
        tree->back().set_r_t(r_t_vector);
        tree->back().set_r_n(r_n_vector);

        // Find the polygonal color for the root by dilating the triangle by varepsilon
        bg::strategy::buffer::distance_symmetric<double> distance_strategy(varepsilon);
        polygon last_polygon = BoostPointToBoostPoly(polygon_to_consider);
        multi_polygon last_polygon_multipolygon_dilated;
        bg::buffer(last_polygon, last_polygon_multipolygon_dilated, distance_strategy, side_strategy, join_strategy, end_strategy, circle_strategy);
        polygon last_polygon_dilated = last_polygon_multipolygon_dilated.front();

        // Compute the intersection with the workspace
        std::deque<polygon> output_1;
        polygon final_polygon;
        bg::intersection(last_polygon_dilated, workspacePolygon, output_1);
        bg::simplify(output_1[0], final_polygon, 0.02);
        std::vector<point> last_polygon_vertices_tilde = BoostPolyToBoostPoint(final_polygon);
        last_polygon_vertices_tilde.pop_back();
        tree->back().set_vertices_tilde(last_polygon_vertices_tilde);

        // Find the tangent and normal vectors for the generated polygonal collar
        std::vector<point> r_tilde_t_vector, r_tilde_n_vector;
        for (size_t i = 0; i < last_polygon_vertices_tilde.size(); i++) {
            size_t j = (i+1)%(last_polygon_vertices_tilde.size());
            double dist_ij = bg::distance(last_polygon_vertices_tilde[i],last_polygon_vertices_tilde[j]);
            r_tilde_t_vector.push_back(point((last_polygon_vertices_tilde[j].get<0>()-last_polygon_vertices_tilde[i].get<0>())/dist_ij, (last_polygon_vertices_tilde[j].get<1>()-last_polygon_vertices_tilde[i].get<1>())/dist_ij));
            r_tilde_n_vector.push_back(point(-(last_polygon_vertices_tilde[j].get<1>()-last_polygon_vertices_tilde[i].get<1>())/dist_ij, (last_polygon_vertices_tilde[j].get<0>()-last_polygon_vertices_tilde[i].get<0>())/dist_ij));
        }
        tree->back().set_r_tilde_t(r_tilde_t_vector);
        tree->back().set_r_tilde_n(r_tilde_n_vector);
    }

    // Identify all the children properties
    for (size_t i = 0; i < tree->size()-1; i++) {
        // Compute the tangent and normal vectors of the child hyperplanes
        // r0 is always the shared edge between the parent and the child, r1 and r2 the rest in CCW order
        std::vector<point> polygon_vertices = (*tree)[i].get_vertices();
        std::vector<point> r_t_vector, r_n_vector;
        for (size_t k = 0; k < polygon_vertices.size(); k++) {
            size_t j = (k+1)%(polygon_vertices.size());
            double dist_jk = bg::distance(polygon_vertices[k],polygon_vertices[j]);
            r_t_vector.push_back(point((polygon_vertices[j].get<0>()-polygon_vertices[k].get<0>())/dist_jk, (polygon_vertices[j].get<1>()-polygon_vertices[k].get<1>())/dist_jk));
            r_n_vector.push_back(point(-(polygon_vertices[j].get<1>()-polygon_vertices[k].get<1>())/dist_jk, (polygon_vertices[j].get<0>()-polygon_vertices[k].get<0>())/dist_jk));
        }
        (*tree)[i].set_r_t(r_t_vector);
        (*tree)[i].set_r_n(r_n_vector);

        // Find the median from the point furthest away from the adjacency edge and from that compute the center for the purging transformation
        // To compute the center, first compute the intersection of the parent polygon with the hyperplanes of the child polygon next to the adjacency edge - this defines the admissible region within which you are allowed to search for a center
        std::vector<double> dot_product_list;
        for (size_t k = 0; k < polygon_vertices.size(); k++) {
            dot_product_list.push_back((polygon_vertices[k].get<0>()-polygon_vertices[0].get<0>())*r_n_vector[0].get<0>() + (polygon_vertices[k].get<1>()-polygon_vertices[0].get<1>())*r_n_vector[0].get<1>());
        }
        size_t max_dot_product_element = std::distance(dot_product_list.begin(), std::max_element(dot_product_list.begin(), dot_product_list.end()));
        std::vector<point> polygon_adj_edge = (*tree)[i].get_adj_edge();
        point median_point(0.5*(polygon_adj_edge[0].get<0>()+polygon_adj_edge[1].get<0>()), 0.5*(polygon_adj_edge[0].get<1>()+polygon_adj_edge[1].get<1>()));
        point median_ray(median_point.get<0>()-polygon_vertices[max_dot_product_element].get<0>(), median_point.get<1>()-polygon_vertices[max_dot_product_element].get<1>());
        median_ray.set<0>(median_ray.get<0>()/bg::distance(median_ray,origin));
        median_ray.set<1>(median_ray.get<1>()/bg::distance(median_ray,origin));
        std::vector<point> predecessor_polygon_vertices = (*tree)[(*tree)[i].get_predecessor()].get_vertices();
        predecessor_polygon_vertices.push_back(predecessor_polygon_vertices[0]);
        polygon intersect_1_cvx = cvxpolyxhplane(BoostPointToBoostPoly(predecessor_polygon_vertices), polygon_adj_edge[0], r_n_vector.back());
        polygon intersect_2_cvx = cvxpolyxhplane(intersect_1_cvx, polygon_adj_edge[1], r_n_vector[1]);
        point intersection_point = polyxray(intersect_2_cvx, point(median_point.get<0>()+0.05*median_ray.get<0>(), median_point.get<1>()+0.05*median_ray.get<1>()), median_ray);
        point polygon_center(0.5*median_point.get<0>()+0.5*intersection_point.get<0>(), 0.5*median_point.get<1>()+0.5*intersection_point.get<1>());
        (*tree)[i].set_center(polygon_center);

        // Find the remaining tangents and normals from vertices 0 and 1 to the center
        std::vector<point> r_center_t_vector = {point((polygon_center.get<0>()-polygon_vertices[0].get<0>())/bg::distance(polygon_vertices[0],polygon_center), (polygon_center.get<1>()-polygon_vertices[0].get<1>())/bg::distance(polygon_vertices[0],polygon_center)), point((polygon_vertices[1].get<0>()-polygon_center.get<0>())/bg::distance(polygon_vertices[1],polygon_center), (polygon_vertices[1].get<1>()-polygon_center.get<1>())/bg::distance(polygon_vertices[1],polygon_center))};
        std::vector<point> r_center_n_vector = {point(-r_center_t_vector[0].get<1>(), r_center_t_vector[0].get<0>()), point(-r_center_t_vector[1].get<1>(), r_center_t_vector[1].get<0>())};
        (*tree)[i].set_r_center_t(r_center_t_vector);
        (*tree)[i].set_r_center_n(r_center_n_vector);

        // Compute the augmented polygon, dilate it and truncate it by the rays emanating from the center
        // Make sure that the intersection of the dilation with the convex pieces succeeding the current convex piece in the transformation does not generate a multipolygon - otherwise reduce the radius of dilation until we have a single piece
        std::vector<point> polygon_vertices_used = polygon_vertices;
        polygon_vertices_used.insert(polygon_vertices_used.begin()+1, polygon_center);
        polygon_vertices_used.push_back(polygon_vertices_used[0]);
        polygon polygon_used = BoostPointToBoostPoly(polygon_vertices_used);
        multi_polygon polygon_used_multipolygon_dilated;
        bool multiple_elements_flag = true;
        double varepsilon_used = varepsilon;
        while (multiple_elements_flag) {
            bg::strategy::buffer::distance_symmetric<double> distance_strategy(varepsilon_used);
            bg::buffer(polygon_used, polygon_used_multipolygon_dilated, distance_strategy, side_strategy, join_strategy, end_strategy, circle_strategy);
            multi_polygon temp_result, output_union;
            std::vector<point> next_polygon = (*tree)[i+1].get_vertices();
            next_polygon.push_back(next_polygon[0]);
            output_union.push_back(BoostPointToBoostPoly(next_polygon));
            for (size_t j = i+2; j < tree->size(); j++) {
                std::vector<point> polygon_to_test_vertices = (*tree)[j].get_vertices();
                polygon_to_test_vertices.push_back(polygon_to_test_vertices[0]);
                bg::union_(output_union, BoostPointToBoostPoly(polygon_to_test_vertices), temp_result);
                output_union = temp_result;
            }
            std::deque<polygon> output_intersection;
            bg::intersection(polygon_used_multipolygon_dilated.front(), output_union, output_intersection);
            if (output_intersection.size() > 1) {
                varepsilon_used = 0.5*varepsilon_used;
            } else {
                multiple_elements_flag = false;
            }
        }
        polygon polygon_used_dilated = polygon_used_multipolygon_dilated.front();
        polygon intersect_1 = cvxpolyxhplane(polygon_used_dilated, polygon_center, (*tree)[i].get_r_center_n().front());
        polygon intersect_2 = cvxpolyxhplane(intersect_1, polygon_center, (*tree)[i].get_r_center_n().back());
        polygon candidate_polygon = intersect_2;

        // Check for collisions with all the polygons that will succeed i in the diffeomorphism construction except for its parent
        for (size_t j = i+1; j < tree->size(); j++) {
            if (j == (*tree)[i].get_predecessor()) {
                continue;
            } else {
                std::vector<point> polygon_to_test_vertices = (*tree)[j].get_vertices();
                polygon_to_test_vertices.push_back(polygon_to_test_vertices[0]);
                std::deque<polygon> difference_output;
                bg::difference(candidate_polygon, BoostPointToBoostPoly(polygon_to_test_vertices), difference_output);

                // If the difference operation created a multipolygon, keep only the polygon that contains the barycenter of the extended triangle
                point point_to_consider((polygon_vertices[0].get<0>()+polygon_vertices[1].get<0>()+polygon_center.get<0>())/3.0, (polygon_vertices[0].get<1>()+polygon_vertices[1].get<1>()+polygon_center.get<1>())/3.0);
                BOOST_FOREACH(polygon const& difference_component, difference_output) {
                    if (bg::within(point_to_consider, difference_component)) {
                        candidate_polygon = difference_component;
                        break;
                    }
                }
            }
        }

        // Extract final vertices
        polygon candidate_polygon_simplified = candidate_polygon;
        bg::simplify(candidate_polygon, candidate_polygon_simplified, 0.02);
        std::vector<std::vector<double>> candidate_polygon_vertices = BoostPointToStd(BoostPolyToBoostPoint(candidate_polygon_simplified));
        candidate_polygon_vertices.pop_back();

        // Decompose the polygon into its convex pieces and find the piece that includes the barycenter of the extended triangle
        CGAL_Polygon_2 cgal_polygon;
        CGAL_Polygon_list partition_polys;
        CGAL_Traits partition_traits;
	    CGAL::set_pretty_mode(std::cout);
        std::vector<point> final_polygon_vertices;
        for (size_t k = 0; k < candidate_polygon_vertices.size(); k++) {
            cgal_polygon.push_back(CGAL_Point_2(candidate_polygon_vertices[k][0], candidate_polygon_vertices[k][1]));
        }
        if (cgal_polygon.is_simple()) {
            // std::cout << cgal_polygon << std::endl;
            CGAL::optimal_convex_partition_2(cgal_polygon.vertices_begin(), 
                                             cgal_polygon.vertices_end(), 
                                             std::back_inserter(partition_polys),
                                             partition_traits);
            assert(CGAL::convex_partition_is_valid_2(cgal_polygon.vertices_begin(),
                                                    cgal_polygon.vertices_end(),
                                                    partition_polys.begin(),
                                                    partition_polys.end(),
                                                    partition_traits));
            CGAL_Point_2 cgal_point_to_consider((polygon_vertices[0].get<0>()+polygon_vertices[1].get<0>()+polygon_center.get<0>())/3.0, (polygon_vertices[0].get<1>()+polygon_vertices[1].get<1>()+polygon_center.get<1>())/3.0);
            for (size_t k = 0; k < partition_polys.size(); k++) {
                auto it = std::next(partition_polys.begin(), k);
                CGAL_Polygon_2 polygon_to_consider = *it;
                final_polygon_vertices = {};
                for (size_t l = 0; l < polygon_to_consider.size(); l++) {
                    final_polygon_vertices.push_back(point(polygon_to_consider.vertex(l).x(), polygon_to_consider.vertex(l).y()));
                }
                final_polygon_vertices.push_back(final_polygon_vertices[0]);
                if (bg::within(point((polygon_vertices[0].get<0>()+polygon_vertices[1].get<0>()+polygon_center.get<0>())/3.0, (polygon_vertices[0].get<1>()+polygon_vertices[1].get<1>()+polygon_center.get<1>())/3.0), BoostPointToBoostPoly(final_polygon_vertices))) {
                    break;
                }
            }
        } else {
            final_polygon_vertices = BoostPolyToBoostPoint(intersect_2);
        }

        // Generate the outer polygonal collar
        polygon final_polygon = BoostPointToBoostPoly(final_polygon_vertices);
        std::deque<polygon> output;
        bg::intersection(final_polygon, workspacePolygon, output);
        std::vector<point> output_vertices = BoostPolyToBoostPoint(output[0]);
        output_vertices.pop_back();
        (*tree)[i].set_vertices_tilde(output_vertices);

        // Find the tangent and normal vectors for the generated polygonal collar
        std::vector<point> r_tilde_t_vector, r_tilde_n_vector;
        for (size_t k = 0; k < output_vertices.size(); k++) {
            size_t l = (k+1)%(output_vertices.size());
            double dist_kl = bg::distance(output_vertices[k],output_vertices[l]);
            r_tilde_t_vector.push_back(point((output_vertices[l].get<0>()-output_vertices[k].get<0>())/dist_kl, (output_vertices[l].get<1>()-output_vertices[k].get<1>())/dist_kl));
            r_tilde_n_vector.push_back(point(-(output_vertices[l].get<1>()-output_vertices[k].get<1>())/dist_kl, (output_vertices[l].get<0>()-output_vertices[k].get<0>())/dist_kl));
        }
        (*tree)[i].set_r_tilde_t(r_tilde_t_vector);
        (*tree)[i].set_r_tilde_n(r_tilde_n_vector);

        // Finally, compute the augmented inner polygon that includes the center of deformation and update
        polygon_vertices.insert(polygon_vertices.begin()+1, polygon_center);
        r_t_vector.erase(r_t_vector.begin());
        r_n_vector.erase(r_n_vector.begin());
        r_t_vector.insert(r_t_vector.begin(), r_center_t_vector[1]);
        r_t_vector.insert(r_t_vector.begin(), r_center_t_vector[0]);
        r_n_vector.insert(r_n_vector.begin(), r_center_n_vector[1]);
        r_n_vector.insert(r_n_vector.begin(), r_center_n_vector[0]);
        (*tree)[i].set_augmented_vertices(polygon_vertices);
        (*tree)[i].set_r_t(r_t_vector);
        (*tree)[i].set_r_n(r_n_vector);
    }

    return;
}


OutputStructVector polygonDiffeoTriangulation(std::vector<double> Position, std::vector<TriangleClass> DiffeoTree, DiffeoParamsClass DiffeoParams) {
    /**
     * Function that computes h(x) (i.e., the position of point x in the model layer), Dh(x) and the derivatives of Dh, when purging a specific polygon whose dual graph and diffeomorphism properties are known (based on the ear clipping method)
     * 
     * Input:
     *  1) Position: Point to consider
     *  2) DiffeoTree: Tree that contains the diffeomorphism properties for a particular polygon
     *  3) DiffeoParams: Options for the diffeomorphism construction
     * 
     * Output:
     *  1) Output: Output struct containing the value, jacobian and derivatives of the jacobian
     */

    // Begin purging process with default values
    OutputStructVector Output;
    std::vector<double> map_h = Position;
    std::vector<std::vector<double>> map_hd = {{1.0, 0.0}, {0.0, 1.0}};
    std::vector<double> map_hdd = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // Iterate through the polygon triangles
    for (size_t i = 0; i < DiffeoTree.size(); i++) {
        OutputStructVector OutputNew = triangleDiffeo(map_h, DiffeoTree[i], DiffeoParams);

        std::vector<double> map_h_new = OutputNew.Value;
        std::vector<std::vector<double>> map_hd_new = OutputNew.Jacobian;
        std::vector<double> map_hdd_new = OutputNew.JacobianD;

        double res1 = map_hd_new[0][0]*map_hdd[0] + map_hd_new[0][1]*map_hdd[4] + map_hd[0][0]*(map_hdd_new[0]*map_hd[0][0] + map_hdd_new[1]*map_hd[1][0]) + map_hd[1][0]*(map_hdd_new[2]*map_hd[0][0] + map_hdd_new[3]*map_hd[1][0]);
        double res2 = map_hd_new[0][0]*map_hdd[1] + map_hd_new[0][1]*map_hdd[5] + map_hd[0][0]*(map_hdd_new[0]*map_hd[0][1] + map_hdd_new[1]*map_hd[1][1]) + map_hd[1][0]*(map_hdd_new[2]*map_hd[0][1] + map_hdd_new[3]*map_hd[1][1]);
        double res3 = map_hd_new[0][0]*map_hdd[2] + map_hd_new[0][1]*map_hdd[6] + map_hd[0][1]*(map_hdd_new[0]*map_hd[0][0] + map_hdd_new[1]*map_hd[1][0]) + map_hd[1][1]*(map_hdd_new[2]*map_hd[0][0] + map_hdd_new[3]*map_hd[1][0]);
        double res4 = map_hd_new[0][0]*map_hdd[3] + map_hd_new[0][1]*map_hdd[7] + map_hd[0][1]*(map_hdd_new[0]*map_hd[0][1] + map_hdd_new[1]*map_hd[1][1]) + map_hd[1][1]*(map_hdd_new[2]*map_hd[0][1] + map_hdd_new[3]*map_hd[1][1]);
        double res5 = map_hd_new[1][0]*map_hdd[0] + map_hd_new[1][1]*map_hdd[4] + map_hd[0][0]*(map_hdd_new[4]*map_hd[0][0] + map_hdd_new[5]*map_hd[1][0]) + map_hd[1][0]*(map_hdd_new[6]*map_hd[0][0] + map_hdd_new[7]*map_hd[1][0]);
        double res6 = map_hd_new[1][0]*map_hdd[1] + map_hd_new[1][1]*map_hdd[5] + map_hd[0][0]*(map_hdd_new[4]*map_hd[0][1] + map_hdd_new[5]*map_hd[1][1]) + map_hd[1][0]*(map_hdd_new[6]*map_hd[0][1] + map_hdd_new[7]*map_hd[1][1]);
        double res7 = map_hd_new[1][0]*map_hdd[2] + map_hd_new[1][1]*map_hdd[6] + map_hd[0][1]*(map_hdd_new[4]*map_hd[0][0] + map_hdd_new[5]*map_hd[1][0]) + map_hd[1][1]*(map_hdd_new[6]*map_hd[0][0] + map_hdd_new[7]*map_hd[1][0]);
        double res8 = map_hd_new[1][0]*map_hdd[3] + map_hd_new[1][1]*map_hdd[7] + map_hd[0][1]*(map_hdd_new[4]*map_hd[0][1] + map_hdd_new[5]*map_hd[1][1]) + map_hd[1][1]*(map_hdd_new[6]*map_hd[0][1] + map_hdd_new[7]*map_hd[1][1]);
        map_hdd[0] = res1;
        map_hdd[1] = res2;
        map_hdd[2] = res3;
        map_hdd[3] = res4;
        map_hdd[4] = res5;
        map_hdd[5] = res6;
        map_hdd[6] = res7;
        map_hdd[7] = res8;
        
        map_hd = MatrixMatrixMultiplication(map_hd_new, map_hd);
        
        map_h = map_h_new;
    }

    // Populate output
    Output.Value = map_h;
    Output.Jacobian = map_hd;
    Output.JacobianD = map_hdd;

    return Output;
}


OutputStructVector polygonDiffeoConvex(std::vector<double> Position, std::vector<PolygonClass> DiffeoTree, DiffeoParamsClass DiffeoParams) {
    /**
     * Function that computes h(x) (i.e., the position of point x in the model layer), Dh(x) and the derivatives of Dh, when purging a specific polygon whose dual graph and diffeomorphism properties are known (based on the convex decomposition method)
     * 
     * Input:
     *  1) Position: Point to consider
     *  2) DiffeoTree: Tree that contains the diffeomorphism properties for a particular polygon
     *  3) DiffeoParams: Options for the diffeomorphism construction
     * 
     * Output:
     *  1) Output: Output struct containing the value, jacobian and derivatives of the jacobian
     */

    // Begin purging process with default values
    OutputStructVector Output;
    std::vector<double> map_h = Position;
    std::vector<std::vector<double>> map_hd = {{1.0, 0.0}, {0.0, 1.0}};
    std::vector<double> map_hdd = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // Iterate through the polygon triangles
    for (size_t i = 0; i < DiffeoTree.size(); i++) {
        OutputStructVector OutputNew = polygonDiffeo(map_h, DiffeoTree[i], DiffeoParams);

        std::vector<double> map_h_new = OutputNew.Value;
        std::vector<std::vector<double>> map_hd_new = OutputNew.Jacobian;
        std::vector<double> map_hdd_new = OutputNew.JacobianD;

        double res1 = map_hd_new[0][0]*map_hdd[0] + map_hd_new[0][1]*map_hdd[4] + map_hd[0][0]*(map_hdd_new[0]*map_hd[0][0] + map_hdd_new[1]*map_hd[1][0]) + map_hd[1][0]*(map_hdd_new[2]*map_hd[0][0] + map_hdd_new[3]*map_hd[1][0]);
        double res2 = map_hd_new[0][0]*map_hdd[1] + map_hd_new[0][1]*map_hdd[5] + map_hd[0][0]*(map_hdd_new[0]*map_hd[0][1] + map_hdd_new[1]*map_hd[1][1]) + map_hd[1][0]*(map_hdd_new[2]*map_hd[0][1] + map_hdd_new[3]*map_hd[1][1]);
        double res3 = map_hd_new[0][0]*map_hdd[2] + map_hd_new[0][1]*map_hdd[6] + map_hd[0][1]*(map_hdd_new[0]*map_hd[0][0] + map_hdd_new[1]*map_hd[1][0]) + map_hd[1][1]*(map_hdd_new[2]*map_hd[0][0] + map_hdd_new[3]*map_hd[1][0]);
        double res4 = map_hd_new[0][0]*map_hdd[3] + map_hd_new[0][1]*map_hdd[7] + map_hd[0][1]*(map_hdd_new[0]*map_hd[0][1] + map_hdd_new[1]*map_hd[1][1]) + map_hd[1][1]*(map_hdd_new[2]*map_hd[0][1] + map_hdd_new[3]*map_hd[1][1]);
        double res5 = map_hd_new[1][0]*map_hdd[0] + map_hd_new[1][1]*map_hdd[4] + map_hd[0][0]*(map_hdd_new[4]*map_hd[0][0] + map_hdd_new[5]*map_hd[1][0]) + map_hd[1][0]*(map_hdd_new[6]*map_hd[0][0] + map_hdd_new[7]*map_hd[1][0]);
        double res6 = map_hd_new[1][0]*map_hdd[1] + map_hd_new[1][1]*map_hdd[5] + map_hd[0][0]*(map_hdd_new[4]*map_hd[0][1] + map_hdd_new[5]*map_hd[1][1]) + map_hd[1][0]*(map_hdd_new[6]*map_hd[0][1] + map_hdd_new[7]*map_hd[1][1]);
        double res7 = map_hd_new[1][0]*map_hdd[2] + map_hd_new[1][1]*map_hdd[6] + map_hd[0][1]*(map_hdd_new[4]*map_hd[0][0] + map_hdd_new[5]*map_hd[1][0]) + map_hd[1][1]*(map_hdd_new[6]*map_hd[0][0] + map_hdd_new[7]*map_hd[1][0]);
        double res8 = map_hd_new[1][0]*map_hdd[3] + map_hd_new[1][1]*map_hdd[7] + map_hd[0][1]*(map_hdd_new[4]*map_hd[0][1] + map_hdd_new[5]*map_hd[1][1]) + map_hd[1][1]*(map_hdd_new[6]*map_hd[0][1] + map_hdd_new[7]*map_hd[1][1]);
        map_hdd[0] = res1;
        map_hdd[1] = res2;
        map_hdd[2] = res3;
        map_hdd[3] = res4;
        map_hdd[4] = res5;
        map_hdd[5] = res6;
        map_hdd[6] = res7;
        map_hdd[7] = res8;
        
        map_hd = MatrixMatrixMultiplication(map_hd_new, map_hd);
        
        map_h = map_h_new;
    }

    // Populate output
    Output.Value = map_h;
    Output.Jacobian = map_hd;
    Output.JacobianD = map_hdd;

    return Output;
}


OutputStructVector triangleDiffeo(std::vector<double> Position, TriangleClass Triangle, DiffeoParamsClass DiffeoParams) {
    /**
     * Function that computes h(x) (i.e., the position of point x in the model layer), Dh(x) and the derivatives of Dh, when purging a specific triangle in space
     * 
     * Input:
     *  1) Position: Point to consider
     *  2) Triangle: Description of the triangle
     *  3) DiffeoParams: Options for the diffeomorphism construction
     * 
     * Output:
     *  1) Output: Output struct containing the value, jacobian and derivatives of the jacobian
     */

    // Compute the triangle switch and its gradient
    OutputStructScalar switchOutput = triangleSwitch(Position, Triangle, DiffeoParams);
    double sigma = switchOutput.Value;
    std::vector<double> sigmad = switchOutput.Gradient;
    std::vector<std::vector<double>> sigmadd = switchOutput.Hessian;

    // Compute the triangle deforming factor
    OutputStructScalar deformingFactorOutput = triangleDeformingFactor(Position, Triangle);
    double nu = deformingFactorOutput.Value;
    std::vector<double> nud = deformingFactorOutput.Gradient;
    std::vector<std::vector<double>> nudd = deformingFactorOutput.Hessian;

    // Extract the center
    std::vector<double> center = {Triangle.get_center().get<0>(), Triangle.get_center().get<1>()};

    // Find the map and its jacobian
    std::vector<double> map_h = {sigma*(center[0] + nu*(Position[0]-center[0])) + (1-sigma)*Position[0],
                                 sigma*(center[1] + nu*(Position[1]-center[1])) + (1-sigma)*Position[1]};
    std::vector<double> PositionMinusCenter = {Position[0]-center[0], Position[1]-center[1]};
    std::vector<std::vector<double>> SigmadOuter = VectorOuterProduct(PositionMinusCenter, sigmad);
    std::vector<std::vector<double>> NudOuter = VectorOuterProduct(PositionMinusCenter, nud);
    std::vector<std::vector<double>> map_hd = {{(nu-1)*SigmadOuter[0][0] + sigma*NudOuter[0][0] + (1+sigma*(nu-1)), (nu-1)*SigmadOuter[0][1] + sigma*NudOuter[0][1]}, {(nu-1)*SigmadOuter[1][0] + sigma*NudOuter[1][0], (nu-1)*SigmadOuter[1][1] + sigma*NudOuter[1][1] + (1+sigma*(nu-1))}};

    // Find the derivatives of the jacobian
    double map_hdd_m0_r0_s0 = 2*sigma*nud[0]+2*(nu-1)*sigmad[0]+2*(Position[0]-center[0])*sigmad[0]*nud[0]+(Position[0]-center[0])*sigma*nudd[0][0]+(Position[0]-center[0])*(nu-1)*sigmadd[0][0];
    double map_hdd_m0_r0_s1 = sigma*nud[1]+(nu-1)*sigmad[1]+(Position[0]-center[0])*sigmad[1]*nud[0]+(Position[0]-center[0])*sigma*nudd[0][1]+(Position[0]-center[0])*sigmad[0]*nud[1]+(Position[0]-center[0])*(nu-1)*sigmadd[0][1];
    double map_hdd_m0_r1_s0 = sigma*nud[1]+(Position[0]-center[0])*sigmad[0]*nud[1]+(Position[0]-center[0])*sigma*nudd[0][1]+(nu-1)*sigmad[1]+(Position[0]-center[0])*sigmad[1]*nud[0]+(Position[0]-center[0])*(nu-1)*sigmadd[0][1];
    double map_hdd_m0_r1_s1 = 2*(Position[0]-center[0])*sigmad[1]*nud[1]+(Position[0]-center[0])*sigma*nudd[1][1]+(Position[0]-center[0])*(nu-1)*sigmadd[1][1];
    double map_hdd_m1_r0_s0 = 2*(Position[1]-center[1])*sigmad[0]*nud[0]+(Position[1]-center[1])*sigma*nudd[0][0]+(Position[1]-center[1])*(nu-1)*sigmadd[0][0];
    double map_hdd_m1_r0_s1 = sigma*nud[0]+(Position[1]-center[1])*sigmad[1]*nud[0]+(Position[1]-center[1])*sigma*nudd[0][1]+(nu-1)*sigmad[0]+(Position[1]-center[1])*sigmad[0]*nud[1]+(Position[1]-center[1])*(nu-1)*sigmadd[0][1];
    double map_hdd_m1_r1_s0 = sigma*nud[0]+(nu-1)*sigmad[0]+(Position[1]-center[1])*sigmad[0]*nud[1]+(Position[1]-center[1])*sigma*nudd[0][1]+(Position[1]-center[1])*sigmad[1]*nud[0]+(Position[1]-center[1])*(nu-1)*sigmadd[0][1];
    double map_hdd_m1_r1_s1 = 2*sigma*nud[1]+2*(nu-1)*sigmad[1]+2*(Position[1]-center[1])*sigmad[1]*nud[1]+(Position[1]-center[1])*sigma*nudd[1][1]+(Position[1]-center[1])*(nu-1)*sigmadd[1][1];
    
    std::vector<double> map_hdd = {map_hdd_m0_r0_s0, map_hdd_m0_r0_s1, map_hdd_m0_r1_s0, map_hdd_m0_r1_s1, map_hdd_m1_r0_s0, map_hdd_m1_r0_s1, map_hdd_m1_r1_s0, map_hdd_m1_r1_s1};

    // Construct the output
    OutputStructVector Output;
    Output.Value = map_h;
    Output.Jacobian = map_hd;
    Output.JacobianD = map_hdd;

    return Output;
}


OutputStructVector polygonDiffeo(std::vector<double> Position, PolygonClass PolygonUsed, DiffeoParamsClass DiffeoParams) {
    /**
     * Function that computes h(x) (i.e., the position of point x in the model layer), Dh(x) and the derivatives of Dh, when purging a specific polygon in space
     * 
     * Input:
     *  1) Position: Point to consider
     *  2) PolygonUsed: Description of the polygon
     *  3) DiffeoParams: Options for the diffeomorphism construction
     * 
     * Output:
     *  1) Output: Output struct containing the value, jacobian and derivatives of the jacobian
     */

    // Compute the polygon switch and its gradient
    OutputStructScalar switchOutput = polygonSwitch(Position, PolygonUsed, DiffeoParams);
    double sigma = switchOutput.Value;
    std::vector<double> sigmad = switchOutput.Gradient;
    std::vector<std::vector<double>> sigmadd = switchOutput.Hessian;

    // Compute the polygon deforming factor
    OutputStructScalar deformingFactorOutput = polygonDeformingFactor(Position, PolygonUsed);
    double nu = deformingFactorOutput.Value;
    std::vector<double> nud = deformingFactorOutput.Gradient;
    std::vector<std::vector<double>> nudd = deformingFactorOutput.Hessian;

    // Extract the center
    std::vector<double> center = {PolygonUsed.get_center().get<0>(), PolygonUsed.get_center().get<1>()};

    // Find the map and its jacobian
    std::vector<double> map_h = {sigma*(center[0] + nu*(Position[0]-center[0])) + (1-sigma)*Position[0],
                                 sigma*(center[1] + nu*(Position[1]-center[1])) + (1-sigma)*Position[1]};
    std::vector<double> PositionMinusCenter = {Position[0]-center[0], Position[1]-center[1]};
    std::vector<std::vector<double>> SigmadOuter = VectorOuterProduct(PositionMinusCenter, sigmad);
    std::vector<std::vector<double>> NudOuter = VectorOuterProduct(PositionMinusCenter, nud);
    std::vector<std::vector<double>> map_hd = {{(nu-1)*SigmadOuter[0][0] + sigma*NudOuter[0][0] + (1+sigma*(nu-1)), (nu-1)*SigmadOuter[0][1] + sigma*NudOuter[0][1]}, {(nu-1)*SigmadOuter[1][0] + sigma*NudOuter[1][0], (nu-1)*SigmadOuter[1][1] + sigma*NudOuter[1][1] + (1+sigma*(nu-1))}};

    // Find the derivatives of the jacobian
    double map_hdd_m0_r0_s0 = 2*sigma*nud[0]+2*(nu-1)*sigmad[0]+2*(Position[0]-center[0])*sigmad[0]*nud[0]+(Position[0]-center[0])*sigma*nudd[0][0]+(Position[0]-center[0])*(nu-1)*sigmadd[0][0];
    double map_hdd_m0_r0_s1 = sigma*nud[1]+(nu-1)*sigmad[1]+(Position[0]-center[0])*sigmad[1]*nud[0]+(Position[0]-center[0])*sigma*nudd[0][1]+(Position[0]-center[0])*sigmad[0]*nud[1]+(Position[0]-center[0])*(nu-1)*sigmadd[0][1];
    double map_hdd_m0_r1_s0 = sigma*nud[1]+(Position[0]-center[0])*sigmad[0]*nud[1]+(Position[0]-center[0])*sigma*nudd[0][1]+(nu-1)*sigmad[1]+(Position[0]-center[0])*sigmad[1]*nud[0]+(Position[0]-center[0])*(nu-1)*sigmadd[0][1];
    double map_hdd_m0_r1_s1 = 2*(Position[0]-center[0])*sigmad[1]*nud[1]+(Position[0]-center[0])*sigma*nudd[1][1]+(Position[0]-center[0])*(nu-1)*sigmadd[1][1];
    double map_hdd_m1_r0_s0 = 2*(Position[1]-center[1])*sigmad[0]*nud[0]+(Position[1]-center[1])*sigma*nudd[0][0]+(Position[1]-center[1])*(nu-1)*sigmadd[0][0];
    double map_hdd_m1_r0_s1 = sigma*nud[0]+(Position[1]-center[1])*sigmad[1]*nud[0]+(Position[1]-center[1])*sigma*nudd[0][1]+(nu-1)*sigmad[0]+(Position[1]-center[1])*sigmad[0]*nud[1]+(Position[1]-center[1])*(nu-1)*sigmadd[0][1];
    double map_hdd_m1_r1_s0 = sigma*nud[0]+(nu-1)*sigmad[0]+(Position[1]-center[1])*sigmad[0]*nud[1]+(Position[1]-center[1])*sigma*nudd[0][1]+(Position[1]-center[1])*sigmad[1]*nud[0]+(Position[1]-center[1])*(nu-1)*sigmadd[0][1];
    double map_hdd_m1_r1_s1 = 2*sigma*nud[1]+2*(nu-1)*sigmad[1]+2*(Position[1]-center[1])*sigmad[1]*nud[1]+(Position[1]-center[1])*sigma*nudd[1][1]+(Position[1]-center[1])*(nu-1)*sigmadd[1][1];
    
    std::vector<double> map_hdd = {map_hdd_m0_r0_s0, map_hdd_m0_r0_s1, map_hdd_m0_r1_s0, map_hdd_m0_r1_s1, map_hdd_m1_r0_s0, map_hdd_m1_r0_s1, map_hdd_m1_r1_s0, map_hdd_m1_r1_s1};

    // Construct the output
    OutputStructVector Output;
    Output.Value = map_h;
    Output.Jacobian = map_hd;
    Output.JacobianD = map_hdd;

    return Output;
}


OutputStructScalar triangleSwitch(std::vector<double> Position, TriangleClass Triangle, DiffeoParamsClass DiffeoParams) {
    /**
     * Function that computes the overall switch value, its gradient and hessian for a point x outside a triangle
     * 
     * Input:
     *  1) Position: Point to consider
     *  2) Triangle: Description of the triangle
     *  3) DiffeoParams: Options for the diffeomorphism construction
     * 
     * Output:
     *  1) Output: Output struct containing the value, gradient and hessian
     */

    // Find the separate switch values, gradients and hessians
    OutputStructScalar BetaOutput = triangleBetaSwitch(Position, Triangle, DiffeoParams);
    OutputStructScalar GammaOutput = triangleGammaSwitch(Position, Triangle, DiffeoParams);

    // Unwrap values
    double sigma_beta = BetaOutput.Value;
    std::vector<double> sigma_betad = BetaOutput.Gradient;
    std::vector<std::vector<double>> sigma_betadd = BetaOutput.Hessian;
    double sigma_gamma = GammaOutput.Value;
    std::vector<double> sigma_gammad = GammaOutput.Gradient;
    std::vector<std::vector<double>> sigma_gammadd = GammaOutput.Hessian;

    // Find the overall switch value gradient and hessian
    double sigma;
    std::vector<double> sigmad;
    std::vector<std::vector<double>> sigmadd;
    if ((sigma_beta == 1.0) && (sigma_gamma == 0.0)) {
        sigma = 1.0;
        sigmad = {0.0, 0.0};
        sigmadd = {{0.0, 0.0}, {0.0, 0.0}};
    } else {
        double nom = sigma_beta*sigma_gamma;
        double denom = sigma_beta*sigma_gamma + (1-sigma_beta);
        sigma = nom/denom;

        std::vector<double> nomd = {sigma_gamma*sigma_betad[0] + sigma_beta*sigma_gammad[0], sigma_gamma*sigma_betad[1] + sigma_beta*sigma_gammad[1]};
        std::vector<double> denomd = {sigma_gamma*sigma_betad[0] + sigma_beta*sigma_gammad[0] - sigma_betad[0], sigma_gamma*sigma_betad[1] + sigma_beta*sigma_gammad[1] - sigma_betad[1]};
        sigmad = {(1/denom)*nomd[0] - (nom/pow(denom,2))*denomd[0], (1/denom)*nomd[1] - (nom/pow(denom,2))*denomd[1]};

        std::vector<std::vector<double>> BetaGammaOuter = VectorOuterProduct(sigma_betad, sigma_gammad);
        std::vector<std::vector<double>> GammaBetaOuter = VectorOuterProduct(sigma_gammad, sigma_betad);
        std::vector<std::vector<double>> nomdd = {{sigma_gamma*sigma_betadd[0][0] + BetaGammaOuter[0][0] + GammaBetaOuter[0][0] + sigma_beta*sigma_gammadd[0][0], sigma_gamma*sigma_betadd[0][1] + BetaGammaOuter[0][1] + GammaBetaOuter[0][1] + sigma_beta*sigma_gammadd[0][1]}, {sigma_gamma*sigma_betadd[1][0] + BetaGammaOuter[1][0] + GammaBetaOuter[1][0] + sigma_beta*sigma_gammadd[1][0], sigma_gamma*sigma_betadd[1][1] + BetaGammaOuter[1][1] + GammaBetaOuter[1][1] + sigma_beta*sigma_gammadd[1][1]}};
        std::vector<std::vector<double>> denomdd = {{sigma_gamma*sigma_betadd[0][0] + BetaGammaOuter[0][0] + GammaBetaOuter[0][0] + sigma_beta*sigma_gammadd[0][0] - sigma_betadd[0][0], sigma_gamma*sigma_betadd[0][1] + BetaGammaOuter[0][1] + GammaBetaOuter[0][1] + sigma_beta*sigma_gammadd[0][1] - sigma_betadd[0][1]}, {sigma_gamma*sigma_betadd[1][0] + BetaGammaOuter[1][0] + GammaBetaOuter[1][0] + sigma_beta*sigma_gammadd[1][0] - sigma_betadd[1][0], sigma_gamma*sigma_betadd[1][1] + BetaGammaOuter[1][1] + GammaBetaOuter[1][1] + sigma_beta*sigma_gammadd[1][1] - sigma_betadd[1][1]}};
        std::vector<std::vector<double>> NomDenomOuter = VectorOuterProduct(nomd, denomd);
        std::vector<std::vector<double>> DenomNomOuter = VectorOuterProduct(denomd, nomd);
        std::vector<std::vector<double>> DenomDenomOuter = VectorOuterProduct(denomd, denomd);
        sigmadd = {{(1/denom)*nomdd[0][0] - (1/pow(denom,2))*(NomDenomOuter[0][0]+DenomNomOuter[0][0]) + 2*(nom/pow(denom,3))*DenomDenomOuter[0][0] - (nom/pow(denom,2))*denomdd[0][0], (1/denom)*nomdd[0][1] - (1/pow(denom,2))*(NomDenomOuter[0][1]+DenomNomOuter[0][1]) + 2*(nom/pow(denom,3))*DenomDenomOuter[0][1] - (nom/pow(denom,2))*denomdd[0][1]}, {(1/denom)*nomdd[1][0] - (1/pow(denom,2))*(NomDenomOuter[1][0]+DenomNomOuter[1][0]) + 2*(nom/pow(denom,3))*DenomDenomOuter[1][0] - (nom/pow(denom,2))*denomdd[1][0], (1/denom)*nomdd[1][1] - (1/pow(denom,2))*(NomDenomOuter[1][1]+DenomNomOuter[1][1]) + 2*(nom/pow(denom,3))*DenomDenomOuter[1][1] - (nom/pow(denom,2))*denomdd[1][1]}};
    }

    // Construct the output
    OutputStructScalar Output;
    Output.Value = sigma;
    Output.Gradient = sigmad;
    Output.Hessian = sigmadd;

    return Output;
}


OutputStructScalar polygonSwitch(std::vector<double> Position, PolygonClass PolygonUsed, DiffeoParamsClass DiffeoParams) {
    /**
     * Function that computes the overall switch value, its gradient and hessian for a point x outside a polygon
     * 
     * Input:
     *  1) Position: Point to consider
     *  2) PolygonUsed: Description of the polygon
     *  3) DiffeoParams: Options for the diffeomorphism construction
     * 
     * Output:
     *  1) Output: Output struct containing the value, gradient and hessian
     */

    // Find the separate switch values, gradients and hessians
    OutputStructScalar BetaOutput = polygonBetaSwitch(Position, PolygonUsed, DiffeoParams);
    OutputStructScalar GammaOutput = polygonGammaSwitch(Position, PolygonUsed, DiffeoParams);

    // Unwrap values
    double sigma_beta = BetaOutput.Value;
    std::vector<double> sigma_betad = BetaOutput.Gradient;
    std::vector<std::vector<double>> sigma_betadd = BetaOutput.Hessian;
    double sigma_gamma = GammaOutput.Value;
    std::vector<double> sigma_gammad = GammaOutput.Gradient;
    std::vector<std::vector<double>> sigma_gammadd = GammaOutput.Hessian;

    // Find the overall switch value gradient and hessian
    double sigma;
    std::vector<double> sigmad;
    std::vector<std::vector<double>> sigmadd;
    if ((sigma_beta == 1.0) && (sigma_gamma == 0.0)) {
        sigma = 1.0;
        sigmad = {0.0, 0.0};
        sigmadd = {{0.0, 0.0}, {0.0, 0.0}};
    } else {
        double nom = sigma_beta*sigma_gamma;
        double denom = sigma_beta*sigma_gamma + (1-sigma_beta);
        sigma = nom/denom;

        std::vector<double> nomd = {sigma_gamma*sigma_betad[0] + sigma_beta*sigma_gammad[0], sigma_gamma*sigma_betad[1] + sigma_beta*sigma_gammad[1]};
        std::vector<double> denomd = {sigma_gamma*sigma_betad[0] + sigma_beta*sigma_gammad[0] - sigma_betad[0], sigma_gamma*sigma_betad[1] + sigma_beta*sigma_gammad[1] - sigma_betad[1]};
        sigmad = {(1/denom)*nomd[0] - (nom/pow(denom,2))*denomd[0], (1/denom)*nomd[1] - (nom/pow(denom,2))*denomd[1]};

        std::vector<std::vector<double>> BetaGammaOuter = VectorOuterProduct(sigma_betad, sigma_gammad);
        std::vector<std::vector<double>> GammaBetaOuter = VectorOuterProduct(sigma_gammad, sigma_betad);
        std::vector<std::vector<double>> nomdd = {{sigma_gamma*sigma_betadd[0][0] + BetaGammaOuter[0][0] + GammaBetaOuter[0][0] + sigma_beta*sigma_gammadd[0][0], sigma_gamma*sigma_betadd[0][1] + BetaGammaOuter[0][1] + GammaBetaOuter[0][1] + sigma_beta*sigma_gammadd[0][1]}, {sigma_gamma*sigma_betadd[1][0] + BetaGammaOuter[1][0] + GammaBetaOuter[1][0] + sigma_beta*sigma_gammadd[1][0], sigma_gamma*sigma_betadd[1][1] + BetaGammaOuter[1][1] + GammaBetaOuter[1][1] + sigma_beta*sigma_gammadd[1][1]}};
        std::vector<std::vector<double>> denomdd = {{sigma_gamma*sigma_betadd[0][0] + BetaGammaOuter[0][0] + GammaBetaOuter[0][0] + sigma_beta*sigma_gammadd[0][0] - sigma_betadd[0][0], sigma_gamma*sigma_betadd[0][1] + BetaGammaOuter[0][1] + GammaBetaOuter[0][1] + sigma_beta*sigma_gammadd[0][1] - sigma_betadd[0][1]}, {sigma_gamma*sigma_betadd[1][0] + BetaGammaOuter[1][0] + GammaBetaOuter[1][0] + sigma_beta*sigma_gammadd[1][0] - sigma_betadd[1][0], sigma_gamma*sigma_betadd[1][1] + BetaGammaOuter[1][1] + GammaBetaOuter[1][1] + sigma_beta*sigma_gammadd[1][1] - sigma_betadd[1][1]}};
        std::vector<std::vector<double>> NomDenomOuter = VectorOuterProduct(nomd, denomd);
        std::vector<std::vector<double>> DenomNomOuter = VectorOuterProduct(denomd, nomd);
        std::vector<std::vector<double>> DenomDenomOuter = VectorOuterProduct(denomd, denomd);
        sigmadd = {{(1/denom)*nomdd[0][0] - (1/pow(denom,2))*(NomDenomOuter[0][0]+DenomNomOuter[0][0]) + 2*(nom/pow(denom,3))*DenomDenomOuter[0][0] - (nom/pow(denom,2))*denomdd[0][0], (1/denom)*nomdd[0][1] - (1/pow(denom,2))*(NomDenomOuter[0][1]+DenomNomOuter[0][1]) + 2*(nom/pow(denom,3))*DenomDenomOuter[0][1] - (nom/pow(denom,2))*denomdd[0][1]}, {(1/denom)*nomdd[1][0] - (1/pow(denom,2))*(NomDenomOuter[1][0]+DenomNomOuter[1][0]) + 2*(nom/pow(denom,3))*DenomDenomOuter[1][0] - (nom/pow(denom,2))*denomdd[1][0], (1/denom)*nomdd[1][1] - (1/pow(denom,2))*(NomDenomOuter[1][1]+DenomNomOuter[1][1]) + 2*(nom/pow(denom,3))*DenomDenomOuter[1][1] - (nom/pow(denom,2))*denomdd[1][1]}};
    }

    // Construct the output
    OutputStructScalar Output;
    Output.Value = sigma;
    Output.Gradient = sigmad;
    Output.Hessian = sigmadd;

    return Output;
}


OutputStructScalar triangleDeformingFactor(std::vector<double> Position, TriangleClass Triangle) {
    /**
     * Function that computes the value, gradient and hessian of the deforming factor for a point x outside a triangle
     * 
     * Input:
     *  1) Position: Point to consider
     *  2) Triangle: Description of the triangle
     * 
     * Output:
     *  1) Output: Output struct containing the value, gradient and hessian
     */

    // Make position a point
    point PositionPoint(Position[0], Position[1]);

    // Unwrap triangle properties
    uint16_t depth = Triangle.get_depth();
    std::vector<point> adj_edge = Triangle.get_adj_edge();
    double radius = Triangle.get_radius();
    point CenterPoint = Triangle.get_center();
    std::vector<double> Center = {CenterPoint.get<0>(), CenterPoint.get<1>()};

    // Get triangle vertices and normal vectors and convert them
    std::vector<point> TriangleVertexListPoint = Triangle.get_vertices();
    std::vector<point> NormalVectorListPoint = Triangle.get_r_n();
    std::vector<std::vector<double>> TriangleVertexList;
    std::vector<std::vector<double>> NormalVectorList;
    for (size_t i = 0; i < Triangle.get_vertices().size(); i++) {
        TriangleVertexList.push_back({TriangleVertexListPoint[i].get<0>(), TriangleVertexListPoint[i].get<1>()});
        NormalVectorList.push_back({NormalVectorListPoint[i].get<0>(), NormalVectorListPoint[i].get<1>()});
    }

    // Distinguish whether the triangle to consider is the root or some child
    double nu;
    std::vector<double> nud;
    std::vector<std::vector<double>> nudd;
    if ((depth == 0) && adj_edge.empty()) {
        nu = radius/bg::distance(PositionPoint, CenterPoint);
        nud = {-(radius/pow(bg::distance(PositionPoint, CenterPoint),3))*(Position[0]-Center[0]), -(radius/pow(bg::distance(PositionPoint, CenterPoint),3))*(Position[1]-Center[1])};
        std::vector<double> PositionMinusCenter = {Position[0]-Center[0], Position[1]-Center[1]};
        std::vector<std::vector<double>> PositionPositionOuter = VectorOuterProduct(PositionMinusCenter, PositionMinusCenter);
        nudd = {{((3*radius)/pow(bg::distance(PositionPoint, CenterPoint),5))*PositionPositionOuter[0][0] - (radius/pow(bg::distance(PositionPoint, CenterPoint),3)), ((3*radius)/pow(bg::distance(PositionPoint, CenterPoint),5))*PositionPositionOuter[0][1]}, {((3*radius)/pow(bg::distance(PositionPoint, CenterPoint),5))*PositionPositionOuter[1][0], ((3*radius)/pow(bg::distance(PositionPoint, CenterPoint),5))*PositionPositionOuter[1][1] - (radius/pow(bg::distance(PositionPoint, CenterPoint),3))}};
    } else {
        nu = ((TriangleVertexList[0][0]-Center[0])*NormalVectorList[0][0]+(TriangleVertexList[0][1]-Center[1])*NormalVectorList[0][1])/((Position[0]-Center[0])*NormalVectorList[0][0]+(Position[1]-Center[1])*NormalVectorList[0][1]);
        nud = {-(((TriangleVertexList[0][0]-Center[0])*NormalVectorList[0][0]+(TriangleVertexList[0][1]-Center[1])*NormalVectorList[0][1])/pow((Position[0]-Center[0])*NormalVectorList[0][0]+(Position[1]-Center[1])*NormalVectorList[0][1],2))*NormalVectorList[0][0], -(((TriangleVertexList[0][0]-Center[0])*NormalVectorList[0][0]+(TriangleVertexList[0][1]-Center[1])*NormalVectorList[0][1])/pow((Position[0]-Center[0])*NormalVectorList[0][0]+(Position[1]-Center[1])*NormalVectorList[0][1],2))*NormalVectorList[0][1]};
        std::vector<double> NormalVector = {NormalVectorList[0][0], NormalVectorList[0][1]};
        std::vector<std::vector<double>> NormalVectorOuter = VectorOuterProduct(NormalVector, NormalVector);
        nudd = {{2*(((TriangleVertexList[0][0]-Center[0])*NormalVectorList[0][0]+(TriangleVertexList[0][1]-Center[1])*NormalVectorList[0][1])/pow((Position[0]-Center[0])*NormalVectorList[0][0]+(Position[1]-Center[1])*NormalVectorList[0][1],3))*NormalVectorOuter[0][0], 2*(((TriangleVertexList[0][0]-Center[0])*NormalVectorList[0][0]+(TriangleVertexList[0][1]-Center[1])*NormalVectorList[0][1])/pow((Position[0]-Center[0])*NormalVectorList[0][0]+(Position[1]-Center[1])*NormalVectorList[0][1],3))*NormalVectorOuter[0][1]}, {2*(((TriangleVertexList[0][0]-Center[0])*NormalVectorList[0][0]+(TriangleVertexList[0][1]-Center[1])*NormalVectorList[0][1])/pow((Position[0]-Center[0])*NormalVectorList[0][0]+(Position[1]-Center[1])*NormalVectorList[0][1],3))*NormalVectorOuter[1][0], 2*(((TriangleVertexList[0][0]-Center[0])*NormalVectorList[0][0]+(TriangleVertexList[0][1]-Center[1])*NormalVectorList[0][1])/pow((Position[0]-Center[0])*NormalVectorList[0][0]+(Position[1]-Center[1])*NormalVectorList[0][1],3))*NormalVectorOuter[1][1]}};
    }

    // Construct the output
    OutputStructScalar Output;
    Output.Value = nu;
    Output.Gradient = nud;
    Output.Hessian = nudd;

    return Output;
}


OutputStructScalar polygonDeformingFactor(std::vector<double> Position, PolygonClass PolygonUsed) {
    /**
     * Function that computes the value, gradient and hessian of the deforming factor for a point x outside a polygon
     * 
     * Input:
     *  1) Position: Point to consider
     *  2) PolygonUsed: Description of the polygon
     * 
     * Output:
     *  1) Output: Output struct containing the value, gradient and hessian
     */

    // Make position a point
    point PositionPoint(Position[0], Position[1]);

    // Unwrap polygon properties
    uint16_t depth = PolygonUsed.get_depth();
    std::vector<point> adj_edge = PolygonUsed.get_adj_edge();
    double radius = PolygonUsed.get_radius();
    point CenterPoint = PolygonUsed.get_center();
    std::vector<double> Center = {CenterPoint.get<0>(), CenterPoint.get<1>()};

    // Get polygon vertices and normal vectors and convert them
    std::vector<point> PolygonVertexListPoint = PolygonUsed.get_augmented_vertices();
    std::vector<point> NormalVectorListPoint = PolygonUsed.get_r_n();
    std::vector<std::vector<double>> PolygonVertexList;
    for (size_t i = 0; i < PolygonUsed.get_augmented_vertices().size(); i++) {
        PolygonVertexList.push_back({PolygonVertexListPoint[i].get<0>(), PolygonVertexListPoint[i].get<1>()});
    }

    // Distinguish whether the polygon to consider is the root or some child
    double nu;
    std::vector<double> nud;
    std::vector<std::vector<double>> nudd;
    if ((depth == 0) && adj_edge.empty()) {
        nu = radius/bg::distance(PositionPoint, CenterPoint);
        nud = {-(radius/pow(bg::distance(PositionPoint, CenterPoint),3))*(Position[0]-Center[0]), -(radius/pow(bg::distance(PositionPoint, CenterPoint),3))*(Position[1]-Center[1])};
        std::vector<double> PositionMinusCenter = {Position[0]-Center[0], Position[1]-Center[1]};
        std::vector<std::vector<double>> PositionPositionOuter = VectorOuterProduct(PositionMinusCenter, PositionMinusCenter);
        nudd = {{((3*radius)/pow(bg::distance(PositionPoint, CenterPoint),5))*PositionPositionOuter[0][0] - (radius/pow(bg::distance(PositionPoint, CenterPoint),3)), ((3*radius)/pow(bg::distance(PositionPoint, CenterPoint),5))*PositionPositionOuter[0][1]}, {((3*radius)/pow(bg::distance(PositionPoint, CenterPoint),5))*PositionPositionOuter[1][0], ((3*radius)/pow(bg::distance(PositionPoint, CenterPoint),5))*PositionPositionOuter[1][1] - (radius/pow(bg::distance(PositionPoint, CenterPoint),3))}};
    } else {
	std::vector<double> adj_edge_tangent = {(adj_edge[1].get<0>()-adj_edge[0].get<0>())/bg::distance(adj_edge[0],adj_edge[1]), (adj_edge[1].get<1>()-adj_edge[0].get<1>())/bg::distance(adj_edge[0],adj_edge[1])};
	std::vector<double> adj_edge_normal = {-adj_edge_tangent[1], adj_edge_tangent[0]};
        nu = ((PolygonVertexList[0][0]-Center[0])*adj_edge_normal[0]+(PolygonVertexList[0][1]-Center[1])*adj_edge_normal[1])/((Position[0]-Center[0])*adj_edge_normal[0]+(Position[1]-Center[1])*adj_edge_normal[1]);
        nud = {-(((PolygonVertexList[0][0]-Center[0])*adj_edge_normal[0]+(PolygonVertexList[0][1]-Center[1])*adj_edge_normal[1])/pow((Position[0]-Center[0])*adj_edge_normal[0]+(Position[1]-Center[1])*adj_edge_normal[1],2))*adj_edge_normal[0], -(((PolygonVertexList[0][0]-Center[0])*adj_edge_normal[0]+(PolygonVertexList[0][1]-Center[1])*adj_edge_normal[1])/pow((Position[0]-Center[0])*adj_edge_normal[0]+(Position[1]-Center[1])*adj_edge_normal[1],2))*adj_edge_normal[1]};
        std::vector<std::vector<double>> NormalVectorOuter = VectorOuterProduct(adj_edge_normal, adj_edge_normal);
        nudd = {{2*(((PolygonVertexList[0][0]-Center[0])*adj_edge_normal[0]+(PolygonVertexList[0][1]-Center[1])*adj_edge_normal[1])/pow((Position[0]-Center[0])*adj_edge_normal[0]+(Position[1]-Center[1])*adj_edge_normal[1],3))*NormalVectorOuter[0][0], 2*(((PolygonVertexList[0][0]-Center[0])*adj_edge_normal[0]+(PolygonVertexList[0][1]-Center[1])*adj_edge_normal[1])/pow((Position[0]-Center[0])*adj_edge_normal[0]+(Position[1]-Center[1])*adj_edge_normal[1],3))*NormalVectorOuter[0][1]}, {2*(((PolygonVertexList[0][0]-Center[0])*adj_edge_normal[0]+(PolygonVertexList[0][1]-Center[1])*adj_edge_normal[1])/pow((Position[0]-Center[0])*adj_edge_normal[0]+(Position[1]-Center[1])*adj_edge_normal[1],3))*NormalVectorOuter[1][0], 2*(((PolygonVertexList[0][0]-Center[0])*adj_edge_normal[0]+(PolygonVertexList[0][1]-Center[1])*adj_edge_normal[1])/pow((Position[0]-Center[0])*adj_edge_normal[0]+(Position[1]-Center[1])*adj_edge_normal[1],3))*NormalVectorOuter[1][1]}};
    }

    // Construct the output
    OutputStructScalar Output;
    Output.Value = nu;
    Output.Gradient = nud;
    Output.Hessian = nudd;

    return Output;
}


OutputStructScalar triangleBetaSwitch(std::vector<double> Position, TriangleClass Triangle, DiffeoParamsClass DiffeoParams) {
    /**
     * Function that computes the value, gradient and hessian of the beta-switch for a point x outside a triangle
     * 
     * Input:
     *  1) Position: Point to consider
     *  2) Triangle: Description of the triangle
     *  3) DiffeoParams: Options for the diffeomorphism construction
     * 
     * Output:
     *  1) Output: Output struct containing the value, gradient and hessian
     */

    // Unwrap parameters
    double mu_1 = DiffeoParams.get_mu_1();
    double epsilon = DiffeoParams.get_epsilon();

    // Compute the value of beta and its gradient and hessian
    OutputStructScalar BetaOutput = triangleOutsideImplicit(Position, Triangle, DiffeoParams);
    double beta = BetaOutput.Value;
    std::vector<double> betad = BetaOutput.Gradient;
    std::vector<std::vector<double>> betadd = BetaOutput.Hessian;

    // Compute the value of the switch
    double sigma;
    std::vector<double> sigmad;
    std::vector<std::vector<double>> sigmadd;
    if (beta >= epsilon-1e-3) {
        sigma = 0.0;
        sigmad = {0.0, 0.0};
        sigmadd = {{0.0, 0.0}, {0.0, 0.0}};
    } else {
        sigma = exp(-mu_1/(epsilon-beta))/exp(-mu_1/epsilon);
        sigmad = {-mu_1*(sigma/pow(epsilon-beta,2))*betad[0], -mu_1*(sigma/pow(epsilon-beta,2))*betad[1]};
        std::vector<std::vector<double>> BetadOuter = VectorOuterProduct(betad, betad);
        sigmadd = {{((pow(mu_1,2)*(sigma/pow(epsilon-beta,4)))-2*mu_1*(sigma/pow(epsilon-beta,3)))*BetadOuter[0][0] - mu_1*(sigma/pow(epsilon-beta,2))*betadd[0][0], ((pow(mu_1,2)*(sigma/pow(epsilon-beta,4)))-2*mu_1*(sigma/pow(epsilon-beta,3)))*BetadOuter[0][1] - mu_1*(sigma/pow(epsilon-beta,2))*betadd[0][1]}, {((pow(mu_1,2)*(sigma/pow(epsilon-beta,4)))-2*mu_1*(sigma/pow(epsilon-beta,3)))*BetadOuter[1][0] - mu_1*(sigma/pow(epsilon-beta,2))*betadd[1][0], ((pow(mu_1,2)*(sigma/pow(epsilon-beta,4)))-2*mu_1*(sigma/pow(epsilon-beta,3)))*BetadOuter[1][1] - mu_1*(sigma/pow(epsilon-beta,2))*betadd[1][1]}};
    }

    // Construct the output
    OutputStructScalar Output;
    Output.Value = sigma;
    Output.Gradient = sigmad;
    Output.Hessian = sigmadd;

    return Output;
}


OutputStructScalar polygonBetaSwitch(std::vector<double> Position, PolygonClass PolygonUsed, DiffeoParamsClass DiffeoParams) {
    /**
     * Function that computes the value, gradient and hessian of the beta-switch for a point x outside a polygon
     * 
     * Input:
     *  1) Position: Point to consider
     *  2) PolygonUsed: Description of the polygon
     *  3) DiffeoParams: Options for the diffeomorphism construction
     * 
     * Output:
     *  1) Output: Output struct containing the value, gradient and hessian
     */

    // Unwrap parameters
    double mu_1 = DiffeoParams.get_mu_1();
    double epsilon = DiffeoParams.get_epsilon();

    // Compute the value of beta and its gradient and hessian
    OutputStructScalar BetaOutput = polygonOutsideImplicit(Position, PolygonUsed, DiffeoParams);
    double beta = BetaOutput.Value;
    std::vector<double> betad = BetaOutput.Gradient;
    std::vector<std::vector<double>> betadd = BetaOutput.Hessian;

    // Compute the value of the switch
    double sigma;
    std::vector<double> sigmad;
    std::vector<std::vector<double>> sigmadd;
    if (beta >= epsilon) {
        sigma = 0.0;
        sigmad = {0.0, 0.0};
        sigmadd = {{0.0, 0.0}, {0.0, 0.0}};
    } else {
        sigma = exp(-mu_1/(epsilon-beta))/exp(-mu_1/epsilon);
        sigmad = {-mu_1*(sigma/pow(epsilon-beta,2))*betad[0], -mu_1*(sigma/pow(epsilon-beta,2))*betad[1]};
        std::vector<std::vector<double>> BetadOuter = VectorOuterProduct(betad, betad);
        sigmadd = {{((pow(mu_1,2)*(sigma/pow(epsilon-beta,4)))-2*mu_1*(sigma/pow(epsilon-beta,3)))*BetadOuter[0][0] - mu_1*(sigma/pow(epsilon-beta,2))*betadd[0][0], ((pow(mu_1,2)*(sigma/pow(epsilon-beta,4)))-2*mu_1*(sigma/pow(epsilon-beta,3)))*BetadOuter[0][1] - mu_1*(sigma/pow(epsilon-beta,2))*betadd[0][1]}, {((pow(mu_1,2)*(sigma/pow(epsilon-beta,4)))-2*mu_1*(sigma/pow(epsilon-beta,3)))*BetadOuter[1][0] - mu_1*(sigma/pow(epsilon-beta,2))*betadd[1][0], ((pow(mu_1,2)*(sigma/pow(epsilon-beta,4)))-2*mu_1*(sigma/pow(epsilon-beta,3)))*BetadOuter[1][1] - mu_1*(sigma/pow(epsilon-beta,2))*betadd[1][1]}};
    }

    // Construct the output
    OutputStructScalar Output;
    Output.Value = sigma;
    Output.Gradient = sigmad;
    Output.Hessian = sigmadd;

    return Output;
}


OutputStructScalar triangleGammaSwitch(std::vector<double> Position, TriangleClass Triangle, DiffeoParamsClass DiffeoParams) {
    /**
     * Function that computes the value, gradient and hessian of the gamma-switch for a point x outside a triangle
     * 
     * Input:
     *  1) Position: Point to consider
     *  2) Triangle: Description of the triangle
     *  3) DiffeoParams: Options for the diffeomorphism construction
     * 
     * Output:
     *  1) Output: Output struct containing the value, gradient and hessian
     */

    // Unwrap parameters
    double mu_2 = DiffeoParams.get_mu_2();

    // Make position a point
    point PositionPoint(Position[0], Position[1]);

    // Unwrap center properties
    point CenterPoint = Triangle.get_center();
    std::vector<double> Center = {CenterPoint.get<0>(), CenterPoint.get<1>()};

    // Compute the value of gamma and its gradient and hessian
    OutputStructScalar GammaOutput = triangleInsideImplicit(Position, Triangle, DiffeoParams);
    double gamma = GammaOutput.Value;
    std::vector<double> gammad = GammaOutput.Gradient;
    std::vector<std::vector<double>> gammadd = GammaOutput.Hessian;

    // Compute the value of the switch
    double sigma;
    std::vector<double> sigmad;
    std::vector<std::vector<double>> sigmadd;
    if (gamma<=0.0) {
        sigma = 0.0;
        sigmad = {0.0, 0.0};
        sigmadd = {{0.0, 0.0}, {0.0, 0.0}};
    } else {
        // Compute the value of alpha and its gradient and hessian
        double nom = gamma;
        double denom = bg::distance(PositionPoint, CenterPoint);
        double alpha = nom/denom;

        std::vector<double> nomd = gammad;
        std::vector<double> denomd = {(1/bg::distance(PositionPoint, CenterPoint))*(Position[0]-Center[0]), (1/bg::distance(PositionPoint, CenterPoint))*(Position[1]-Center[1])};
        std::vector<double> alphad = {(1/denom)*nomd[0]-(nom/pow(denom,2))*denomd[0], (1/denom)*nomd[1]-(nom/pow(denom,2))*denomd[1]};

        std::vector<std::vector<double>> nomdd = gammadd;
        std::vector<double> PositionMinusCenter = {Position[0]-Center[0], Position[1]-Center[1]};
        std::vector<std::vector<double>> PositionPositionOuter = VectorOuterProduct(PositionMinusCenter, PositionMinusCenter);
        std::vector<std::vector<double>> denomdd = {{(1/bg::distance(PositionPoint, CenterPoint)) - (1/pow(bg::distance(PositionPoint, CenterPoint),3))*PositionPositionOuter[0][0], -(1/pow(bg::distance(PositionPoint, CenterPoint),3))*PositionPositionOuter[0][1]}, {-(1/pow(bg::distance(PositionPoint, CenterPoint),3))*PositionPositionOuter[1][0], (1/bg::distance(PositionPoint, CenterPoint)) - (1/pow(bg::distance(PositionPoint, CenterPoint),3))*PositionPositionOuter[1][1]}};
        std::vector<std::vector<double>> NomDenomOuter = VectorOuterProduct(nomd, denomd);
        std::vector<std::vector<double>> DenomNomOuter = VectorOuterProduct(denomd, nomd);
        std::vector<std::vector<double>> DenomDenomOuter = VectorOuterProduct(denomd, denomd);
        std::vector<std::vector<double>> alphadd = {{(1/denom)*nomdd[0][0] - (1/pow(denom,2))*(NomDenomOuter[0][0]+DenomNomOuter[0][0]) + 2*(nom/pow(denom,3))*DenomDenomOuter[0][0] - (nom/pow(denom,2))*denomdd[0][0], (1/denom)*nomdd[0][1] - (1/pow(denom,2))*(NomDenomOuter[0][1]+DenomNomOuter[0][1]) + 2*(nom/pow(denom,3))*DenomDenomOuter[0][1] - (nom/pow(denom,2))*denomdd[0][1]}, {(1/denom)*nomdd[1][0] - (1/pow(denom,2))*(NomDenomOuter[1][0]+DenomNomOuter[1][0]) + 2*(nom/pow(denom,3))*DenomDenomOuter[1][0] - (nom/pow(denom,2))*denomdd[1][0], (1/denom)*nomdd[1][1] - (1/pow(denom,2))*(NomDenomOuter[1][1]+DenomNomOuter[1][1]) + 2*(nom/pow(denom,3))*DenomDenomOuter[1][1] - (nom/pow(denom,2))*denomdd[1][1]}};

        sigma = exp(-mu_2/alpha);
        sigmad = {mu_2*(sigma/pow(alpha,2))*alphad[0], mu_2*(sigma/pow(alpha,2))*alphad[1]};
        std::vector<std::vector<double>> AlphadOuter = VectorOuterProduct(alphad, alphad);
        sigmadd = {{(pow(mu_2,2)*(sigma/pow(alpha,4))-2*mu_2*(sigma/pow(alpha,3)))*AlphadOuter[0][0]+mu_2*(sigma/pow(alpha,2))*alphadd[0][0], (pow(mu_2,2)*(sigma/pow(alpha,4))-2*mu_2*(sigma/pow(alpha,3)))*AlphadOuter[0][1]+mu_2*(sigma/pow(alpha,2))*alphadd[0][1]}, {(pow(mu_2,2)*(sigma/pow(alpha,4))-2*mu_2*(sigma/pow(alpha,3)))*AlphadOuter[1][0]+mu_2*(sigma/pow(alpha,2))*alphadd[1][0], (pow(mu_2,2)*(sigma/pow(alpha,4))-2*mu_2*(sigma/pow(alpha,3)))*AlphadOuter[1][1]+mu_2*(sigma/pow(alpha,2))*alphadd[1][1]}};
    }

    // Construct the output
    OutputStructScalar Output;
    Output.Value = sigma;
    Output.Gradient = sigmad;
    Output.Hessian = sigmadd;

    return Output;
}


OutputStructScalar polygonGammaSwitch(std::vector<double> Position, PolygonClass PolygonUsed, DiffeoParamsClass DiffeoParams) {
    /**
     * Function that computes the value, gradient and hessian of the gamma-switch for a point x outside a polygon
     * 
     * Input:
     *  1) Position: Point to consider
     *  2) PolygonUsed: Description of the polygon
     *  3) DiffeoParams: Options for the diffeomorphism construction
     * 
     * Output:
     *  1) Output: Output struct containing the value, gradient and hessian
     */

    // Unwrap parameters
    double mu_2 = DiffeoParams.get_mu_2();

    // Make position a point
    point PositionPoint(Position[0], Position[1]);

    // Unwrap center properties
    point CenterPoint = PolygonUsed.get_center();
    std::vector<double> Center = {CenterPoint.get<0>(), CenterPoint.get<1>()};

    // Compute the value of gamma and its gradient and hessian
    OutputStructScalar GammaOutput = polygonInsideImplicit(Position, PolygonUsed, DiffeoParams);
    double gamma = GammaOutput.Value;
    std::vector<double> gammad = GammaOutput.Gradient;
    std::vector<std::vector<double>> gammadd = GammaOutput.Hessian;

    // Compute the value of the switch
    double sigma;
    std::vector<double> sigmad;
    std::vector<std::vector<double>> sigmadd;
    if (gamma<=0.01) {
        sigma = 0.0;
        sigmad = {0.0, 0.0};
        sigmadd = {{0.0, 0.0}, {0.0, 0.0}};
    } else {
        // Compute the value of alpha and its gradient and hessian
        double nom = gamma;
        double denom = bg::distance(PositionPoint, CenterPoint);
        double alpha = nom/denom;

        std::vector<double> nomd = gammad;
        std::vector<double> denomd = {(1/bg::distance(PositionPoint, CenterPoint))*(Position[0]-Center[0]), (1/bg::distance(PositionPoint, CenterPoint))*(Position[1]-Center[1])};
        std::vector<double> alphad = {(1/denom)*nomd[0]-(nom/pow(denom,2))*denomd[0], (1/denom)*nomd[1]-(nom/pow(denom,2))*denomd[1]};

        std::vector<std::vector<double>> nomdd = gammadd;
        std::vector<double> PositionMinusCenter = {Position[0]-Center[0], Position[1]-Center[1]};
        std::vector<std::vector<double>> PositionPositionOuter = VectorOuterProduct(PositionMinusCenter, PositionMinusCenter);
        std::vector<std::vector<double>> denomdd = {{(1/bg::distance(PositionPoint, CenterPoint)) - (1/pow(bg::distance(PositionPoint, CenterPoint),3))*PositionPositionOuter[0][0], -(1/pow(bg::distance(PositionPoint, CenterPoint),3))*PositionPositionOuter[0][1]}, {-(1/pow(bg::distance(PositionPoint, CenterPoint),3))*PositionPositionOuter[1][0], (1/bg::distance(PositionPoint, CenterPoint)) - (1/pow(bg::distance(PositionPoint, CenterPoint),3))*PositionPositionOuter[1][1]}};
        std::vector<std::vector<double>> NomDenomOuter = VectorOuterProduct(nomd, denomd);
        std::vector<std::vector<double>> DenomNomOuter = VectorOuterProduct(denomd, nomd);
        std::vector<std::vector<double>> DenomDenomOuter = VectorOuterProduct(denomd, denomd);
        std::vector<std::vector<double>> alphadd = {{(1/denom)*nomdd[0][0] - (1/pow(denom,2))*(NomDenomOuter[0][0]+DenomNomOuter[0][0]) + 2*(nom/pow(denom,3))*DenomDenomOuter[0][0] - (nom/pow(denom,2))*denomdd[0][0], (1/denom)*nomdd[0][1] - (1/pow(denom,2))*(NomDenomOuter[0][1]+DenomNomOuter[0][1]) + 2*(nom/pow(denom,3))*DenomDenomOuter[0][1] - (nom/pow(denom,2))*denomdd[0][1]}, {(1/denom)*nomdd[1][0] - (1/pow(denom,2))*(NomDenomOuter[1][0]+DenomNomOuter[1][0]) + 2*(nom/pow(denom,3))*DenomDenomOuter[1][0] - (nom/pow(denom,2))*denomdd[1][0], (1/denom)*nomdd[1][1] - (1/pow(denom,2))*(NomDenomOuter[1][1]+DenomNomOuter[1][1]) + 2*(nom/pow(denom,3))*DenomDenomOuter[1][1] - (nom/pow(denom,2))*denomdd[1][1]}};

        sigma = exp(-mu_2/alpha);
        sigmad = {mu_2*(sigma/pow(alpha,2))*alphad[0], mu_2*(sigma/pow(alpha,2))*alphad[1]};
        std::vector<std::vector<double>> AlphadOuter = VectorOuterProduct(alphad, alphad);
        sigmadd = {{(pow(mu_2,2)*(sigma/pow(alpha,4))-2*mu_2*(sigma/pow(alpha,3)))*AlphadOuter[0][0]+mu_2*(sigma/pow(alpha,2))*alphadd[0][0], (pow(mu_2,2)*(sigma/pow(alpha,4))-2*mu_2*(sigma/pow(alpha,3)))*AlphadOuter[0][1]+mu_2*(sigma/pow(alpha,2))*alphadd[0][1]}, {(pow(mu_2,2)*(sigma/pow(alpha,4))-2*mu_2*(sigma/pow(alpha,3)))*AlphadOuter[1][0]+mu_2*(sigma/pow(alpha,2))*alphadd[1][0], (pow(mu_2,2)*(sigma/pow(alpha,4))-2*mu_2*(sigma/pow(alpha,3)))*AlphadOuter[1][1]+mu_2*(sigma/pow(alpha,2))*alphadd[1][1]}};
    }

    // Construct the output
    OutputStructScalar Output;
    Output.Value = sigma;
    Output.Gradient = sigmad;
    Output.Hessian = sigmadd;

    return Output;
}


OutputStructScalar triangleOutsideImplicit(std::vector<double> Position, TriangleClass Triangle, DiffeoParamsClass DiffeoParams) {
    /**
     * Function that computes beta(x) (i.e., the R-function) for a point x outside a triangle, and its gradient and hessian
     * 
     * Input:
     *  1) Position: Point to consider
     *  2) Triangle: Description of the triangle
     *  3) DiffeoParams: Options for the diffeomorphism construction
     * 
     * Output:
     *  1) Output: Output struct containing the value, gradient and hessian
     */

    // Unwrap parameters
    double p = DiffeoParams.get_p();

    // Make position a point
    point PositionPoint(Position[0], Position[1]);

    // Unwrap triangle properties
    uint16_t depth = Triangle.get_depth();
    std::vector<point> adj_edge = Triangle.get_adj_edge();
    point CenterPoint = Triangle.get_center();
    std::vector<double> Center = {CenterPoint.get<0>(), CenterPoint.get<1>()};

    // Get triangle vertices and normal vectors and convert them
    std::vector<point> TriangleVertexListPoint = Triangle.get_vertices();
    std::vector<point> NormalVectorListPoint = Triangle.get_r_n();
    std::vector<std::vector<double>> TriangleVertexList;
    std::vector<std::vector<double>> NormalVectorList;
    for (size_t i = 0; i < Triangle.get_vertices().size(); i++) {
        TriangleVertexList.push_back({TriangleVertexListPoint[i].get<0>(), TriangleVertexListPoint[i].get<1>()});
        NormalVectorList.push_back({NormalVectorListPoint[i].get<0>(), NormalVectorListPoint[i].get<1>()});
    }

    // Distinguish between the root triangle and all other triangles
    double beta;
    std::vector<double> betad;
    std::vector<std::vector<double>> betadd;
    if ((depth == 0) && adj_edge.empty()) {
        // Find the hyperplane values
        double hyperplane_1 = (Position[0]-TriangleVertexList[2][0])*NormalVectorList[1][0] + (Position[1]-TriangleVertexList[2][1])*NormalVectorList[1][1];
        double hyperplane_2 = (Position[0]-TriangleVertexList[2][0])*NormalVectorList[2][0] + (Position[1]-TriangleVertexList[2][1])*NormalVectorList[2][1];
        double hyperplane_3 = (Position[0]-TriangleVertexList[0][0])*NormalVectorList[0][0] + (Position[1]-TriangleVertexList[0][1])*NormalVectorList[0][1];

        // Compute the R-function
        double hyperplane_12 = hyperplane_1 + hyperplane_2 - pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(1/p));
        double hyperplane_123 = hyperplane_12 + hyperplane_3 - pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(1/p));
        beta = -hyperplane_123;

        // Compute the gradients
        std::vector<double> hyperplane_1d = {NormalVectorList[1][0], NormalVectorList[1][1]};
        std::vector<double> hyperplane_2d = {NormalVectorList[2][0], NormalVectorList[2][1]};
        std::vector<double> hyperplane_3d = {NormalVectorList[0][0], NormalVectorList[0][1]};
        std::vector<double> hyperplane_12d = {(1-((pow(hyperplane_1,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_1d[0] + (1-((pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_2d[0], (1-((pow(hyperplane_1,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_1d[1] + (1-((pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_2d[1]};
        std::vector<double> hyperplane_123d = {(1-((pow(hyperplane_12,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p))))*hyperplane_12d[0] + (1-((pow(hyperplane_3,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p))))*hyperplane_3d[0], (1-((pow(hyperplane_12,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p))))*hyperplane_12d[1] + (1-((pow(hyperplane_3,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p))))*hyperplane_3d[1]};
        betad = {-hyperplane_123d[0], -hyperplane_123d[1]};

        // Compute the hessian
        std::vector<std::vector<double>> hyperplane_1dd = {{0.0, 0.0}, {0.0, 0.0}};
        std::vector<std::vector<double>> hyperplane_2dd = {{0.0, 0.0}, {0.0, 0.0}};
        std::vector<std::vector<double>> hyperplane_3dd = {{0.0, 0.0}, {0.0, 0.0}};
        std::vector<std::vector<double>> Hyperplane1Hyperplane1Outer = VectorOuterProduct(hyperplane_1d, hyperplane_1d);
        std::vector<std::vector<double>> Hyperplane1Hyperplane2Outer = VectorOuterProduct(hyperplane_1d, hyperplane_2d);
        std::vector<std::vector<double>> Hyperplane2Hyperplane1Outer = VectorOuterProduct(hyperplane_2d, hyperplane_1d);
        std::vector<std::vector<double>> Hyperplane2Hyperplane2Outer = VectorOuterProduct(hyperplane_2d, hyperplane_2d);
        std::vector<std::vector<double>> hyperplane_12dd = {{(-(p-1)*((pow(hyperplane_1,p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p)))+(p-1)*((pow(hyperplane_1,2*p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane1Hyperplane1Outer[0][0] + ((p-1)*((pow(hyperplane_1,p-1)*pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane1Hyperplane2Outer[0][0] + (1-((pow(hyperplane_1,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_1dd[0][0] + (-(p-1)*((pow(hyperplane_2,p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p)))+(p-1)*((pow(hyperplane_2,2*p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane2Hyperplane2Outer[0][0] + ((p-1)*((pow(hyperplane_1,p-1)*pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane2Hyperplane1Outer[0][0]  + (1-((pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_2dd[0][0], (-(p-1)*((pow(hyperplane_1,p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p)))+(p-1)*((pow(hyperplane_1,2*p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane1Hyperplane1Outer[0][1] + ((p-1)*((pow(hyperplane_1,p-1)*pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane1Hyperplane2Outer[0][1] + (1-((pow(hyperplane_1,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_1dd[0][1] + (-(p-1)*((pow(hyperplane_2,p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p)))+(p-1)*((pow(hyperplane_2,2*p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane2Hyperplane2Outer[0][1] + ((p-1)*((pow(hyperplane_1,p-1)*pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane2Hyperplane1Outer[0][1]  + (1-((pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_2dd[0][1]}, {(-(p-1)*((pow(hyperplane_1,p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p)))+(p-1)*((pow(hyperplane_1,2*p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane1Hyperplane1Outer[1][0] + ((p-1)*((pow(hyperplane_1,p-1)*pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane1Hyperplane2Outer[1][0] + (1-((pow(hyperplane_1,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_1dd[1][0] + (-(p-1)*((pow(hyperplane_2,p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p)))+(p-1)*((pow(hyperplane_2,2*p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane2Hyperplane2Outer[1][0] + ((p-1)*((pow(hyperplane_1,p-1)*pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane2Hyperplane1Outer[1][0]  + (1-((pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_2dd[1][0], (-(p-1)*((pow(hyperplane_1,p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p)))+(p-1)*((pow(hyperplane_1,2*p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane1Hyperplane1Outer[1][1] + ((p-1)*((pow(hyperplane_1,p-1)*pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane1Hyperplane2Outer[1][1] + (1-((pow(hyperplane_1,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_1dd[1][1] + (-(p-1)*((pow(hyperplane_2,p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p)))+(p-1)*((pow(hyperplane_2,2*p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane2Hyperplane2Outer[1][1] + ((p-1)*((pow(hyperplane_1,p-1)*pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane2Hyperplane1Outer[1][1]  + (1-((pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_2dd[1][1]}};
        std::vector<std::vector<double>> Hyperplane12Hyperplane12Outer = VectorOuterProduct(hyperplane_12d, hyperplane_12d);
        std::vector<std::vector<double>> Hyperplane12Hyperplane3Outer = VectorOuterProduct(hyperplane_12d, hyperplane_3d);
        std::vector<std::vector<double>> Hyperplane3Hyperplane12Outer = VectorOuterProduct(hyperplane_3d, hyperplane_12d);
        std::vector<std::vector<double>> Hyperplane3Hyperplane3Outer = VectorOuterProduct(hyperplane_3d, hyperplane_3d);
        std::vector<std::vector<double>> hyperplane_123dd = {{(-(p-1)*((pow(hyperplane_12,p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p)))+(p-1)*((pow(hyperplane_12,2*p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),1+(p-1)/p))))*Hyperplane12Hyperplane12Outer[0][0] + ((p-1)*((pow(hyperplane_12,p-1)*pow(hyperplane_3,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),1+(p-1)/p))))*Hyperplane12Hyperplane3Outer[0][0] + (1-((pow(hyperplane_12,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p))))*hyperplane_12dd[0][0] + (-(p-1)*((pow(hyperplane_3,p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p)))+(p-1)*((pow(hyperplane_3,2*p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),1+(p-1)/p))))*Hyperplane3Hyperplane3Outer[0][0] + ((p-1)*((pow(hyperplane_12,p-1)*pow(hyperplane_3,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),1+(p-1)/p))))*Hyperplane3Hyperplane12Outer[0][0]  + (1-((pow(hyperplane_3,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p))))*hyperplane_3dd[0][0], (-(p-1)*((pow(hyperplane_12,p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p)))+(p-1)*((pow(hyperplane_12,2*p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),1+(p-1)/p))))*Hyperplane12Hyperplane12Outer[0][1] + ((p-1)*((pow(hyperplane_12,p-1)*pow(hyperplane_3,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),1+(p-1)/p))))*Hyperplane12Hyperplane3Outer[0][1] + (1-((pow(hyperplane_12,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p))))*hyperplane_12dd[0][1] + (-(p-1)*((pow(hyperplane_3,p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p)))+(p-1)*((pow(hyperplane_3,2*p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),1+(p-1)/p))))*Hyperplane3Hyperplane3Outer[0][1] + ((p-1)*((pow(hyperplane_12,p-1)*pow(hyperplane_3,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),1+(p-1)/p))))*Hyperplane3Hyperplane12Outer[0][1]  + (1-((pow(hyperplane_3,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p))))*hyperplane_3dd[0][1]}, {(-(p-1)*((pow(hyperplane_12,p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p)))+(p-1)*((pow(hyperplane_12,2*p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),1+(p-1)/p))))*Hyperplane12Hyperplane12Outer[1][0] + ((p-1)*((pow(hyperplane_12,p-1)*pow(hyperplane_3,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),1+(p-1)/p))))*Hyperplane12Hyperplane3Outer[1][0] + (1-((pow(hyperplane_12,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p))))*hyperplane_12dd[1][0] + (-(p-1)*((pow(hyperplane_3,p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p)))+(p-1)*((pow(hyperplane_3,2*p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),1+(p-1)/p))))*Hyperplane3Hyperplane3Outer[1][0] + ((p-1)*((pow(hyperplane_12,p-1)*pow(hyperplane_3,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),1+(p-1)/p))))*Hyperplane3Hyperplane12Outer[1][0]  + (1-((pow(hyperplane_3,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p))))*hyperplane_3dd[1][0], (-(p-1)*((pow(hyperplane_12,p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p)))+(p-1)*((pow(hyperplane_12,2*p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),1+(p-1)/p))))*Hyperplane12Hyperplane12Outer[1][1] + ((p-1)*((pow(hyperplane_12,p-1)*pow(hyperplane_3,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),1+(p-1)/p))))*Hyperplane12Hyperplane3Outer[1][1] + (1-((pow(hyperplane_12,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p))))*hyperplane_12dd[1][1] + (-(p-1)*((pow(hyperplane_3,p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p)))+(p-1)*((pow(hyperplane_3,2*p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),1+(p-1)/p))))*Hyperplane3Hyperplane3Outer[1][1] + ((p-1)*((pow(hyperplane_12,p-1)*pow(hyperplane_3,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),1+(p-1)/p))))*Hyperplane3Hyperplane12Outer[1][1]  + (1-((pow(hyperplane_3,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_3,p),(p-1)/p))))*hyperplane_3dd[1][1]}};
        betadd = {{-hyperplane_123dd[0][0], -hyperplane_123dd[0][1]}, {-hyperplane_123dd[1][0], -hyperplane_123dd[1][1]}};
    } else {
        // Find the center hyperplane normals
        std::vector<point> RCenterNormalVectorListPoint = Triangle.get_r_center_n();
        std::vector<std::vector<double>> RCenterNormalVectorList;
        for (size_t i = 0; i < RCenterNormalVectorListPoint.size(); i++) {
            RCenterNormalVectorList.push_back({RCenterNormalVectorListPoint[i].get<0>(), RCenterNormalVectorListPoint[i].get<1>()});
        }

        // Find the hyperplane values
        double hyperplane_1 = (Position[0]-TriangleVertexList[2][0])*NormalVectorList[1][0] + (Position[1]-TriangleVertexList[2][1])*NormalVectorList[1][1];
        double hyperplane_2 = (Position[0]-TriangleVertexList[2][0])*NormalVectorList[2][0] + (Position[1]-TriangleVertexList[2][1])*NormalVectorList[2][1];
        double hyperplane_3 = (Position[0]-Center[0])*RCenterNormalVectorList[0][0] + (Position[1]-Center[1])*RCenterNormalVectorList[0][1];
        double hyperplane_4 = (Position[0]-Center[0])*RCenterNormalVectorList[1][0] + (Position[1]-Center[1])*RCenterNormalVectorList[1][1];

        // Compute the R-function
        double hyperplane_12 = hyperplane_1 + hyperplane_2 - pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(1/p));
        double hyperplane_34 = hyperplane_3 + hyperplane_4 - pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(1/p));
        double hyperplane_1234 = hyperplane_12 + hyperplane_34 - pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(1/p));
        beta = -hyperplane_1234;

        // Compute the gradients
        std::vector<double> hyperplane_1d = NormalVectorList[1];
        std::vector<double> hyperplane_2d = NormalVectorList[2];
        std::vector<double> hyperplane_3d = RCenterNormalVectorList[0];
        std::vector<double> hyperplane_4d = RCenterNormalVectorList[1];
        std::vector<double> hyperplane_12d = {(1-((pow(hyperplane_1,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_1d[0] + (1-((pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_2d[0], (1-((pow(hyperplane_1,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_1d[1] + (1-((pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_2d[1]};
        std::vector<double> hyperplane_34d = {(1-((pow(hyperplane_3,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p))))*hyperplane_3d[0] + (1-((pow(hyperplane_4,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p))))*hyperplane_4d[0], (1-((pow(hyperplane_3,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p))))*hyperplane_3d[1] + (1-((pow(hyperplane_4,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p))))*hyperplane_4d[1]};
        std::vector<double> hyperplane_1234d = {(1-((pow(hyperplane_12,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p))))*hyperplane_12d[0] + (1-((pow(hyperplane_34,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p))))*hyperplane_34d[0], (1-((pow(hyperplane_12,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p))))*hyperplane_12d[1] + (1-((pow(hyperplane_34,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p))))*hyperplane_34d[1]};
        betad = {-hyperplane_1234d[0], -hyperplane_1234d[1]};

        // Compute the hessian
        std::vector<std::vector<double>> hyperplane_1dd = {{0.0, 0.0}, {0.0, 0.0}};
        std::vector<std::vector<double>> hyperplane_2dd = {{0.0, 0.0}, {0.0, 0.0}};
        std::vector<std::vector<double>> hyperplane_3dd = {{0.0, 0.0}, {0.0, 0.0}};
        std::vector<std::vector<double>> hyperplane_4dd = {{0.0, 0.0}, {0.0, 0.0}};
        std::vector<std::vector<double>> Hyperplane1Hyperplane1Outer = VectorOuterProduct(hyperplane_1d, hyperplane_1d);
        std::vector<std::vector<double>> Hyperplane1Hyperplane2Outer = VectorOuterProduct(hyperplane_1d, hyperplane_2d);
        std::vector<std::vector<double>> Hyperplane2Hyperplane1Outer = VectorOuterProduct(hyperplane_2d, hyperplane_1d);
        std::vector<std::vector<double>> Hyperplane2Hyperplane2Outer = VectorOuterProduct(hyperplane_2d, hyperplane_2d);
        std::vector<std::vector<double>> hyperplane_12dd = {{(-(p-1)*((pow(hyperplane_1,p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p)))+(p-1)*((pow(hyperplane_1,2*p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane1Hyperplane1Outer[0][0] + ((p-1)*((pow(hyperplane_1,p-1)*pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane1Hyperplane2Outer[0][0] + (1-((pow(hyperplane_1,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_1dd[0][0] + (-(p-1)*((pow(hyperplane_2,p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p)))+(p-1)*((pow(hyperplane_2,2*p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane2Hyperplane2Outer[0][0] + ((p-1)*((pow(hyperplane_1,p-1)*pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane2Hyperplane1Outer[0][0]  + (1-((pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_2dd[0][0], (-(p-1)*((pow(hyperplane_1,p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p)))+(p-1)*((pow(hyperplane_1,2*p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane1Hyperplane1Outer[0][1] + ((p-1)*((pow(hyperplane_1,p-1)*pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane1Hyperplane2Outer[0][1] + (1-((pow(hyperplane_1,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_1dd[0][1] + (-(p-1)*((pow(hyperplane_2,p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p)))+(p-1)*((pow(hyperplane_2,2*p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane2Hyperplane2Outer[0][1] + ((p-1)*((pow(hyperplane_1,p-1)*pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane2Hyperplane1Outer[0][1]  + (1-((pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_2dd[0][1]}, {(-(p-1)*((pow(hyperplane_1,p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p)))+(p-1)*((pow(hyperplane_1,2*p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane1Hyperplane1Outer[1][0] + ((p-1)*((pow(hyperplane_1,p-1)*pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane1Hyperplane2Outer[1][0] + (1-((pow(hyperplane_1,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_1dd[1][0] + (-(p-1)*((pow(hyperplane_2,p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p)))+(p-1)*((pow(hyperplane_2,2*p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane2Hyperplane2Outer[1][0] + ((p-1)*((pow(hyperplane_1,p-1)*pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane2Hyperplane1Outer[1][0]  + (1-((pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_2dd[1][0], (-(p-1)*((pow(hyperplane_1,p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p)))+(p-1)*((pow(hyperplane_1,2*p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane1Hyperplane1Outer[1][1] + ((p-1)*((pow(hyperplane_1,p-1)*pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane1Hyperplane2Outer[1][1] + (1-((pow(hyperplane_1,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_1dd[1][1] + (-(p-1)*((pow(hyperplane_2,p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p)))+(p-1)*((pow(hyperplane_2,2*p-2))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane2Hyperplane2Outer[1][1] + ((p-1)*((pow(hyperplane_1,p-1)*pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),1+(p-1)/p))))*Hyperplane2Hyperplane1Outer[1][1]  + (1-((pow(hyperplane_2,p-1))/(pow(pow(hyperplane_1,p)+pow(hyperplane_2,p),(p-1)/p))))*hyperplane_2dd[1][1]}};

        std::vector<std::vector<double>> Hyperplane3Hyperplane3Outer = VectorOuterProduct(hyperplane_3d, hyperplane_3d);
        std::vector<std::vector<double>> Hyperplane3Hyperplane4Outer = VectorOuterProduct(hyperplane_3d, hyperplane_4d);
        std::vector<std::vector<double>> Hyperplane4Hyperplane3Outer = VectorOuterProduct(hyperplane_4d, hyperplane_3d);
        std::vector<std::vector<double>> Hyperplane4Hyperplane4Outer = VectorOuterProduct(hyperplane_4d, hyperplane_4d);
        std::vector<std::vector<double>> hyperplane_34dd = {{(-(p-1)*((pow(hyperplane_3,p-2))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p)))+(p-1)*((pow(hyperplane_3,2*p-2))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),1+(p-1)/p))))*Hyperplane3Hyperplane3Outer[0][0] + ((p-1)*((pow(hyperplane_3,p-1)*pow(hyperplane_4,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),1+(p-1)/p))))*Hyperplane3Hyperplane4Outer[0][0] + (1-((pow(hyperplane_3,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p))))*hyperplane_3dd[0][0] + (-(p-1)*((pow(hyperplane_4,p-2))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p)))+(p-1)*((pow(hyperplane_4,2*p-2))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),1+(p-1)/p))))*Hyperplane4Hyperplane4Outer[0][0] + ((p-1)*((pow(hyperplane_3,p-1)*pow(hyperplane_4,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),1+(p-1)/p))))*Hyperplane4Hyperplane3Outer[0][0]  + (1-((pow(hyperplane_4,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p))))*hyperplane_4dd[0][0], (-(p-1)*((pow(hyperplane_3,p-2))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p)))+(p-1)*((pow(hyperplane_3,2*p-2))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),1+(p-1)/p))))*Hyperplane3Hyperplane3Outer[0][1] + ((p-1)*((pow(hyperplane_3,p-1)*pow(hyperplane_4,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),1+(p-1)/p))))*Hyperplane3Hyperplane4Outer[0][1] + (1-((pow(hyperplane_3,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p))))*hyperplane_3dd[0][1] + (-(p-1)*((pow(hyperplane_4,p-2))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p)))+(p-1)*((pow(hyperplane_4,2*p-2))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),1+(p-1)/p))))*Hyperplane4Hyperplane4Outer[0][1] + ((p-1)*((pow(hyperplane_3,p-1)*pow(hyperplane_4,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),1+(p-1)/p))))*Hyperplane4Hyperplane3Outer[0][1]  + (1-((pow(hyperplane_4,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p))))*hyperplane_4dd[0][1]}, {(-(p-1)*((pow(hyperplane_3,p-2))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p)))+(p-1)*((pow(hyperplane_3,2*p-2))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),1+(p-1)/p))))*Hyperplane3Hyperplane3Outer[1][0] + ((p-1)*((pow(hyperplane_3,p-1)*pow(hyperplane_4,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),1+(p-1)/p))))*Hyperplane3Hyperplane4Outer[1][0] + (1-((pow(hyperplane_3,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p))))*hyperplane_3dd[1][0] + (-(p-1)*((pow(hyperplane_4,p-2))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p)))+(p-1)*((pow(hyperplane_4,2*p-2))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),1+(p-1)/p))))*Hyperplane4Hyperplane4Outer[1][0] + ((p-1)*((pow(hyperplane_3,p-1)*pow(hyperplane_4,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),1+(p-1)/p))))*Hyperplane4Hyperplane3Outer[1][0]  + (1-((pow(hyperplane_4,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p))))*hyperplane_4dd[1][0], (-(p-1)*((pow(hyperplane_3,p-2))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p)))+(p-1)*((pow(hyperplane_3,2*p-2))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),1+(p-1)/p))))*Hyperplane3Hyperplane3Outer[1][1] + ((p-1)*((pow(hyperplane_3,p-1)*pow(hyperplane_4,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),1+(p-1)/p))))*Hyperplane3Hyperplane4Outer[1][1] + (1-((pow(hyperplane_3,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p))))*hyperplane_3dd[1][1] + (-(p-1)*((pow(hyperplane_4,p-2))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p)))+(p-1)*((pow(hyperplane_4,2*p-2))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),1+(p-1)/p))))*Hyperplane4Hyperplane4Outer[1][1] + ((p-1)*((pow(hyperplane_3,p-1)*pow(hyperplane_4,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),1+(p-1)/p))))*Hyperplane4Hyperplane3Outer[1][1]  + (1-((pow(hyperplane_4,p-1))/(pow(pow(hyperplane_3,p)+pow(hyperplane_4,p),(p-1)/p))))*hyperplane_4dd[1][1]}};

        std::vector<std::vector<double>> Hyperplane12Hyperplane12Outer = VectorOuterProduct(hyperplane_12d, hyperplane_12d);
        std::vector<std::vector<double>> Hyperplane12Hyperplane34Outer = VectorOuterProduct(hyperplane_12d, hyperplane_34d);
        std::vector<std::vector<double>> Hyperplane34Hyperplane12Outer = VectorOuterProduct(hyperplane_34d, hyperplane_12d);
        std::vector<std::vector<double>> Hyperplane34Hyperplane34Outer = VectorOuterProduct(hyperplane_34d, hyperplane_34d);
        std::vector<std::vector<double>> hyperplane_1234dd = {{(-(p-1)*((pow(hyperplane_12,p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p)))+(p-1)*((pow(hyperplane_12,2*p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),1+(p-1)/p))))*Hyperplane12Hyperplane12Outer[0][0] + ((p-1)*((pow(hyperplane_12,p-1)*pow(hyperplane_34,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),1+(p-1)/p))))*Hyperplane12Hyperplane34Outer[0][0] + (1-((pow(hyperplane_12,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p))))*hyperplane_12dd[0][0] + (-(p-1)*((pow(hyperplane_34,p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p)))+(p-1)*((pow(hyperplane_34,2*p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),1+(p-1)/p))))*Hyperplane34Hyperplane34Outer[0][0] + ((p-1)*((pow(hyperplane_12,p-1)*pow(hyperplane_34,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),1+(p-1)/p))))*Hyperplane34Hyperplane12Outer[0][0]  + (1-((pow(hyperplane_34,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p))))*hyperplane_34dd[0][0], (-(p-1)*((pow(hyperplane_12,p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p)))+(p-1)*((pow(hyperplane_12,2*p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),1+(p-1)/p))))*Hyperplane12Hyperplane12Outer[0][1] + ((p-1)*((pow(hyperplane_12,p-1)*pow(hyperplane_34,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),1+(p-1)/p))))*Hyperplane12Hyperplane34Outer[0][1] + (1-((pow(hyperplane_12,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p))))*hyperplane_12dd[0][1] + (-(p-1)*((pow(hyperplane_34,p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p)))+(p-1)*((pow(hyperplane_34,2*p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),1+(p-1)/p))))*Hyperplane34Hyperplane34Outer[0][1] + ((p-1)*((pow(hyperplane_12,p-1)*pow(hyperplane_34,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),1+(p-1)/p))))*Hyperplane34Hyperplane12Outer[0][1]  + (1-((pow(hyperplane_34,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p))))*hyperplane_34dd[0][1]}, {(-(p-1)*((pow(hyperplane_12,p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p)))+(p-1)*((pow(hyperplane_12,2*p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),1+(p-1)/p))))*Hyperplane12Hyperplane12Outer[1][0] + ((p-1)*((pow(hyperplane_12,p-1)*pow(hyperplane_34,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),1+(p-1)/p))))*Hyperplane12Hyperplane34Outer[1][0] + (1-((pow(hyperplane_12,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p))))*hyperplane_12dd[1][0] + (-(p-1)*((pow(hyperplane_34,p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p)))+(p-1)*((pow(hyperplane_34,2*p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),1+(p-1)/p))))*Hyperplane34Hyperplane34Outer[1][0] + ((p-1)*((pow(hyperplane_12,p-1)*pow(hyperplane_34,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),1+(p-1)/p))))*Hyperplane34Hyperplane12Outer[1][0]  + (1-((pow(hyperplane_34,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p))))*hyperplane_34dd[1][0], (-(p-1)*((pow(hyperplane_12,p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p)))+(p-1)*((pow(hyperplane_12,2*p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),1+(p-1)/p))))*Hyperplane12Hyperplane12Outer[1][1] + ((p-1)*((pow(hyperplane_12,p-1)*pow(hyperplane_34,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),1+(p-1)/p))))*Hyperplane12Hyperplane34Outer[1][1] + (1-((pow(hyperplane_12,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p))))*hyperplane_12dd[1][1] + (-(p-1)*((pow(hyperplane_34,p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p)))+(p-1)*((pow(hyperplane_34,2*p-2))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),1+(p-1)/p))))*Hyperplane34Hyperplane34Outer[1][1] + ((p-1)*((pow(hyperplane_12,p-1)*pow(hyperplane_34,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),1+(p-1)/p))))*Hyperplane34Hyperplane12Outer[1][1]  + (1-((pow(hyperplane_34,p-1))/(pow(pow(hyperplane_12,p)+pow(hyperplane_34,p),(p-1)/p))))*hyperplane_34dd[1][1]}};
        betadd = {{-hyperplane_1234dd[0][0], -hyperplane_1234dd[0][1]}, {-hyperplane_1234dd[1][0], -hyperplane_1234dd[1][1]}};
    }

    // Construct the output
    OutputStructScalar Output;
    Output.Value = beta;
    Output.Gradient = betad;
    Output.Hessian = betadd;

    return Output;
}


OutputStructScalar polygonOutsideImplicit(std::vector<double> Position, PolygonClass PolygonUsed, DiffeoParamsClass DiffeoParams) {
    /**
     * Function that computes beta(x) (i.e., the R-function) for a point x outside a polygon, and its gradient and hessian
     * 
     * Input:
     *  1) Position: Point to consider
     *  2) PolygonUsed: Description of the polygon
     *  3) DiffeoParams: Options for the diffeomorphism construction
     * 
     * Output:
     *  1) Output: Output struct containing the value, gradient and hessian
     */

    // Unwrap parameters
    double p = DiffeoParams.get_p();

    // Make position a point
    point PositionPoint(Position[0], Position[1]);

    // Get polygon vertices and normal vectors and convert them
    std::vector<point> PolygonVertexListPoint = PolygonUsed.get_augmented_vertices();
    std::vector<point> NormalVectorListPoint = PolygonUsed.get_r_n();
    std::vector<std::vector<double>> PolygonVertexList;
    std::vector<std::vector<double>> NormalVectorList;
    for (size_t i = 0; i < PolygonVertexListPoint.size(); i++) {
        PolygonVertexList.push_back({PolygonVertexListPoint[i].get<0>(), PolygonVertexListPoint[i].get<1>()});
        NormalVectorList.push_back({NormalVectorListPoint[i].get<0>(), NormalVectorListPoint[i].get<1>()});
    }

    // Compute hyperplane functions
    std::vector<double> hyperplane(PolygonVertexList.size(), 0.0);
    for (size_t i = 0; i < PolygonVertexList.size(); i++) {
        hyperplane[i] = (Position[0]-PolygonVertexList[i][0])*NormalVectorList[i][0] + (Position[1]-PolygonVertexList[i][1])*NormalVectorList[i][1];
    }

    // Compute the R-function and its gradient and hessian
    std::vector<std::vector<double>> N0N0Outer = VectorOuterProduct(NormalVectorList[0], NormalVectorList[0]);
    std::vector<std::vector<double>> N0N1Outer = VectorOuterProduct(NormalVectorList[0], NormalVectorList[1]);
    std::vector<std::vector<double>> N1N0Outer = VectorOuterProduct(NormalVectorList[1], NormalVectorList[0]);
    std::vector<std::vector<double>> N1N1Outer = VectorOuterProduct(NormalVectorList[1], NormalVectorList[1]);
    std::vector<std::vector<double>> betadd = {{(-(p-1)*((pow(hyperplane[0],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[0],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N0Outer[0][0] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N1Outer[0][0] + (1-((pow(hyperplane[0],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0 + (-(p-1)*((pow(hyperplane[1],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[1],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N1Outer[0][0] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N0Outer[0][0] + (1-((pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0, (-(p-1)*((pow(hyperplane[0],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[0],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N0Outer[0][1] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N1Outer[0][1] + (1-((pow(hyperplane[0],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0 + (-(p-1)*((pow(hyperplane[1],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[1],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N1Outer[0][1] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N0Outer[0][1] + (1-((pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0}, {(-(p-1)*((pow(hyperplane[0],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[0],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N0Outer[1][0] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N1Outer[1][0] + (1-((pow(hyperplane[0],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0 + (-(p-1)*((pow(hyperplane[1],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[1],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N1Outer[1][0] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N0Outer[1][0] + (1-((pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0, (-(p-1)*((pow(hyperplane[0],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[0],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N0Outer[1][1] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N1Outer[1][1] + (1-((pow(hyperplane[0],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0 + (-(p-1)*((pow(hyperplane[1],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[1],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N1Outer[1][1] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N0Outer[1][1] + (1-((pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0}};

    std::vector<double> betad = {(1-((pow(hyperplane[0],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*NormalVectorList[0][0] + (1-((pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*NormalVectorList[1][0], (1-((pow(hyperplane[0],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*NormalVectorList[0][1] + (1-((pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*NormalVectorList[1][1]};

    double beta = hyperplane[0] + hyperplane[1] - pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1/p);
    
    for (size_t i = 2; i < hyperplane.size(); i++) {
        std::vector<std::vector<double>> BetadBetadOuter = VectorOuterProduct(betad, betad);
        std::vector<std::vector<double>> BetaNiOuter = VectorOuterProduct(betad, NormalVectorList[i]);
        std::vector<std::vector<double>> NiBetadOuter = VectorOuterProduct(NormalVectorList[i], betad);
        std::vector<std::vector<double>> NiNiOuter = VectorOuterProduct(NormalVectorList[i], NormalVectorList[i]);
        betadd = {{(-(p-1)*((pow(beta,p-2))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(beta,2*p-2))/(pow(pow(beta,p)+pow(hyperplane[i],p),1+(p-1)/p))))*BetadBetadOuter[0][0] + ((p-1)*((pow(beta,p-1)*pow(hyperplane[i],p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),1+(p-1)/p))))*BetaNiOuter[0][0] + (1-((pow(beta,p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p))))*betadd[0][0] + (-(p-1)*((pow(hyperplane[i],p-2))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(hyperplane[i],2*p-2))/(pow(pow(beta,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiNiOuter[0][0] + ((p-1)*((pow(beta,p-1)*pow(hyperplane[i],p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiBetadOuter[0][0] + (1-((pow(hyperplane[i],p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p))))*0.0, (-(p-1)*((pow(beta,p-2))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(beta,2*p-2))/(pow(pow(beta,p)+pow(hyperplane[i],p),1+(p-1)/p))))*BetadBetadOuter[0][1] + ((p-1)*((pow(beta,p-1)*pow(hyperplane[i],p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),1+(p-1)/p))))*BetaNiOuter[0][1] + (1-((pow(beta,p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p))))*betadd[0][1] + (-(p-1)*((pow(hyperplane[i],p-2))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(hyperplane[i],2*p-2))/(pow(pow(beta,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiNiOuter[0][1] + ((p-1)*((pow(beta,p-1)*pow(hyperplane[i],p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiBetadOuter[0][1] + (1-((pow(hyperplane[i],p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p))))*0.0}, {(-(p-1)*((pow(beta,p-2))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(beta,2*p-2))/(pow(pow(beta,p)+pow(hyperplane[i],p),1+(p-1)/p))))*BetadBetadOuter[1][0] + ((p-1)*((pow(beta,p-1)*pow(hyperplane[i],p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),1+(p-1)/p))))*BetaNiOuter[1][0] + (1-((pow(beta,p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p))))*betadd[1][0] + (-(p-1)*((pow(hyperplane[i],p-2))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(hyperplane[i],2*p-2))/(pow(pow(beta,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiNiOuter[1][0] + ((p-1)*((pow(beta,p-1)*pow(hyperplane[i],p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiBetadOuter[1][0] + (1-((pow(hyperplane[i],p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p))))*0.0, (-(p-1)*((pow(beta,p-2))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(beta,2*p-2))/(pow(pow(beta,p)+pow(hyperplane[i],p),1+(p-1)/p))))*BetadBetadOuter[1][1] + ((p-1)*((pow(beta,p-1)*pow(hyperplane[i],p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),1+(p-1)/p))))*BetaNiOuter[1][1] + (1-((pow(beta,p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p))))*betadd[1][1] + (-(p-1)*((pow(hyperplane[i],p-2))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(hyperplane[i],2*p-2))/(pow(pow(beta,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiNiOuter[1][1] + ((p-1)*((pow(beta,p-1)*pow(hyperplane[i],p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiBetadOuter[1][1] + (1-((pow(hyperplane[i],p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p))))*0.0}};

        betad = {(1-((pow(beta,p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p))))*betad[0] + (1-((pow(hyperplane[i],p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p))))*NormalVectorList[i][0], (1-((pow(beta,p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p))))*betad[1] + (1-((pow(hyperplane[i],p-1))/(pow(pow(beta,p)+pow(hyperplane[i],p),(p-1)/p))))*NormalVectorList[i][1]};

        beta = beta + hyperplane[i] - pow(pow(beta,p)+pow(hyperplane[i],p),1/p);
    }

    // Construct the output
    OutputStructScalar Output;
    Output.Value = -beta;
    Output.Gradient = {-betad[0], -betad[1]};;
    Output.Hessian = {{-betadd[0][0], -betadd[0][1]}, {-betadd[1][0], -betadd[1][1]}};;

    return Output;
}


OutputStructScalar triangleInsideImplicit(std::vector<double> Position, TriangleClass Triangle, DiffeoParamsClass DiffeoParams) {
    /**
     * Function that computes gamma(x) (i.e., the R-function) for a point x inside an enclosing polygon, and its gradient and hessian
     * 
     * Input:
     *  1) Position: Point to consider
     *  2) Triangle: Description of the triangle
     *  3) DiffeoParams: Options for the diffeomorphism construction
     * 
     * Output:
     *  1) Output: Output struct containing the value, gradient and hessian
     */

    // Unwrap parameters
    double p = DiffeoParams.get_p();

    // Make position a point
    point PositionPoint(Position[0], Position[1]);

    // Get collar vertices and normal vectors and convert them
    std::vector<point> TriangleVertexTildeListPoint = Triangle.get_vertices_tilde();
    std::vector<point> NormalVectorTildeListPoint = Triangle.get_r_tilde_n();
    std::vector<std::vector<double>> TriangleVertexTildeList;
    std::vector<std::vector<double>> NormalVectorTildeList;
    for (size_t i = 0; i < TriangleVertexTildeListPoint.size(); i++) {
        TriangleVertexTildeList.push_back({TriangleVertexTildeListPoint[i].get<0>(), TriangleVertexTildeListPoint[i].get<1>()});
        NormalVectorTildeList.push_back({NormalVectorTildeListPoint[i].get<0>(), NormalVectorTildeListPoint[i].get<1>()});
    }

    // Compute hyperplane functions
    std::vector<double> hyperplane(TriangleVertexTildeList.size(), 0.0);
    for (size_t i = 0; i < TriangleVertexTildeList.size(); i++) {
        hyperplane[i] = (Position[0]-TriangleVertexTildeList[i][0])*NormalVectorTildeList[i][0] + (Position[1]-TriangleVertexTildeList[i][1])*NormalVectorTildeList[i][1];
    }

    // Compute the R-function and its gradient and hessian
    std::vector<std::vector<double>> N0N0Outer = VectorOuterProduct(NormalVectorTildeList[0], NormalVectorTildeList[0]);
    std::vector<std::vector<double>> N0N1Outer = VectorOuterProduct(NormalVectorTildeList[0], NormalVectorTildeList[1]);
    std::vector<std::vector<double>> N1N0Outer = VectorOuterProduct(NormalVectorTildeList[1], NormalVectorTildeList[0]);
    std::vector<std::vector<double>> N1N1Outer = VectorOuterProduct(NormalVectorTildeList[1], NormalVectorTildeList[1]);
    std::vector<std::vector<double>> gammadd = {{(-(p-1)*((pow(hyperplane[0],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[0],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N0Outer[0][0] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N1Outer[0][0] + (1-((pow(hyperplane[0],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0 + (-(p-1)*((pow(hyperplane[1],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[1],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N1Outer[0][0] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N0Outer[0][0] + (1-((pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0, (-(p-1)*((pow(hyperplane[0],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[0],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N0Outer[0][1] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N1Outer[0][1] + (1-((pow(hyperplane[0],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0 + (-(p-1)*((pow(hyperplane[1],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[1],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N1Outer[0][1] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N0Outer[0][1] + (1-((pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0}, {(-(p-1)*((pow(hyperplane[0],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[0],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N0Outer[1][0] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N1Outer[1][0] + (1-((pow(hyperplane[0],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0 + (-(p-1)*((pow(hyperplane[1],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[1],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N1Outer[1][0] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N0Outer[1][0] + (1-((pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0, (-(p-1)*((pow(hyperplane[0],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[0],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N0Outer[1][1] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N1Outer[1][1] + (1-((pow(hyperplane[0],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0 + (-(p-1)*((pow(hyperplane[1],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[1],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N1Outer[1][1] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N0Outer[1][1] + (1-((pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0}};

    std::vector<double> gammad = {(1-((pow(hyperplane[0],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*NormalVectorTildeList[0][0] + (1-((pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*NormalVectorTildeList[1][0], (1-((pow(hyperplane[0],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*NormalVectorTildeList[0][1] + (1-((pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*NormalVectorTildeList[1][1]};

    double gamma = hyperplane[0] + hyperplane[1] - pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1/p);
    
    for (size_t i = 2; i < hyperplane.size(); i++) {
        std::vector<std::vector<double>> GammadGammadOuter = VectorOuterProduct(gammad, gammad);
        std::vector<std::vector<double>> GammadNiOuter = VectorOuterProduct(gammad, NormalVectorTildeList[i]);
        std::vector<std::vector<double>> NiGammadOuter = VectorOuterProduct(NormalVectorTildeList[i], gammad);
        std::vector<std::vector<double>> NiNiOuter = VectorOuterProduct(NormalVectorTildeList[i], NormalVectorTildeList[i]);
        gammadd = {{(-(p-1)*((pow(gamma,p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(gamma,2*p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*GammadGammadOuter[0][0] + ((p-1)*((pow(gamma,p-1)*pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*GammadNiOuter[0][0] + (1-((pow(gamma,p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*gammadd[0][0] + (-(p-1)*((pow(hyperplane[i],p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(hyperplane[i],2*p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiNiOuter[0][0] + ((p-1)*((pow(gamma,p-1)*pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiGammadOuter[0][0] + (1-((pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*0.0, (-(p-1)*((pow(gamma,p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(gamma,2*p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*GammadGammadOuter[0][1] + ((p-1)*((pow(gamma,p-1)*pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*GammadNiOuter[0][1] + (1-((pow(gamma,p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*gammadd[0][1] + (-(p-1)*((pow(hyperplane[i],p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(hyperplane[i],2*p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiNiOuter[0][1] + ((p-1)*((pow(gamma,p-1)*pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiGammadOuter[0][1] + (1-((pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*0.0}, {(-(p-1)*((pow(gamma,p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(gamma,2*p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*GammadGammadOuter[1][0] + ((p-1)*((pow(gamma,p-1)*pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*GammadNiOuter[1][0] + (1-((pow(gamma,p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*gammadd[1][0] + (-(p-1)*((pow(hyperplane[i],p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(hyperplane[i],2*p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiNiOuter[1][0] + ((p-1)*((pow(gamma,p-1)*pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiGammadOuter[1][0] + (1-((pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*0.0, (-(p-1)*((pow(gamma,p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(gamma,2*p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*GammadGammadOuter[1][1] + ((p-1)*((pow(gamma,p-1)*pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*GammadNiOuter[1][1] + (1-((pow(gamma,p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*gammadd[1][1] + (-(p-1)*((pow(hyperplane[i],p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(hyperplane[i],2*p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiNiOuter[1][1] + ((p-1)*((pow(gamma,p-1)*pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiGammadOuter[1][1] + (1-((pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*0.0}};

        gammad = {(1-((pow(gamma,p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*gammad[0] + (1-((pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*NormalVectorTildeList[i][0], (1-((pow(gamma,p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*gammad[1] + (1-((pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*NormalVectorTildeList[i][1]};

        gamma = gamma + hyperplane[i] - pow(pow(gamma,p)+pow(hyperplane[i],p),1/p);
    }

    // Construct the output
    OutputStructScalar Output;
    Output.Value = gamma;
    Output.Gradient = gammad;
    Output.Hessian = gammadd;

    return Output;
}


OutputStructScalar polygonInsideImplicit(std::vector<double> Position, PolygonClass PolygonUsed, DiffeoParamsClass DiffeoParams) {
    /**
     * Function that computes gamma(x) (i.e., the R-function) for a point x inside an enclosing polygon, and its gradient and hessian
     * 
     * Input:
     *  1) Position: Point to consider
     *  2) PolygonUsed: Description of the polygon
     *  3) DiffeoParams: Options for the diffeomorphism construction
     * 
     * Output:
     *  1) Output: Output struct containing the value, gradient and hessian
     */

    // Unwrap parameters
    double p = DiffeoParams.get_p();

    // Make position a point
    point PositionPoint(Position[0], Position[1]);

    // Get collar vertices and normal vectors and convert them
    std::vector<point> PolygonVertexTildeListPoint = PolygonUsed.get_vertices_tilde();
    std::vector<point> NormalVectorTildeListPoint = PolygonUsed.get_r_tilde_n();
    std::vector<std::vector<double>> PolygonVertexTildeList;
    std::vector<std::vector<double>> NormalVectorTildeList;
    for (size_t i = 0; i < PolygonVertexTildeListPoint.size(); i++) {
        PolygonVertexTildeList.push_back({PolygonVertexTildeListPoint[i].get<0>(), PolygonVertexTildeListPoint[i].get<1>()});
        NormalVectorTildeList.push_back({NormalVectorTildeListPoint[i].get<0>(), NormalVectorTildeListPoint[i].get<1>()});
    }

    // Compute hyperplane functions
    std::vector<double> hyperplane(PolygonVertexTildeList.size(), 0.0);
    for (size_t i = 0; i < PolygonVertexTildeList.size(); i++) {
        hyperplane[i] = (Position[0]-PolygonVertexTildeList[i][0])*NormalVectorTildeList[i][0] + (Position[1]-PolygonVertexTildeList[i][1])*NormalVectorTildeList[i][1];
    }

    // Compute the R-function and its gradient and hessian
    std::vector<std::vector<double>> N0N0Outer = VectorOuterProduct(NormalVectorTildeList[0], NormalVectorTildeList[0]);
    std::vector<std::vector<double>> N0N1Outer = VectorOuterProduct(NormalVectorTildeList[0], NormalVectorTildeList[1]);
    std::vector<std::vector<double>> N1N0Outer = VectorOuterProduct(NormalVectorTildeList[1], NormalVectorTildeList[0]);
    std::vector<std::vector<double>> N1N1Outer = VectorOuterProduct(NormalVectorTildeList[1], NormalVectorTildeList[1]);
    std::vector<std::vector<double>> gammadd = {{(-(p-1)*((pow(hyperplane[0],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[0],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N0Outer[0][0] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N1Outer[0][0] + (1-((pow(hyperplane[0],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0 + (-(p-1)*((pow(hyperplane[1],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[1],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N1Outer[0][0] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N0Outer[0][0] + (1-((pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0, (-(p-1)*((pow(hyperplane[0],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[0],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N0Outer[0][1] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N1Outer[0][1] + (1-((pow(hyperplane[0],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0 + (-(p-1)*((pow(hyperplane[1],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[1],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N1Outer[0][1] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N0Outer[0][1] + (1-((pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0}, {(-(p-1)*((pow(hyperplane[0],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[0],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N0Outer[1][0] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N1Outer[1][0] + (1-((pow(hyperplane[0],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0 + (-(p-1)*((pow(hyperplane[1],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[1],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N1Outer[1][0] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N0Outer[1][0] + (1-((pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0, (-(p-1)*((pow(hyperplane[0],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[0],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N0Outer[1][1] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N0N1Outer[1][1] + (1-((pow(hyperplane[0],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0 + (-(p-1)*((pow(hyperplane[1],p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p)))+(p-1)*((pow(hyperplane[1],2*p-2))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N1Outer[1][1] + ((p-1)*((pow(hyperplane[0],p-1)*pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1+(p-1)/p))))*N1N0Outer[1][1] + (1-((pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*0.0}};

    std::vector<double> gammad = {(1-((pow(hyperplane[0],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*NormalVectorTildeList[0][0] + (1-((pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*NormalVectorTildeList[1][0], (1-((pow(hyperplane[0],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*NormalVectorTildeList[0][1] + (1-((pow(hyperplane[1],p-1))/(pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),(p-1)/p))))*NormalVectorTildeList[1][1]};

    double gamma = hyperplane[0] + hyperplane[1] - pow(pow(hyperplane[0],p)+pow(hyperplane[1],p),1/p);
    
    for (size_t i = 2; i < hyperplane.size(); i++) {
        std::vector<std::vector<double>> GammadGammadOuter = VectorOuterProduct(gammad, gammad);
        std::vector<std::vector<double>> GammadNiOuter = VectorOuterProduct(gammad, NormalVectorTildeList[i]);
        std::vector<std::vector<double>> NiGammadOuter = VectorOuterProduct(NormalVectorTildeList[i], gammad);
        std::vector<std::vector<double>> NiNiOuter = VectorOuterProduct(NormalVectorTildeList[i], NormalVectorTildeList[i]);
        gammadd = {{(-(p-1)*((pow(gamma,p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(gamma,2*p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*GammadGammadOuter[0][0] + ((p-1)*((pow(gamma,p-1)*pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*GammadNiOuter[0][0] + (1-((pow(gamma,p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*gammadd[0][0] + (-(p-1)*((pow(hyperplane[i],p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(hyperplane[i],2*p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiNiOuter[0][0] + ((p-1)*((pow(gamma,p-1)*pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiGammadOuter[0][0] + (1-((pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*0.0, (-(p-1)*((pow(gamma,p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(gamma,2*p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*GammadGammadOuter[0][1] + ((p-1)*((pow(gamma,p-1)*pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*GammadNiOuter[0][1] + (1-((pow(gamma,p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*gammadd[0][1] + (-(p-1)*((pow(hyperplane[i],p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(hyperplane[i],2*p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiNiOuter[0][1] + ((p-1)*((pow(gamma,p-1)*pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiGammadOuter[0][1] + (1-((pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*0.0}, {(-(p-1)*((pow(gamma,p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(gamma,2*p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*GammadGammadOuter[1][0] + ((p-1)*((pow(gamma,p-1)*pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*GammadNiOuter[1][0] + (1-((pow(gamma,p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*gammadd[1][0] + (-(p-1)*((pow(hyperplane[i],p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(hyperplane[i],2*p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiNiOuter[1][0] + ((p-1)*((pow(gamma,p-1)*pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiGammadOuter[1][0] + (1-((pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*0.0, (-(p-1)*((pow(gamma,p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(gamma,2*p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*GammadGammadOuter[1][1] + ((p-1)*((pow(gamma,p-1)*pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*GammadNiOuter[1][1] + (1-((pow(gamma,p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*gammadd[1][1] + (-(p-1)*((pow(hyperplane[i],p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p)))+(p-1)*((pow(hyperplane[i],2*p-2))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiNiOuter[1][1] + ((p-1)*((pow(gamma,p-1)*pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),1+(p-1)/p))))*NiGammadOuter[1][1] + (1-((pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*0.0}};

        gammad = {(1-((pow(gamma,p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*gammad[0] + (1-((pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*NormalVectorTildeList[i][0], (1-((pow(gamma,p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*gammad[1] + (1-((pow(hyperplane[i],p-1))/(pow(pow(gamma,p)+pow(hyperplane[i],p),(p-1)/p))))*NormalVectorTildeList[i][1]};

        gamma = gamma + hyperplane[i] - pow(pow(gamma,p)+pow(hyperplane[i],p),1/p);
    }

    // Construct the output
    OutputStructScalar Output;
    Output.Value = gamma;
    Output.Gradient = gammad;
    Output.Hessian = gammadd;

    return Output;
}
