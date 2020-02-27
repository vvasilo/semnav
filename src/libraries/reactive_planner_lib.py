#!/usr/bin/env python

"""
MIT License (modified)

Copyright (c) 2020 The Trustees of the University of Pennsylvania
Authors:
Vasileios Vasilopoulos <vvasilo@seas.upenn.edu>
Omur Arslan <omur@seas.upenn.edu>

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
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import time
import itertools
import shapely as sp
from matplotlib.collections import PatchCollection
from shapely.geometry import Polygon, Point, LineString, LinearRing
from shapely.ops import cascaded_union
from scipy.spatial import ConvexHull
from scipy.signal import butter, lfilter
from operator import itemgetter

# Geometry imports
from polygeom_lib import cvxpolyxhplane, cvxpolyerode, polydist, polyxline, inpolygon, cvxpolyintersect, polyxray
from polygeom_lib import polytriangulation, polycvxdecomp, polyconvexdecomposition


class LIDARClass:
    """
    Class that describes a LIDAR object and is updated as new measurements are received

    Properties:
        1) RangeMeasurements: Range measurements received
        2) Range: Range of the sensor
        3) Infinity: Range to be considered as infinity
        4) MinAngle: Minimum angle of the sensor
        5) MaxAngle: Maximum angle of the sensor
        6) Resolution: Sensor angular resolution
    """
    def __init__(self, RangeMeasurements, Range, Infinity, MinAngle, MaxAngle, Resolution):
        self.RangeMeasurements = RangeMeasurements # 1D numpy array
        self.Range = Range # float
        self.Infinity = Infinity # float
        self.MinAngle = MinAngle # float
        self.MaxAngle = MaxAngle # float
        self.Resolution = Resolution # float
        self.NumSample = int(1+round((MaxAngle-MinAngle)/Resolution)) # integer
        self.Angle = np.linspace(MinAngle, MaxAngle, self.NumSample) # 1D numpy array


def completeLIDAR2D(LIDAR):
    """
    Function that completes missing LIDAR data due to limited field of view

    Input: 
        1) LIDAR: Incomplete LIDAR object

    Output:
        1) LIDAR: Complete modified LIDAR object
    """
    if ((LIDAR.MaxAngle-LIDAR.MinAngle) < 2*np.pi):
        tempR = LIDAR.RangeMeasurements
        tempAngle = LIDAR.Angle
        tempResolution = LIDAR.Resolution

        # Updated LIDAR model
        LIDAR.MaxAngle = LIDAR.MinAngle + 2*np.pi
        LIDAR.NumSample = int(round((2*np.pi/tempResolution)+1))
        LIDAR.Resolution = (LIDAR.MaxAngle-LIDAR.MinAngle)/(LIDAR.NumSample-1)
        LIDAR.Angle = np.linspace(LIDAR.MinAngle, LIDAR.MaxAngle, LIDAR.NumSample)

        # Completed Range Data
        R = LIDAR.Infinity*np.ones(LIDAR.NumSample)
        Indices = np.floor(((tempAngle-LIDAR.MinAngle+LIDAR.Resolution/2)%(2*np.pi))/LIDAR.Resolution)
        Indices = Indices.astype(int)
        R[Indices] = tempR
        R[R>LIDAR.Range] = LIDAR.Range
        LIDAR.RangeMeasurements = R
    
    return LIDAR


def constructLIDAR2D(DataLIDAR, CutoffRange, AllowableRange, Pitch=0.):
    """
    Function that constructs a LIDAR object to be used by the reactive planner
    
    Input:
        1) DataLIDAR: Received LIDAR data
        2) CutoffRange: Cutoff range below which the range measurement is ignored
        3) AllowableRange: Maximum allowable LIDAR range
        4) Pitch: Robot pitch to be considered for range compensation (default is 0)
    
    Output:
        1) LIDAR: Constructed LIDAR object
    """
    # LIDAR operations
    Range = AllowableRange
    Infinity = AllowableRange+10.
    MinAngle = float(DataLIDAR.angle_min)
    MaxAngle = float(DataLIDAR.angle_max)
    Resolution = float(DataLIDAR.angle_increment)
    RangeMeasurements = np.array(DataLIDAR.ranges)

    # Project on the robot plane
    RangeMeasurements = RangeMeasurements*np.cos(Pitch)

    # Reject NaN
    Inan = np.nonzero(np.isnan(RangeMeasurements))
    Inan = Inan[0]
    for k in Inan:
    	RangeMeasurements[k] = Infinity
    
    # Cutoff LIDAR range
    Icutoff = np.nonzero(RangeMeasurements <= CutoffRange)
    for k in Icutoff:
    	RangeMeasurements[k] = Infinity
    
    # Construct LIDAR object
    LIDAR = LIDARClass(RangeMeasurements, Range, Infinity, MinAngle, MaxAngle, Resolution)
    
    return LIDAR


def obstaclePointsLIDAR2D(RobotState, LIDAR):
    """
    Function that returns the coordinates of observed obstacle points from the LIDAR
    
    Input:
        1) RobotState: Robot position and orientation
        2) LIDAR: Current LIDAR object
    
    Output:
        1) PointsAll: Coordinates of the observed obstacle points from the LIDAR - Nx2 numpy.array
    """
    # Find robot position and orientation
    RobotPosition = RobotState[0:2]
    RobotOrientation = RobotState[2]

    # Range measurements
    R = np.array(LIDAR.RangeMeasurements)

    # Find observed LIDAR points coordinates
    PointsAll = np.zeros((R.shape[0], 2))
    for i in range(PointsAll.shape[0]):
        PointsAll[i][0] = RobotPosition[0]+R[i]*np.cos(LIDAR.Angle[i]+RobotOrientation)
        PointsAll[i][1] = RobotPosition[1]+R[i]*np.sin(LIDAR.Angle[i]+RobotOrientation)
    
    return PointsAll


def compensateObstacleLIDAR2D(RobotState, Obstacle, LIDAR):
    """
    Function that checks if the LIDAR hits a specific obstacle whose polygon is known

    Input:
        1) RobotState: Robot position and orientation
        2) Obstacle: shapely.geometry.Polygon obstacle defining the polygonal obstacle
        3) LIDAR: Current LIDAR object

    Output:
        1) LIDAR: Final LIDAR object
    """
    # Find robot position and orientation
    RobotPosition = RobotState[0:2]
    RobotOrientation = RobotState[2]
    
    # Find the indices that correspond to a LIDAR range less than the maximum range
    lidar_indices = np.where(LIDAR.RangeMeasurements < LIDAR.Range)
    for k in lidar_indices[0]:
        lidar_point_x = RobotPosition[0]+LIDAR.RangeMeasurements[k]*np.cos(LIDAR.Angle[k]+RobotOrientation)
        lidar_point_y = RobotPosition[1]+LIDAR.RangeMeasurements[k]*np.sin(LIDAR.Angle[k]+RobotOrientation)
        if Obstacle.contains(Point(lidar_point_x, lidar_point_y)):
            LIDAR.RangeMeasurements[k] = LIDAR.Range

    return LIDAR


def readLIDAR2D(RobotState,Obstacles,Range,MinAngle,MaxAngle,NumSample):
    """
    Function that generates a virtual LIDAR object, based on the position of the robot and the surrounding obstacles

    Input: 
        1) RobotState: Robot position and orientation
        2) Obstacles: shapely.geometry.Polygon obstacle array defining the polygonal obstacles
        3) Range: Range of the LIDAR object to be generated
        4) MinAngle: Minimum angle of the LIDAR object to be generated
        5) MaxAngle: Maximum angle of the LIDAR object to be generated
        6) NumSample: Number of samples used in the process

    Output:
        1) virtualLIDAR: Complete virtual LIDAR object
    """
    # Find robot position and orientation
    RobotPosition = RobotState[0:2]
    RobotOrientation = RobotState[2]

    # Initialize LIDAR object
    Resolution = (MaxAngle-MinAngle)/(NumSample-1)
    RangeMeasurements = np.zeros(NumSample)
    Infinity = Range + 20.0
    virtualLIDAR = LIDARClass(RangeMeasurements, Range, Infinity, MinAngle, MaxAngle, Resolution)

    # Rotation matrix from the global frame to the local sensor frame
    RotMat = np.array([[np.cos(-RobotOrientation),-np.sin(-RobotOrientation)],[np.sin(-RobotOrientation),np.cos(-RobotOrientation)]])

    # Determine distance to the workspace boundary and obstacles
    virtualLIDAR.RangeMeasurements = virtualLIDAR.Infinity*np.ones(virtualLIDAR.NumSample)
    for co in range(len(Obstacles)):
        # Obstacle in the local sensor frame
        Obs = np.vstack((Obstacles[co].exterior.coords.xy[0],Obstacles[co].exterior.coords.xy[1])).transpose() - RobotPosition
        Obs = np.matmul(Obs,RotMat.transpose())

        # Compute distance to every obstacle edge
        for cv in range(Obs.shape[0]):
            cn = ((cv+1)%Obs.shape[0])
            vc = Obs[cv][:] # current vertex
            vn = Obs[cn][:] # next vertex

            # Compute the distance to the origin
            dist = np.min([np.linalg.norm(vc),np.linalg.norm(vn)])
            w = (np.dot(vn,(vn-vc).transpose()))/(np.linalg.norm(vn-vc)**2)
            if (w>=0.0) and (w<=1.0):
                vx = w*vc + (1-w)*vn
                dist = np.min([dist, np.linalg.norm(vx)])
            
            ac = np.arctan2(vc[1],vc[0]) # Relative angle of the current vertex
            an = np.arctan2(vn[1],vn[0]) # Relative angle of the next vertex

            flagDist = (dist <= virtualLIDAR.Range)
            flagAngle = (np.min([ac,an]) <= virtualLIDAR.MaxAngle) and (np.max([ac,an]) >= virtualLIDAR.MinAngle)
            # Compute LIDAR range if the obstacle segment is in the sensing region
            if flagDist and flagAngle:
                # Closest LIDAR ray index
                I = int(round((np.max([np.min([ac,an]),virtualLIDAR.MinAngle])-virtualLIDAR.MinAngle)/virtualLIDAR.Resolution))
                I = (I%virtualLIDAR.NumSample)

                # Compute the intersection of the LIDAR ray with the sensor footprint
                vtemp = np.array([np.cos(virtualLIDAR.Angle[I]), np.sin(virtualLIDAR.Angle[I])])
                vRtemp = np.array([-np.sin(virtualLIDAR.Angle[I]), np.cos(virtualLIDAR.Angle[I])])
                w = -(np.dot(vn,vRtemp.transpose()))/(np.dot(vc-vn,vRtemp.transpose()))
                if (w>=0.0) and (w<=1.0):
                    xtemp = w*vc + (1-w)*vn
                    if (np.dot(xtemp,vtemp.transpose()) >= 0):
                        virtualLIDAR.RangeMeasurements[I] = np.min([virtualLIDAR.RangeMeasurements[I],np.linalg.norm(xtemp)])
                
                # Compute the intersection of adjacent LIDAR rays
                J = ((I+1)%virtualLIDAR.NumSample)
                flagValid = True
                while flagValid and (J is not I):
                    vtemp = np.array([np.cos(virtualLIDAR.Angle[J]),np.sin(virtualLIDAR.Angle[J])])
                    vRtemp = np.array([-np.sin(virtualLIDAR.Angle[J]),np.cos(virtualLIDAR.Angle[J])])
                    w = -(np.dot(vn,vRtemp.transpose()))/(np.dot(vc-vn,vRtemp.transpose()))
                    if (w>=0.0) and (w<=1.0):
                        xtemp = w*vc + (1-w)*vn
                        if (np.dot(xtemp,vtemp.transpose()) >= 0):
                            virtualLIDAR.RangeMeasurements[J] = np.min([virtualLIDAR.RangeMeasurements[J],np.linalg.norm(xtemp)])
                            J = ((J+1)%virtualLIDAR.NumSample)
                        else:
                            flagValid = False
                    else:
                        flagValid = False
                
                J = (I-1)%virtualLIDAR.NumSample
                flagValid = True
                while flagValid and (J is not I):
                    vtemp = np.array([np.cos(virtualLIDAR.Angle[J]),np.sin(virtualLIDAR.Angle[J])])
                    vRtemp = np.array([-np.sin(virtualLIDAR.Angle[J]),np.cos(virtualLIDAR.Angle[J])])
                    w = -(np.dot(vn,vRtemp.transpose()))/(np.dot(vc-vn,vRtemp.transpose()))
                    if (w>=0.0) and (w<=1.0):
                        xtemp = w*vc + (1-w)*vn
                        if (np.dot(xtemp,vtemp.transpose()) >= 0):
                            virtualLIDAR.RangeMeasurements[J] = np.min([virtualLIDAR.RangeMeasurements[J],np.linalg.norm(xtemp)])
                            J = (J-1)%virtualLIDAR.NumSample
                        else:
                            flagValid = False
                    else:
                        flagValid = False
    
    # Check sensor range
    virtualLIDAR.RangeMeasurements[virtualLIDAR.RangeMeasurements > virtualLIDAR.Range] = virtualLIDAR.Range
    
    return virtualLIDAR


def translateLIDAR2D(RobotState, RobotStateTransformed, RobotRadius, LIDAR):
    """
    Rebase LIDAR readings from RobotState to RobotStateTransformed

    Input:
        1) RobotState: Original robot position and orientation
        2) RobotStateTransformed: Transformed robot position and orientation
        3) RobotRadius: Robot radius
        4) LIDAR: Current LIDAR object

    Output:
        1) newLIDAR: Transformed LIDAR object
    """
    # Find robot position and orientation
    RobotPosition = RobotState[0:2]
    RobotOrientation = RobotState[2]

    # Find transformed robot position and orientation
    RobotPositionTransformed = RobotStateTransformed[0:2]
    RobotOrientationTransformed = RobotStateTransformed[2]

    # Account for the robot radius
    LIDAR.RangeMeasurements = LIDAR.RangeMeasurements - RobotRadius
    
    # Find obstacle points
    obstacle_points = obstaclePointsLIDAR2D(RobotState, LIDAR) - RobotPosition

    # Rotation matrix from the global frame to the local sensor frame
    RotMat = np.array([[np.cos(-RobotOrientation),-np.sin(-RobotOrientation)],[np.sin(-RobotOrientation),np.cos(-RobotOrientation)]])

    # Points in the local sensor frame
    obstacle_points = np.matmul(obstacle_points,RotMat.transpose())

    # Rotation matrix from the local sensor frame to the global frame in the transformed space
    RotMatTransformed = np.array([[np.cos(RobotOrientationTransformed),-np.sin(RobotOrientationTransformed)],[np.sin(RobotOrientationTransformed),np.cos(RobotOrientationTransformed)]])

    # Points in the transformed space
    obstacle_points_transformed = np.matmul(obstacle_points,RotMatTransformed.transpose())
    obstacle_points_transformed = obstacle_points_transformed + RobotPositionTransformed

    # Find new ranges
    obstacle_points_transformedT = obstacle_points_transformed.transpose()
    obstacle_points_transformed_x = obstacle_points_transformedT[0][:]
    obstacle_points_transformed_y = obstacle_points_transformedT[1][:]
    R = np.sqrt((obstacle_points_transformed_x-RobotPositionTransformed[0])**2+(obstacle_points_transformed_y-RobotPositionTransformed[1])**2)
    newLIDAR = LIDARClass(R, LIDAR.Range-np.linalg.norm(RobotPositionTransformed-RobotPosition), LIDAR.Infinity, LIDAR.MinAngle, LIDAR.MaxAngle, LIDAR.Resolution)

    return newLIDAR


def localminLIDAR2D(LIDAR):
    """
    Function that finds the indices of local minima of the LIDAR range data
    
    Input:
        1) LIDAR: Current LIDAR object
    
    Output:
        1) Imin: Indices of local minima of the LIDAR range data
    """
    R = LIDAR.RangeMeasurements

    # Compute the indices of strictly local minima of the LIDAR range data
    if ((LIDAR.MaxAngle-LIDAR.MinAngle)<2*np.pi):
        # Assume that the region outside the angular range of LIDAR is free
        Rp = np.hstack((np.array([LIDAR.Range]), R[0:-1]))
        Rn = np.hstack((R[1:], np.array([LIDAR.Range])))
    else:
        Rp = np.hstack((np.array([R[-2]]), R[0:-1]))
        Rn = np.hstack((R[1:], np.array([R[1]])))
    
    # Logical tests
    logical_test = np.logical_or(np.logical_and(R<=Rp, R<Rn), np.logical_and(R<Rp, R<=Rn))
    Imin = logical_test+0
    return Imin


def localworkspaceLIDAR2D(RobotState, RobotRadius, LIDAR):
    """
    Function that returns the local workspace
    
    Input:
        1) RobotState: Robot position and orientation
        2) RobotRadius: Robot radius
        3) LIDAR: Current LIDAR object
    
    Output:
        1) LW: Local workspace polygon array
    """
    X = RobotState
    epsilon = 0.000000001

    R = LIDAR.RangeMeasurements
    if R.min(0)<epsilon:
        LW = np.array([[]])
        return LW
    else:
        # Complete missing data due to the LIDAR's angular range limits
        LIDAR = completeLIDAR2D(LIDAR)

        # Modified range data defining the local workspace
        R = 0.5*(LIDAR.RangeMeasurements+RobotRadius)

        # Initialize the local workspace with the minimum square that respects
        # the LIDAR's sensing range
        LW = (0.5*(LIDAR.Range+RobotRadius))*np.array([[-1,-1], [-1,1], [1,1], [1,-1]])
        Imin = np.nonzero(localminLIDAR2D(LIDAR))
        Imin = Imin[0]

        for k in Imin:
            if not LW.any():
                return LW
            else:
                # Local minimum parameters
                Ak = LIDAR.Angle[k] # Angle
                Rk = R[k] # Range

                # Separating hyperplane parameters
                n = np.array([-np.cos(Ak+X[2]), -np.sin(Ak+X[2])]) # Hyperplane normal
                m = -Rk*n # A point on the separating hyperplane

                # Update the local workspace by taking its intersection with the associated halfplane
                LW = cvxpolyxhplane(LW, m, n)
        
        # Local workspace footprint
        LocalFootprint = np.vstack((R*np.cos(LIDAR.Angle+X[2]), R*np.sin(LIDAR.Angle+X[2])))
        LocalFootprint = LocalFootprint.transpose()

        # Update local workspace
        if Polygon(LW).is_valid and Polygon(LocalFootprint).is_valid:
            LW = cvxpolyintersect(LW, LocalFootprint)
            if LW.any():
                LW = LW + np.array([[X[0], X[1]]])
            else:
                LW = np.array([[]])
                return LW

            # Make local workspace convex
            convhullind = ConvexHull(LW)
            LW = LW[convhullind.vertices]
            return LW
        else:
            LW = np.array([[]])
            return LW


def localfreespaceLIDAR2D(RobotState, RobotRadius, LIDAR):
    """
    Function that returns the local freespace
    
    Input: 
        1) RobotState: Robot position and orientation
        2) RobotRadius: Robot radius
        3) LIDAR: Current LIDAR object
    
    Output: 
        1) LF: Local freespace polygon array
    """
    X = RobotState
    epsilon = 0.000000001

    R = LIDAR.RangeMeasurements
    if R.min(0)<epsilon:
        LF = np.array([[]])
        return LF
    else:
        # Complete missing data due to the LIDAR's angular range limits
        LIDAR = completeLIDAR2D(LIDAR)

        # Modified range data defining the local freespace
        R = 0.5*(LIDAR.RangeMeasurements-RobotRadius)

        # Initialize the local freespace with the minimum square that respects
        # the LIDAR's sensing range
        LF = (0.5*(LIDAR.Range-RobotRadius))*np.array([[-1,-1], [-1,1], [1,1], [1,-1]])
        Imin = np.nonzero(localminLIDAR2D(LIDAR))
        Imin = Imin[0]

        for k in Imin:
            if not LF.any():
                return LF
            else:
                # Local minimum parameters
                Ak = LIDAR.Angle[k] # Angle
                Rk = R[k] # Range

                # Separating hyperplane parameters
                n = np.array([-np.cos(Ak+X[2]), -np.sin(Ak+X[2])]) # Hyperplane normal
                m = -Rk*n # A point on the separating hyperplane

                # Update the local freespace by taking its intersection with the associated halfplane
                LF = cvxpolyxhplane(LF, m, n)

        LocalFootprint = np.vstack((R*np.cos(LIDAR.Angle+X[2]), R*np.sin(LIDAR.Angle+X[2])))
        LocalFootprint = LocalFootprint.transpose()

        # Update local freespace
        if Polygon(LF).is_valid and Polygon(LocalFootprint).is_valid:
            LF = cvxpolyintersect(LF, LocalFootprint)
            if LF.any():
                LF = LF + np.array([[X[0], X[1]]])
            else:
                LF = np.array([[]])
                return LF

            # Make local freespace convex
            convhullind = ConvexHull(LF)
            LF = LF[convhullind.vertices]
            return LF
        else:
            LF = np.array([[]])
            return LF


def localfreespace_linearLIDAR2D(RobotState, LF):
    """
    Function that returns the linear local freespace as the intersection of the local freespace with the current heading line
    
    Input:
        1) RobotState: Robot position and orientation
        2) LF: Local freespace
    
    Output:
        1) LFL: Linear freespace
    """
    X = RobotState
    
    if not LF.any():
        LFL = np.array([[]])
        return LFL
    else:
        RobotPosition = np.array([X[0], X[1]])
        RobotDirection = np.array([np.cos(X[2]), np.sin(X[2])])
        LFL = polyxray(LF, RobotPosition, RobotDirection)
        LFL = np.vstack((RobotPosition, LFL))
        return LFL


def localfreespace_angularLIDAR2D(RobotState, LF, Goal):
    """
    Function that returns the angular local freespace as the intersection of the local freespace with the line connecting the robot to the goal
    
    Input:
        1) RobotState: Robot position and orientation
        2) LF: Local freespace

    Output:
        1) LFA: Angular freespace
    """
    X = RobotState

    if not LF.any():
        LFA = np.array([[]])
        return LFA
    else:
        RobotPosition = np.array([X[0], X[1]])
        GoalDirection = np.array([Goal[0]-RobotPosition[0], Goal[1]-RobotPosition[1]])
        LFA = polyxray(LF, RobotPosition, GoalDirection)
        LFA = np.vstack((RobotPosition, LFA))
        return LFA


def localgoalLIDAR2D(LF, Goal):
    """
    Function that computes the local goal as the projection of the global goal on the local freespace
    
    Input:
        1) LF: Local freespace
        2) Goal: Global goal
    
    Output:
        1) LGA1: Local goal
    """
    # Compute local goal --- the closest point of local free space to the global goal
    if not LF.any():
        LGA1 = np.array([Goal])
    else:
        if inpolygon(LF, Goal):
            LGA1 = np.array([Goal])
        else:
            D, LGA1 = polydist(LF, Goal)
    return LGA1


def localgoal_linearLIDAR2D(RobotState, LF, Goal):
    """
    Function that computes the linear local goal as the projection of the global goal on the linear local freespace
    
    Input:
        1) RobotState: Robot position and orientation
        2) LF: Local freespace
        3) Goal: Global goal
    
    Output:
        1) LGL: Local linear goal
    """
    # Compute linear local free space
    LFL = localfreespace_linearLIDAR2D(RobotState, LF)

    # Compute local goal for unicycle
    if not LFL.any():
        LGL = np.array([[RobotState[0], RobotState[1]]])
    else:
        D, LGL = polydist(LFL, Goal)
    return LGL


def localgoal_angularLIDAR2D(RobotState, LF, Goal):
    """
    Function that computes the angular local goal as the projection of the global goal on the angular local freespace
    
    Input:
        1) RobotState: Robot position and orientation
        2) LF: Local freespace
        3) Goal: Global goal
    
    Output:
        1) LGA2: Local angular goal
    """
    # Compute angular local free space
    LFA = localfreespace_angularLIDAR2D(RobotState, LF, Goal)

    # Compute local goal for unicycle
    if not LFA.any():
        LGA2 = np.array([[Goal[0], Goal[1]]])
    else:
        D, LGA2 = polydist(LFA, Goal)
    return LGA2


def diffeoTreeTriangulation(PolygonVertices, DiffeoParams):
    """
    Function that calculates the triangulation tree of a polygon and augments it with properties used in semantic navigation

    Input:
        1) PolygonVertices: Vertex Coordinates of input polygon - Nx2 numpy.array (start and end vertices must be the same)
        2) DiffeoParams: Options for the diffeomorphism construction
    
    Output:
        1) tree: Modified tree with added properties
                 a) For the root:
                                 i) 'radius': the radius of the final sphere to be constructed
                 b) For the children:
                                 i) 'r_center_t': The tangents from vertices 0 and 1 to the center - 2x2 numpy array
                                 2) 'r_center_n': The normals corresponding to 'r_center_t'
                 c) For the root and the children:
                                     i) 'vertices': the vertices of the triangle - 3x2 numpy array in CCW order
                                     ii) 'vertices_tilde': the vertices of the polygonal collar that encompasses the triangle - Mx2 numpy array in CCW order starting from the center in the parent
                                     iii) 'r_t': the unit tangents for the triangle to be deformed in a 3x2 numpy array in CCW order
                                                 (the 1st row is the tangent shared between the parent and the child)
                                     iv) 'r_n': the unit normals for the triangle to be deformed in a 3x2 array corresponding to 'r_t'
                                     v) 'r_tilde_t': the unit tangents for the polygonal collar in a Mx2 numpy array in CCW order
                                     vi) 'r_tilde_n': the unit normals for the polygonal collar in a Mx2 numpy array corresponding to 'r_tilde_t'
                                     vii) 'center': the center in the parent, used for the purging transformation, or the center of the root used for the final transformation
    """
    # Unpack diffeomorphism parameters
    varepsilon = DiffeoParams['varepsilon']
    workspace = DiffeoParams['workspace']

    # Check if the polygon intersects the workspace boundary
    if Polygon(PolygonVertices).intersects(LineString(workspace)):
        # Compute the intersection with the workspace
        polygon_to_use = sp.geometry.polygon.orient(Polygon(PolygonVertices).intersection(Polygon(workspace)), 1.0)
        PolygonVertices = np.vstack((polygon_to_use.exterior.coords.xy[0], polygon_to_use.exterior.coords.xy[1])).transpose()
        
        # Compute the triangulation tree of the polygon with its dual (adjacency) graph
        tree = polytriangulation(PolygonVertices, workspace, True)

        # Find the center and the adjacency edge to the boundary
        D, C = polydist(workspace,tree[-1]['vertices'])
        inds = D.argsort()
        tree[-1]['vertices'] = tree[-1]['vertices'][inds]
        root_polygon_coords = np.vstack((tree[-1]['vertices'], tree[-1]['vertices'][0]))
        if not LinearRing(root_polygon_coords).is_ccw:
            tree[-1]['vertices'][[0,1]] = tree[-1]['vertices'][[1,0]]
        tree[-1]['adj_edge'] = np.vstack((tree[-1]['vertices'][0], tree[-1]['vertices'][1]))
        median_point = 0.5*np.array([[tree[-1]['adj_edge'][1][0]+tree[-1]['adj_edge'][0][0], tree[-1]['adj_edge'][1][1]+tree[-1]['adj_edge'][0][1]]])
        median_ray = np.array([[median_point[0][0]-tree[-1]['vertices'][2][0], median_point[0][1]-tree[-1]['vertices'][2][1]]])
        median_ray = median_ray/np.linalg.norm(median_ray[0])
        tree[-1]['center'] = np.array([[median_point[0][0]+1.0*median_ray[0][0], median_point[0][1]+1.0*median_ray[0][1]]])

        # Compute the tangent and normal vectors of the root triangle
        tree[-1]['r_t'] = np.array(tree[-1]['vertices'][1]-tree[-1]['vertices'][0])/np.linalg.norm(tree[-1]['vertices'][1]-tree[-1]['vertices'][0])
        tree[-1]['r_t'] = np.vstack((tree[-1]['r_t'], np.array(tree[-1]['vertices'][2]-tree[-1]['vertices'][1])/np.linalg.norm(tree[-1]['vertices'][2]-tree[-1]['vertices'][1])))
        tree[-1]['r_t'] = np.vstack((tree[-1]['r_t'], np.array(tree[-1]['vertices'][0]-tree[-1]['vertices'][2])/np.linalg.norm(tree[-1]['vertices'][0]-tree[-1]['vertices'][2])))
        tree[-1]['r_n'] = np.array([-tree[-1]['r_t'][0][1], tree[-1]['r_t'][0][0]])
        tree[-1]['r_n'] = np.vstack((tree[-1]['r_n'], np.array([-tree[-1]['r_t'][1][1],tree[-1]['r_t'][1][0]])))
        tree[-1]['r_n'] = np.vstack((tree[-1]['r_n'], np.array([-tree[-1]['r_t'][2][1],tree[-1]['r_t'][2][0]])))

        # Find the remaining tangents and normals from vertices 0 and 1 to the center
        tree[-1]['r_center_t'] = (tree[-1]['center'][0]-tree[-1]['vertices'][0])/np.linalg.norm(tree[-1]['center'][0]-tree[-1]['vertices'][0])
        tree[-1]['r_center_n'] = np.array([-tree[-1]['r_center_t'][1], tree[-1]['r_center_t'][0]])
        tree[-1]['r_center_t'] = np.vstack((tree[-1]['r_center_t'], (tree[-1]['vertices'][1]-tree[-1]['center'][0])/np.linalg.norm(tree[-1]['vertices'][1]-tree[-1]['center'][0])))
        tree[-1]['r_center_n'] = np.vstack((tree[-1]['r_center_n'], np.array([-tree[-1]['r_center_t'][1][1],tree[-1]['r_center_t'][1][0]])))

        # Compute the dilated polygon and truncate it by the rays emanating from the center
        original_polygon = np.array([tree[-1]['center'][0], tree[-1]['vertices'][1], tree[-1]['vertices'][2], tree[-1]['vertices'][0], tree[-1]['center'][0]])
        polygon_tilde = sp.geometry.polygon.orient(Polygon(original_polygon).buffer(varepsilon, join_style=1), 1.0)
        dilation = np.vstack((polygon_tilde.exterior.coords.xy[0], polygon_tilde.exterior.coords.xy[1])).transpose()
        intersect_1 = cvxpolyxhplane(dilation[0:-1], tree[-1]['center'][0], tree[-1]['r_center_n'][0])
        intersect_2 = cvxpolyxhplane(intersect_1, tree[-1]['center'][0], tree[-1]['r_center_n'][1])
        polygon_tilde_vertices = np.vstack((intersect_2,intersect_2[0]))

        # Compute the intersection with the workspace
        final_polygon = sp.geometry.polygon.orient(Polygon(polygon_tilde_vertices).intersection(Polygon(workspace).union(Polygon(np.array([tree[-1]['center'][0], tree[-1]['vertices'][1], tree[-1]['vertices'][2], tree[-1]['vertices'][0], tree[-1]['center'][0]])))).simplify(0.01), 1.0)
        tree[-1]['vertices_tilde'] = np.vstack((final_polygon.exterior.coords.xy[0][0:-1], final_polygon.exterior.coords.xy[1][0:-1])).transpose()

        # Find the tangent and normal vectors for the generated polygonal collar
        vertices_to_consider = np.vstack((tree[-1]['vertices_tilde'],tree[-1]['vertices_tilde'][0]))
        tree[-1]['r_tilde_t'] = (vertices_to_consider[1]-vertices_to_consider[0])/np.linalg.norm(vertices_to_consider[1]-vertices_to_consider[0])
        tree[-1]['r_tilde_n'] = np.array([-tree[-1]['r_tilde_t'][1],tree[-1]['r_tilde_t'][0]])
        for j in range(1,vertices_to_consider.shape[0]-1):
            tree[-1]['r_tilde_t'] = np.vstack((tree[-1]['r_tilde_t'],(vertices_to_consider[j+1]-vertices_to_consider[j])/np.linalg.norm(vertices_to_consider[j+1]-vertices_to_consider[j])))
            tree[-1]['r_tilde_n'] = np.vstack((tree[-1]['r_tilde_n'],np.array([-tree[-1]['r_tilde_t'][j][1],tree[-1]['r_tilde_t'][j][0]])))
        
        # Add a dummy radius
        tree[-1]['radius'] = 0.0
    else:
        # Compute the triangulation tree of the polygon with its dual (adjacency) graph
        tree = polytriangulation(PolygonVertices, workspace, False)

        # Start with the root and find the center and the radius
        root_coords = tree[-1]['vertices'].transpose()
        tree[-1]['center'] = np.array([[sum(root_coords[0])/len(root_coords[0]), sum(root_coords[1])/len(root_coords[1])]])
        D, closest_point = polydist(tree[-1]['vertices'], tree[-1]['center'])
        tree[-1]['radius'] = 0.8*D[0]

        # Compute the tangent and normal vectors of the root triangle
        tree[-1]['r_t'] = np.array(tree[-1]['vertices'][1]-tree[-1]['vertices'][0])/np.linalg.norm(tree[-1]['vertices'][1]-tree[-1]['vertices'][0])
        tree[-1]['r_t'] = np.vstack((tree[-1]['r_t'], np.array(tree[-1]['vertices'][2]-tree[-1]['vertices'][1])/np.linalg.norm(tree[-1]['vertices'][2]-tree[-1]['vertices'][1])))
        tree[-1]['r_t'] = np.vstack((tree[-1]['r_t'], np.array(tree[-1]['vertices'][0]-tree[-1]['vertices'][2])/np.linalg.norm(tree[-1]['vertices'][0]-tree[-1]['vertices'][2])))
        tree[-1]['r_n'] = np.array([-tree[-1]['r_t'][0][1], tree[-1]['r_t'][0][0]])
        tree[-1]['r_n'] = np.vstack((tree[-1]['r_n'], np.array([-tree[-1]['r_t'][1][1],tree[-1]['r_t'][1][0]])))
        tree[-1]['r_n'] = np.vstack((tree[-1]['r_n'], np.array([-tree[-1]['r_t'][2][1],tree[-1]['r_t'][2][0]])))

        # Find the polygonal collar for the root by dilating the triangle by varepsilon
        polygon_tilde = sp.geometry.polygon.orient(Polygon(tree[-1]['vertices']).buffer(varepsilon, join_style=1).intersection(Polygon(workspace)).simplify(0.01), 1.0)
        tree[-1]['vertices_tilde'] = np.vstack((polygon_tilde.exterior.coords.xy[0][0:-1], polygon_tilde.exterior.coords.xy[1][0:-1])).transpose()

        # Find the tangent and normal vectors for the generated polygonal collar
        vertices_to_consider = np.vstack((tree[-1]['vertices_tilde'],tree[-1]['vertices_tilde'][0]))
        tree[-1]['r_tilde_t'] = (vertices_to_consider[1]-vertices_to_consider[0])/np.linalg.norm(vertices_to_consider[1]-vertices_to_consider[0])
        tree[-1]['r_tilde_n'] = np.array([-tree[-1]['r_tilde_t'][1],tree[-1]['r_tilde_t'][0]])
        for j in range(1,vertices_to_consider.shape[0]-1):
            tree[-1]['r_tilde_t'] = np.vstack((tree[-1]['r_tilde_t'],(vertices_to_consider[j+1]-vertices_to_consider[j])/np.linalg.norm(vertices_to_consider[j+1]-vertices_to_consider[j])))
            tree[-1]['r_tilde_n'] = np.vstack((tree[-1]['r_tilde_n'],np.array([-tree[-1]['r_tilde_t'][j][1],tree[-1]['r_tilde_t'][j][0]])))

    # Identify all the children properties
    for i in range(0,len(tree)-1):
        # Compute the tangent and normal vectors of the child hyperplanes
        # r0 is always the shared edge between the parent and the child, r1 and r2 the rest in CCW order
        tree[i]['r_t'] = np.array(tree[i]['vertices'][1]-tree[i]['vertices'][0])/np.linalg.norm(tree[i]['vertices'][1]-tree[i]['vertices'][0])
        tree[i]['r_t'] = np.vstack((tree[i]['r_t'], np.array(tree[i]['vertices'][2]-tree[i]['vertices'][1])/np.linalg.norm(tree[i]['vertices'][2]-tree[i]['vertices'][1])))
        tree[i]['r_t'] = np.vstack((tree[i]['r_t'], np.array(tree[i]['vertices'][0]-tree[i]['vertices'][2])/np.linalg.norm(tree[i]['vertices'][0]-tree[i]['vertices'][2])))
        tree[i]['r_n'] = np.array([-tree[i]['r_t'][0][1], tree[i]['r_t'][0][0]])
        tree[i]['r_n'] = np.vstack((tree[i]['r_n'], np.array([-tree[i]['r_t'][1][1],tree[i]['r_t'][1][0]])))
        tree[i]['r_n'] = np.vstack((tree[i]['r_n'], np.array([-tree[i]['r_t'][2][1],tree[i]['r_t'][2][0]])))

        # Find the median from the 3rd point to the shared edge and from that compute the center for the purging transformation
        median_point = 0.5*np.array([[tree[i]['adj_edge'][1][0]+tree[i]['adj_edge'][0][0], tree[i]['adj_edge'][1][1]+tree[i]['adj_edge'][0][1]]])
        median_ray = np.array([[median_point[0][0]-tree[i]['vertices'][2][0], median_point[0][1]-tree[i]['vertices'][2][1]]])
        median_ray = median_ray/np.linalg.norm(median_ray[0])
        intersection_point = polyxray(tree[tree[i]['predecessor']]['vertices'], median_point[0], median_ray[0]) # offset median point by a little bit to avoid numerical problems
        tree[i]['center'] = np.array([[0.2*median_point[0][0]+0.8*intersection_point[0], 0.2*median_point[0][1]+0.8*intersection_point[1]]])

        # Find the remaining tangents and normals from vertices 0 and 1 to the center
        tree[i]['r_center_t'] = (tree[i]['center'][0]-tree[i]['vertices'][0])/np.linalg.norm(tree[i]['center'][0]-tree[i]['vertices'][0])
        tree[i]['r_center_n'] = np.array([-tree[i]['r_center_t'][1], tree[i]['r_center_t'][0]])
        tree[i]['r_center_t'] = np.vstack((tree[i]['r_center_t'], (tree[i]['vertices'][1]-tree[i]['center'][0])/np.linalg.norm(tree[i]['vertices'][1]-tree[i]['center'][0])))
        tree[i]['r_center_n'] = np.vstack((tree[i]['r_center_n'], np.array([-tree[i]['r_center_t'][1][1],tree[i]['r_center_t'][1][0]])))

        # Compute the dilated polygon and truncate it by the rays emanating from the center
        original_polygon = np.array([tree[i]['center'][0], tree[i]['vertices'][1], tree[i]['vertices'][2], tree[i]['vertices'][0], tree[i]['center'][0]])
        polygon_tilde = sp.geometry.polygon.orient(Polygon(original_polygon).buffer(varepsilon, join_style=1).simplify(0.01), 1.0)
        dilation = np.vstack((polygon_tilde.exterior.coords.xy[0], polygon_tilde.exterior.coords.xy[1])).transpose()
        intersect_1 = cvxpolyxhplane(dilation[0:-1], tree[i]['center'][0], tree[i]['r_center_n'][0])
        intersect_2 = cvxpolyxhplane(intersect_1, tree[i]['center'][0], tree[i]['r_center_n'][1])
        candidate_polygon_vertices = np.vstack((intersect_2,intersect_2[0]))
        candidate_polygon = Polygon(candidate_polygon_vertices)

        # Check for collisions with all the triangles that will succeed i in the diffeomorphism construction except for its parent
        for j in range(i+1,len(tree)):
            if (j == tree[i]['predecessor']):
                continue
            else:
                polygon_to_test = Polygon(tree[j]['vertices'])
                candidate_polygon = (candidate_polygon.buffer(0)).difference(polygon_to_test.buffer(0))
                # If the difference operation created a multipolygon, keep only the polygon that contains the barycenter of the extended triangle
                if candidate_polygon.geom_type == 'MultiPolygon':
                    point_to_consider = Point((tree[i]['vertices'][0][0]+tree[i]['vertices'][1][0]+tree[i]['center'][0][0])/3.0, (tree[i]['vertices'][0][1]+tree[i]['vertices'][1][1]+tree[i]['center'][0][1])/3.0)
                    for k in range(len(candidate_polygon)):
                        if candidate_polygon[k].contains(point_to_consider):
                            candidate_polygon = candidate_polygon[k]
                            break
        
        # Extract final vertices
        candidate_polygon = sp.geometry.polygon.orient(candidate_polygon.simplify(0.01), 1.0)
        candidate_polygon_vertices = np.vstack((candidate_polygon.exterior.coords.xy[0], candidate_polygon.exterior.coords.xy[1])).transpose()

        # Decompose the polygon into its convex pieces and find the piece that includes the barycenter of the extended triangle
        decomposition = polycvxdecomp(candidate_polygon_vertices.tolist())
        for j in range(len(decomposition)):
            point_to_consider = Point((tree[i]['vertices'][0][0]+tree[i]['vertices'][1][0]+tree[i]['center'][0][0])/3.0, (tree[i]['vertices'][0][1]+tree[i]['vertices'][1][1]+tree[i]['center'][0][1])/3.0)
            polygon_to_consider = Polygon(decomposition[j])
            if polygon_to_consider.buffer(0.01).contains(point_to_consider):
                final_polygon_vertices = np.vstack((polygon_to_consider.exterior.coords.xy[0], polygon_to_consider.exterior.coords.xy[1])).transpose()
                break

        # Generate the outer polygonal collar
        final_polygon = sp.geometry.polygon.orient(Polygon(final_polygon_vertices).intersection(Polygon(workspace)), 1.0)
        tree[i]['vertices_tilde'] = np.vstack((final_polygon.exterior.coords.xy[0][0:-1], final_polygon.exterior.coords.xy[1][0:-1])).transpose()

        # Find the tangent and normal vectors for the generated polygonal collar
        vertices_to_consider = np.vstack((tree[i]['vertices_tilde'],tree[i]['vertices_tilde'][0]))
        tree[i]['r_tilde_t'] = (vertices_to_consider[1]-vertices_to_consider[0])/np.linalg.norm(vertices_to_consider[1]-vertices_to_consider[0])
        tree[i]['r_tilde_n'] = np.array([-tree[i]['r_tilde_t'][1],tree[i]['r_tilde_t'][0]])
        for j in range(1,vertices_to_consider.shape[0]-1):
            tree[i]['r_tilde_t'] = np.vstack((tree[i]['r_tilde_t'],(vertices_to_consider[j+1]-vertices_to_consider[j])/np.linalg.norm(vertices_to_consider[j+1]-vertices_to_consider[j])))
            tree[i]['r_tilde_n'] = np.vstack((tree[i]['r_tilde_n'],np.array([-tree[i]['r_tilde_t'][j][1],tree[i]['r_tilde_t'][j][0]])))
        
    return tree


def diffeoTreeConvex(PolygonVertices, DiffeoParams):
    """
    Function that calculates the convex decomposition of a polygon and augments it with properties used in semantic navigation

    Input:
        1) PolygonVertices: Vertex Coordinates of input polygon - Nx2 numpy.array (start and end vertices must be the same)
        2) DiffeoParams: Options for the diffeomorphism construction
    
    Output:
        1) tree: Modified tree with added properties
                 a) For the root:
                                 i) 'radius': the radius of the final sphere to be constructed
                 b) For the children:
                                 i) 'r_center_t': The tangents from vertices 0 and 1 to the center - 2x2 numpy array
                                 2) 'r_center_n': The normals corresponding to 'r_center_t'
                 c) For the root and the children:
                                     i) 'vertices': the vertices of the polygon - Nx2 numpy array in CCW order
                                     ii) 'augmented_vertices': the vertices of the polygon including the center of deformation (the second element in this array) - (N+1)x2 numpy array in CCW order
                                     ii) 'vertices_tilde': the vertices of the polygonal collar that encompasses the polygon - Mx2 numpy array in CCW order starting from the center in the parent
                                     iii) 'r_t': the unit tangents corresponding to augmented_vertices in CCW order
                                     iv) 'r_n': the unit normals for the polygon to be deformed in an array corresponding to 'r_t'
                                     v) 'r_tilde_t': the unit tangents for the polygonal collar in a Mx2 numpy array in CCW order
                                     vi) 'r_tilde_n': the unit normals for the polygonal collar in a Mx2 numpy array corresponding to 'r_tilde_t'
                                     vii) 'center': the center in the parent, used for the purging transformation, or the center of the root used for the final transformation
    """
    # Unpack diffeomorphism parameters
    varepsilon = DiffeoParams['varepsilon']
    workspace = DiffeoParams['workspace']

    # Check if the polygon intersects the workspace boundary
    if Polygon(PolygonVertices).intersects(LineString(workspace)):
        # Compute the intersection with the workspace
        polygon_to_use = sp.geometry.polygon.orient(Polygon(PolygonVertices).intersection(Polygon(workspace)), 1.0)
        PolygonVertices = np.vstack((polygon_to_use.exterior.coords.xy[0], polygon_to_use.exterior.coords.xy[1])).transpose()
        
        # Compute the convex decomposition tree of the polygon with its dual (adjacency) graph
        tree = polyconvexdecomposition(PolygonVertices, workspace, True)

        # Find the center and the adjacency edge to the boundary
        D, C = polydist(workspace,tree[-1]['vertices'])
        inds = D.argsort()
        if D[(inds[0]+1)%tree[-1]['vertices'].shape[0]] >= D[(inds[0]-1)%tree[-1]['vertices'].shape[0]]:
            tree[-1]['vertices'] = np.roll(tree[-1]['vertices'],-(inds[0]-1)%tree[-1]['vertices'].shape[0],axis=0)
        else:
            tree[-1]['vertices'] = np.roll(tree[-1]['vertices'],-(inds[0])%tree[-1]['vertices'].shape[0],axis=0)
        tree[-1]['adj_edge'] = np.vstack((tree[-1]['vertices'][0], tree[-1]['vertices'][1]))
        median_point = 0.5*np.array([[tree[-1]['adj_edge'][1][0]+tree[-1]['adj_edge'][0][0], tree[-1]['adj_edge'][1][1]+tree[-1]['adj_edge'][0][1]]])
        median_ray = np.array([[median_point[0][0]-tree[-1]['vertices'][2][0], median_point[0][1]-tree[-1]['vertices'][2][1]]])
        median_ray = median_ray/np.linalg.norm(median_ray[0])
        tree[-1]['center'] = np.array([[median_point[0][0]+0.3*median_ray[0][0], median_point[0][1]+0.3*median_ray[0][1]]])

        # Compute the tangent and normal vectors of the child hyperplanes
        # r0 is always the shared edge between the parent and the child, the rest in CCW order
        tree[-1]['r_t'] = []
        for j in range(0,tree[-1]['vertices'].shape[0]):
            tree[-1]['r_t'].append(np.array(tree[-1]['vertices'][(j+1)%tree[-1]['vertices'].shape[0]]-tree[-1]['vertices'][j%tree[-1]['vertices'].shape[0]])/np.linalg.norm(tree[-1]['vertices'][(j+1)%tree[-1]['vertices'].shape[0]]-tree[-1]['vertices'][j%tree[-1]['vertices'].shape[0]]))
        tree[-1]['r_t'] = np.array(tree[-1]['r_t'])
        tree[-1]['r_n'] = np.zeros((tree[-1]['r_t'].shape[0],2))
        for j in range(0,tree[-1]['r_n'].shape[0]):
            tree[-1]['r_n'][j][0] = -tree[-1]['r_t'][j][1]
            tree[-1]['r_n'][j][1] = tree[-1]['r_t'][j][0]

        # Find the remaining tangents and normals from vertices 0 and 1 to the center
        tree[-1]['r_center_t'] = (tree[-1]['center'][0]-tree[-1]['vertices'][0])/np.linalg.norm(tree[-1]['center'][0]-tree[-1]['vertices'][0])
        tree[-1]['r_center_n'] = np.array([-tree[-1]['r_center_t'][1], tree[-1]['r_center_t'][0]])
        tree[-1]['r_center_t'] = np.vstack((tree[-1]['r_center_t'], (tree[-1]['vertices'][1]-tree[-1]['center'][0])/np.linalg.norm(tree[-1]['vertices'][1]-tree[-1]['center'][0])))
        tree[-1]['r_center_n'] = np.vstack((tree[-1]['r_center_n'], np.array([-tree[-1]['r_center_t'][1][1],tree[-1]['r_center_t'][1][0]])))

        # Compute the dilated polygon and truncate it by the rays emanating from the center
        original_polygon = np.vstack((tree[-1]['vertices'][0], tree[-1]['center'], tree[-1]['vertices'][1:]))
        polygon_tilde = sp.geometry.polygon.orient(Polygon(original_polygon).buffer(varepsilon, join_style=1), 1.0)
        dilation = np.vstack((polygon_tilde.exterior.coords.xy[0], polygon_tilde.exterior.coords.xy[1])).transpose()
        intersect_1 = cvxpolyxhplane(dilation[0:-1], tree[-1]['center'][0], tree[-1]['r_center_n'][0])
        intersect_2 = cvxpolyxhplane(intersect_1, tree[-1]['center'][0], tree[-1]['r_center_n'][1])
        polygon_tilde_vertices = np.vstack((intersect_2,intersect_2[0]))

        # Compute the intersection with the workspace
        final_polygon = sp.geometry.polygon.orient(((Polygon(polygon_tilde_vertices).intersection(Polygon(workspace))).union(Polygon(np.vstack((tree[-1]['center'][0], tree[-1]['vertices'][1:], tree[-1]['vertices'][0], tree[-1]['center'][0]))))).simplify(0.01), 1.0)
        tree[-1]['vertices_tilde'] = np.vstack((final_polygon.exterior.coords.xy[0][0:-1], final_polygon.exterior.coords.xy[1][0:-1])).transpose()

        # Find the tangent and normal vectors for the generated polygonal collar
        tree[-1]['r_tilde_t'] = []
        for j in range(0,tree[-1]['vertices_tilde'].shape[0]):
            tree[-1]['r_tilde_t'].append(np.array(tree[-1]['vertices_tilde'][(j+1)%tree[-1]['vertices_tilde'].shape[0]]-tree[-1]['vertices_tilde'][j%tree[-1]['vertices_tilde'].shape[0]])/np.linalg.norm(tree[-1]['vertices_tilde'][(j+1)%tree[-1]['vertices_tilde'].shape[0]]-tree[-1]['vertices_tilde'][j%tree[-1]['vertices_tilde'].shape[0]]))
        tree[-1]['r_tilde_t'] = np.array(tree[-1]['r_tilde_t'])
        tree[-1]['r_tilde_n'] = np.zeros((tree[-1]['r_tilde_t'].shape[0],2))
        for j in range(0,tree[-1]['r_tilde_n'].shape[0]):
            tree[-1]['r_tilde_n'][j][0] = -tree[-1]['r_tilde_t'][j][1]
            tree[-1]['r_tilde_n'][j][1] = tree[-1]['r_tilde_t'][j][0]
        
        # Finally, compute the augmented inner polygon that includes the center of deformation and update
        tree[-1]['augmented_vertices'] = np.vstack((tree[-1]['vertices'][0], tree[-1]['center'], tree[-1]['vertices'][1:]))
        tree[-1]['r_t'] = np.vstack((tree[-1]['r_center_t'][0], tree[-1]['r_center_t'][1], tree[-1]['r_t'][1:]))
        tree[-1]['r_n'] = np.vstack((tree[-1]['r_center_n'][0], tree[-1]['r_center_n'][1], tree[-1]['r_n'][1:]))
        
        # Add a dummy radius
        tree[-1]['radius'] = 0.0
    else:
        # Compute the convex decomposition tree of the polygon with its dual (adjacency) graph
        tree = polyconvexdecomposition(PolygonVertices, workspace, False)

        # Start with the root and find the center and the radius
        root_coords = tree[-1]['vertices'].transpose()
        tree[-1]['center'] = np.array([[sum(root_coords[0])/len(root_coords[0]), sum(root_coords[1])/len(root_coords[1])]])
        D, closest_point = polydist(tree[-1]['vertices'], tree[-1]['center'])
        tree[-1]['radius'] = 0.8*D[0]

        # Compute the tangent and normal vectors of the root polygon
        tree[-1]['r_t'] = []
        for j in range(0,tree[-1]['vertices'].shape[0]):
            tree[-1]['r_t'].append(np.array(tree[-1]['vertices'][(j+1)%tree[-1]['vertices'].shape[0]]-tree[-1]['vertices'][j%tree[-1]['vertices'].shape[0]])/np.linalg.norm(tree[-1]['vertices'][(j+1)%tree[-1]['vertices'].shape[0]]-tree[-1]['vertices'][j%tree[-1]['vertices'].shape[0]]))
        tree[-1]['r_t'] = np.array(tree[-1]['r_t'])
        tree[-1]['r_n'] = np.zeros((tree[-1]['r_t'].shape[0],2))
        for j in range(0,tree[-1]['r_n'].shape[0]):
            tree[-1]['r_n'][j][0] = -tree[-1]['r_t'][j][1]
            tree[-1]['r_n'][j][1] = tree[-1]['r_t'][j][0]

        # Find the polygonal collar for the root by dilating the polygon by varepsilon
        polygon_tilde = sp.geometry.polygon.orient(Polygon(tree[-1]['vertices']).buffer(varepsilon, join_style=1).intersection(Polygon(workspace)).simplify(0.01), 1.0)
        tree[-1]['vertices_tilde'] = np.vstack((polygon_tilde.exterior.coords.xy[0][0:-1], polygon_tilde.exterior.coords.xy[1][0:-1])).transpose()

        # Find the tangent and normal vectors for the generated polygonal collar
        tree[-1]['r_tilde_t'] = []
        for j in range(0,tree[-1]['vertices_tilde'].shape[0]):
            tree[-1]['r_tilde_t'].append(np.array(tree[-1]['vertices_tilde'][(j+1)%tree[-1]['vertices_tilde'].shape[0]]-tree[-1]['vertices_tilde'][j%tree[-1]['vertices_tilde'].shape[0]])/np.linalg.norm(tree[-1]['vertices_tilde'][(j+1)%tree[-1]['vertices_tilde'].shape[0]]-tree[-1]['vertices_tilde'][j%tree[-1]['vertices_tilde'].shape[0]]))
        tree[-1]['r_tilde_t'] = np.array(tree[-1]['r_tilde_t'])
        tree[-1]['r_tilde_n'] = np.zeros((tree[-1]['r_tilde_t'].shape[0],2))
        for j in range(0,tree[-1]['r_tilde_n'].shape[0]):
            tree[-1]['r_tilde_n'][j][0] = -tree[-1]['r_tilde_t'][j][1]
            tree[-1]['r_tilde_n'][j][1] = tree[-1]['r_tilde_t'][j][0]
        
        # In this case the augmented vertices are the same
        tree[-1]['augmented_vertices'] = tree[-1]['vertices']

    # Identify all the children properties
    for i in range(0,len(tree)-1):
        # Compute the tangent and normal vectors of the child hyperplanes
        # r0 is always the shared edge between the parent and the child, the rest in CCW order
        tree[i]['r_t'] = []
        for j in range(0,tree[i]['vertices'].shape[0]):
            tree[i]['r_t'].append(np.array(tree[i]['vertices'][(j+1)%tree[i]['vertices'].shape[0]]-tree[i]['vertices'][j%tree[i]['vertices'].shape[0]])/np.linalg.norm(tree[i]['vertices'][(j+1)%tree[i]['vertices'].shape[0]]-tree[i]['vertices'][j%tree[i]['vertices'].shape[0]]))
        tree[i]['r_t'] = np.array(tree[i]['r_t'])
        tree[i]['r_n'] = np.zeros((tree[i]['r_t'].shape[0],2))
        for j in range(0,tree[i]['r_n'].shape[0]):
            tree[i]['r_n'][j][0] = -tree[i]['r_t'][j][1]
            tree[i]['r_n'][j][1] = tree[i]['r_t'][j][0]
        
        # Find the median from the point furthest away from the adjacency edge and from that compute the center for the purging transformation
        # To compute the center, first compute the intersection of the parent polygon with the hyperplanes of the child polygon next to the adjacency edge - this defines the admissible region within which you are allowed to search for a center
        dot_product_list = []
        for j in range(0,tree[i]['vertices'].shape[0]):
            dot_product_list.append(np.dot(tree[i]['vertices'][j]-tree[i]['vertices'][0],tree[i]['r_n'][0]))
        inds = np.array(dot_product_list).argsort()
        median_point = 0.5*np.array([[tree[i]['adj_edge'][1][0]+tree[i]['adj_edge'][0][0], tree[i]['adj_edge'][1][1]+tree[i]['adj_edge'][0][1]]])
        median_ray = np.array([[median_point[0][0]-tree[i]['vertices'][inds[-1]][0], median_point[0][1]-tree[i]['vertices'][inds[-1]][1]]])
        median_ray = median_ray/np.linalg.norm(median_ray[0])
        intersect_1_cvx = cvxpolyxhplane(tree[tree[i]['predecessor']]['vertices'], tree[i]['adj_edge'][0], tree[i]['r_n'][-1])
        intersect_2_cvx = cvxpolyxhplane(intersect_1_cvx, tree[i]['adj_edge'][1], tree[i]['r_n'][1])
        intersection_point = polyxray(intersect_2_cvx, median_point[0]+0.05*median_ray[0], median_ray[0]) # offset median point by a little bit to avoid numerical problems
        tree[i]['center'] = np.array([[0.5*median_point[0][0]+0.5*intersection_point[0], 0.5*median_point[0][1]+0.5*intersection_point[1]]])

        # Find the remaining tangents and normals from vertices 0 and 1 to the center
        tree[i]['r_center_t'] = (tree[i]['center'][0]-tree[i]['vertices'][0])/np.linalg.norm(tree[i]['center'][0]-tree[i]['vertices'][0])
        tree[i]['r_center_n'] = np.array([-tree[i]['r_center_t'][1], tree[i]['r_center_t'][0]])
        tree[i]['r_center_t'] = np.vstack((tree[i]['r_center_t'], (tree[i]['vertices'][1]-tree[i]['center'][0])/np.linalg.norm(tree[i]['vertices'][1]-tree[i]['center'][0])))
        tree[i]['r_center_n'] = np.vstack((tree[i]['r_center_n'], np.array([-tree[i]['r_center_t'][1][1],tree[i]['r_center_t'][1][0]])))

        # Compute the dilated polygon and truncate it by the rays emanating from the center
        # Make sure that the intersection of the dilation with the convex pieces succeeding the current convex piece in the transformation does not generate a multipolygon - otherwise reduce the radius of dilation until we have a single piece
        succeeding_polygons = []
        for j in range(i+1,len(tree)):
            succeeding_polygons.append(Polygon(tree[j]['vertices']))
        succeeding_polygons_union = cascaded_union(succeeding_polygons)
        original_polygon = np.vstack((tree[i]['vertices'][0], tree[i]['center'], tree[i]['vertices'][1:]))
        varepsilon_used = varepsilon
        polygon_tilde = sp.geometry.polygon.orient(Polygon(original_polygon).buffer(varepsilon_used, join_style=1).simplify(0.01), 1.0)
        while (polygon_tilde.intersection(succeeding_polygons_union)).geom_type == 'MultiPolygon':
            varepsilon_used = 0.5*varepsilon_used
            polygon_tilde = sp.geometry.polygon.orient(Polygon(original_polygon).buffer(varepsilon_used, join_style=1).simplify(0.01), 1.0)
        dilation = np.vstack((polygon_tilde.exterior.coords.xy[0], polygon_tilde.exterior.coords.xy[1])).transpose()
        intersect_1 = cvxpolyxhplane(dilation[0:-1], tree[i]['center'][0], tree[i]['r_center_n'][0])
        intersect_2 = cvxpolyxhplane(intersect_1, tree[i]['center'][0], tree[i]['r_center_n'][1])
        candidate_polygon_vertices = np.vstack((intersect_2,intersect_2[0]))
        candidate_polygon = Polygon(candidate_polygon_vertices)

        # Check for collisions with all the polygons that will succeed i in the diffeomorphism construction except for its parent
        for j in range(i+1,len(tree)):
            if (j == tree[i]['predecessor']):
                continue
            else:
                polygon_to_test = Polygon(tree[j]['vertices'])
                candidate_polygon = (candidate_polygon.buffer(0)).difference(polygon_to_test.buffer(0))
                # If the difference operation created a multipolygon, keep only the polygon that contains the barycenter of the extended triangle
                if candidate_polygon.geom_type == 'MultiPolygon':
                    point_to_consider = Point((tree[i]['vertices'][0][0]+tree[i]['vertices'][1][0]+tree[i]['center'][0][0])/3.0, (tree[i]['vertices'][0][1]+tree[i]['vertices'][1][1]+tree[i]['center'][0][1])/3.0)
                    for k in range(len(candidate_polygon)):
                        if candidate_polygon[k].contains(point_to_consider):
                            candidate_polygon = candidate_polygon[k]
                            break
        
        # Extract final vertices
        candidate_polygon = sp.geometry.polygon.orient(candidate_polygon.simplify(0.01), 1.0)
        candidate_polygon_vertices = np.vstack((candidate_polygon.exterior.coords.xy[0], candidate_polygon.exterior.coords.xy[1])).transpose()

        # Decompose the polygon into its convex pieces and find the piece that includes the barycenter of the extended triangle
        decomposition = polycvxdecomp(candidate_polygon_vertices.tolist())
        for j in range(len(decomposition)):
            point_to_consider = Point((tree[i]['vertices'][0][0]+tree[i]['vertices'][1][0]+tree[i]['center'][0][0])/3.0, (tree[i]['vertices'][0][1]+tree[i]['vertices'][1][1]+tree[i]['center'][0][1])/3.0)
            polygon_to_consider = Polygon(decomposition[j])
            if polygon_to_consider.buffer(0.01).contains(point_to_consider):
                final_polygon_vertices = np.vstack((polygon_to_consider.exterior.coords.xy[0], polygon_to_consider.exterior.coords.xy[1])).transpose()
                break

        # Generate the outer polygonal collar
        final_polygon = sp.geometry.polygon.orient(Polygon(final_polygon_vertices).intersection(Polygon(workspace)), 1.0)
        tree[i]['vertices_tilde'] = np.vstack((final_polygon.exterior.coords.xy[0][0:-1], final_polygon.exterior.coords.xy[1][0:-1])).transpose()

        # Find the tangent and normal vectors for the generated polygonal collar
        tree[i]['r_tilde_t'] = []
        for j in range(0,tree[i]['vertices_tilde'].shape[0]):
            tree[i]['r_tilde_t'].append(np.array(tree[i]['vertices_tilde'][(j+1)%tree[i]['vertices_tilde'].shape[0]]-tree[i]['vertices_tilde'][j%tree[i]['vertices_tilde'].shape[0]])/np.linalg.norm(tree[i]['vertices_tilde'][(j+1)%tree[i]['vertices_tilde'].shape[0]]-tree[i]['vertices_tilde'][j%tree[i]['vertices_tilde'].shape[0]]))
        tree[i]['r_tilde_t'] = np.array(tree[i]['r_tilde_t'])
        tree[i]['r_tilde_n'] = np.zeros((tree[i]['r_tilde_t'].shape[0],2))
        for j in range(0,tree[i]['r_tilde_n'].shape[0]):
            tree[i]['r_tilde_n'][j][0] = -tree[i]['r_tilde_t'][j][1]
            tree[i]['r_tilde_n'][j][1] = tree[i]['r_tilde_t'][j][0]
        
        # Finally, compute the augmented inner polygon that includes the center of deformation and update
        tree[i]['augmented_vertices'] = np.vstack((tree[i]['vertices'][0], tree[i]['center'], tree[i]['vertices'][1:]))
        tree[i]['r_t'] = np.vstack((tree[i]['r_center_t'][0], tree[i]['r_center_t'][1], tree[i]['r_t'][1:]))
        tree[i]['r_n'] = np.vstack((tree[i]['r_center_n'][0], tree[i]['r_center_n'][1], tree[i]['r_n'][1:]))
        
    return tree


def polygonDiffeoTriangulation(Position, DiffeoTree, DiffeoParams):
    """
    Function that computes h(x) (i.e., the position of point x in the model layer), Dh(x) and the derivatives of Dh, when purging a specific polygon whose dual graph and diffeomorphism properties are known, using the ear clipping method
    
    Input:
        1) Position: Point to consider - 1x2 numpy.array
        2) DiffeoTree: Tree that contains the diffeomorphism properties for a particular polygon
        3) DiffeoParams: Options for the diffeomorphism construction
        
    Output:
        1) map_h: Value of the map
        2) map_hd: Value of the map differential - 1x2 numpy.array
        3) map_hdd: Values of the derivatives of the differential - 8-element numpy array
    """
    # Begin purging process with default values
    if Position.shape == (2,):
        Position = np.array([Position])
    map_h = Position
    map_hd = np.eye(2)
    map_hdd = np.zeros(8)
    
    # Iterate through the polygon triangles
    for i in range(len(DiffeoTree)):
        map_h_new, map_hd_new, map_hdd_new = triangleDiffeo(map_h, DiffeoTree[i], DiffeoParams)

        res1 = map_hd_new[0][0]*map_hdd[0] + map_hd_new[0][1]*map_hdd[4] + map_hd[0][0]*(map_hdd_new[0]*map_hd[0][0] + map_hdd_new[1]*map_hd[1][0]) + map_hd[1][0]*(map_hdd_new[2]*map_hd[0][0] + map_hdd_new[3]*map_hd[1][0])
        res2 = map_hd_new[0][0]*map_hdd[1] + map_hd_new[0][1]*map_hdd[5] + map_hd[0][0]*(map_hdd_new[0]*map_hd[0][1] + map_hdd_new[1]*map_hd[1][1]) + map_hd[1][0]*(map_hdd_new[2]*map_hd[0][1] + map_hdd_new[3]*map_hd[1][1])
        res3 = map_hd_new[0][0]*map_hdd[2] + map_hd_new[0][1]*map_hdd[6] + map_hd[0][1]*(map_hdd_new[0]*map_hd[0][0] + map_hdd_new[1]*map_hd[1][0]) + map_hd[1][1]*(map_hdd_new[2]*map_hd[0][0] + map_hdd_new[3]*map_hd[1][0])
        res4 = map_hd_new[0][0]*map_hdd[3] + map_hd_new[0][1]*map_hdd[7] + map_hd[0][1]*(map_hdd_new[0]*map_hd[0][1] + map_hdd_new[1]*map_hd[1][1]) + map_hd[1][1]*(map_hdd_new[2]*map_hd[0][1] + map_hdd_new[3]*map_hd[1][1])
        res5 = map_hd_new[1][0]*map_hdd[0] + map_hd_new[1][1]*map_hdd[4] + map_hd[0][0]*(map_hdd_new[4]*map_hd[0][0] + map_hdd_new[5]*map_hd[1][0]) + map_hd[1][0]*(map_hdd_new[6]*map_hd[0][0] + map_hdd_new[7]*map_hd[1][0])
        res6 = map_hd_new[1][0]*map_hdd[1] + map_hd_new[1][1]*map_hdd[5] + map_hd[0][0]*(map_hdd_new[4]*map_hd[0][1] + map_hdd_new[5]*map_hd[1][1]) + map_hd[1][0]*(map_hdd_new[6]*map_hd[0][1] + map_hdd_new[7]*map_hd[1][1])
        res7 = map_hd_new[1][0]*map_hdd[2] + map_hd_new[1][1]*map_hdd[6] + map_hd[0][1]*(map_hdd_new[4]*map_hd[0][0] + map_hdd_new[5]*map_hd[1][0]) + map_hd[1][1]*(map_hdd_new[6]*map_hd[0][0] + map_hdd_new[7]*map_hd[1][0])
        res8 = map_hd_new[1][0]*map_hdd[3] + map_hd_new[1][1]*map_hdd[7] + map_hd[0][1]*(map_hdd_new[4]*map_hd[0][1] + map_hdd_new[5]*map_hd[1][1]) + map_hd[1][1]*(map_hdd_new[6]*map_hd[0][1] + map_hdd_new[7]*map_hd[1][1])
        map_hdd[0] = res1
        map_hdd[1] = res2
        map_hdd[2] = res3
        map_hdd[3] = res4
        map_hdd[4] = res5
        map_hdd[5] = res6
        map_hdd[6] = res7
        map_hdd[7] = res8

        map_hd = np.matmul(map_hd_new, map_hd)

        map_h = map_h_new
        
    return map_h, map_hd, map_hdd


def polygonDiffeoConvex(Position, DiffeoTree, DiffeoParams):
    """
    Function that computes h(x) (i.e., the position of point x in the model layer), Dh(x) and the derivatives of Dh, when purging a specific polygon whose dual graph and diffeomorphism properties are known, using the convex decomposition method
    
    Input:
        1) Position: Point to consider - 1x2 numpy.array
        2) DiffeoTree: Tree that contains the diffeomorphism properties for a particular polygon
        3) DiffeoParams: Options for the diffeomorphism construction
        
    Output:
        1) map_h: Value of the map
        2) map_hd: Value of the map differential - 1x2 numpy.array
        3) map_hdd: Values of the derivatives of the differential - 8-element numpy array
    """
    # Begin purging process with default values
    if Position.shape == (2,):
        Position = np.array([Position])
    map_h = Position
    map_hd = np.eye(2)
    map_hdd = np.zeros(8)
    
    # Iterate through the polygon triangles
    for i in range(len(DiffeoTree)):
        map_h_new, map_hd_new, map_hdd_new = polygonDiffeo(map_h, DiffeoTree[i], DiffeoParams)

        res1 = map_hd_new[0][0]*map_hdd[0] + map_hd_new[0][1]*map_hdd[4] + map_hd[0][0]*(map_hdd_new[0]*map_hd[0][0] + map_hdd_new[1]*map_hd[1][0]) + map_hd[1][0]*(map_hdd_new[2]*map_hd[0][0] + map_hdd_new[3]*map_hd[1][0])
        res2 = map_hd_new[0][0]*map_hdd[1] + map_hd_new[0][1]*map_hdd[5] + map_hd[0][0]*(map_hdd_new[0]*map_hd[0][1] + map_hdd_new[1]*map_hd[1][1]) + map_hd[1][0]*(map_hdd_new[2]*map_hd[0][1] + map_hdd_new[3]*map_hd[1][1])
        res3 = map_hd_new[0][0]*map_hdd[2] + map_hd_new[0][1]*map_hdd[6] + map_hd[0][1]*(map_hdd_new[0]*map_hd[0][0] + map_hdd_new[1]*map_hd[1][0]) + map_hd[1][1]*(map_hdd_new[2]*map_hd[0][0] + map_hdd_new[3]*map_hd[1][0])
        res4 = map_hd_new[0][0]*map_hdd[3] + map_hd_new[0][1]*map_hdd[7] + map_hd[0][1]*(map_hdd_new[0]*map_hd[0][1] + map_hdd_new[1]*map_hd[1][1]) + map_hd[1][1]*(map_hdd_new[2]*map_hd[0][1] + map_hdd_new[3]*map_hd[1][1])
        res5 = map_hd_new[1][0]*map_hdd[0] + map_hd_new[1][1]*map_hdd[4] + map_hd[0][0]*(map_hdd_new[4]*map_hd[0][0] + map_hdd_new[5]*map_hd[1][0]) + map_hd[1][0]*(map_hdd_new[6]*map_hd[0][0] + map_hdd_new[7]*map_hd[1][0])
        res6 = map_hd_new[1][0]*map_hdd[1] + map_hd_new[1][1]*map_hdd[5] + map_hd[0][0]*(map_hdd_new[4]*map_hd[0][1] + map_hdd_new[5]*map_hd[1][1]) + map_hd[1][0]*(map_hdd_new[6]*map_hd[0][1] + map_hdd_new[7]*map_hd[1][1])
        res7 = map_hd_new[1][0]*map_hdd[2] + map_hd_new[1][1]*map_hdd[6] + map_hd[0][1]*(map_hdd_new[4]*map_hd[0][0] + map_hdd_new[5]*map_hd[1][0]) + map_hd[1][1]*(map_hdd_new[6]*map_hd[0][0] + map_hdd_new[7]*map_hd[1][0])
        res8 = map_hd_new[1][0]*map_hdd[3] + map_hd_new[1][1]*map_hdd[7] + map_hd[0][1]*(map_hdd_new[4]*map_hd[0][1] + map_hdd_new[5]*map_hd[1][1]) + map_hd[1][1]*(map_hdd_new[6]*map_hd[0][1] + map_hdd_new[7]*map_hd[1][1])
        map_hdd[0] = res1
        map_hdd[1] = res2
        map_hdd[2] = res3
        map_hdd[3] = res4
        map_hdd[4] = res5
        map_hdd[5] = res6
        map_hdd[6] = res7
        map_hdd[7] = res8

        map_hd = np.matmul(map_hd_new, map_hd)

        map_h = map_h_new
        
    return map_h, map_hd, map_hdd


def triangleDiffeo(Position, Triangle, DiffeoParams):
    """
    Function that computes h(x) (i.e., the position of point x in the model layer), Dh(x) and the derivatives of Dh, when purging a specific triangle in space
    
    Input:
        1) Position: Point to consider - 1x2 numpy.array
        2) Triangle: Description of the triangle - Dictionary
        3) DiffeoParams: Options for the diffeomorphism construction
        
    Output:
        1) map_h: Value of the map
        2) map_hd: Value of the map differential - 1x2 numpy.array
        3) map_hdd: Values of the derivatives of the differential - 8-element numpy array
    """
    # Compute the triangle switch and its gradient
    sigma, sigmad, sigmadd = triangleSwitch(Position, Triangle, DiffeoParams)

    # Compute the triangle deforming factor
    nu, nud, nudd = deformingFactor(Position, Triangle)

    # Find the map and its differential
    map_h = sigma*(Triangle['center']+nu*(Position-Triangle['center'])) + (1-sigma)*Position
    map_hd = (nu-1)*np.dot((Position-Triangle['center']).transpose(),sigmad) + sigma*np.dot((Position-Triangle['center']).transpose(),nud) + (1+sigma*(nu-1))*np.eye(2)

    # Find the derivatives of the jacobian
    map_hdd_m0_r0_s0 = 2*sigma*nud[0][0]+2*(nu-1)*sigmad[0][0]+2*(Position[0][0]-Triangle['center'][0][0])*sigmad[0][0]*nud[0][0]+(Position[0][0]-Triangle['center'][0][0])*sigma*nudd[0][0]+(Position[0][0]-Triangle['center'][0][0])*(nu-1)*sigmadd[0][0]
    map_hdd_m0_r0_s1 = sigma*nud[0][1]+(nu-1)*sigmad[0][1]+(Position[0][0]-Triangle['center'][0][0])*sigmad[0][1]*nud[0][0]+(Position[0][0]-Triangle['center'][0][0])*sigma*nudd[0][1]+(Position[0][0]-Triangle['center'][0][0])*sigmad[0][0]*nud[0][1]+(Position[0][0]-Triangle['center'][0][0])*(nu-1)*sigmadd[0][1]
    map_hdd_m0_r1_s0 = sigma*nud[0][1]+(Position[0][0]-Triangle['center'][0][0])*sigmad[0][0]*nud[0][1]+(Position[0][0]-Triangle['center'][0][0])*sigma*nudd[0][1]+(nu-1)*sigmad[0][1]+(Position[0][0]-Triangle['center'][0][0])*sigmad[0][1]*nud[0][0]+(Position[0][0]-Triangle['center'][0][0])*(nu-1)*sigmadd[0][1]
    map_hdd_m0_r1_s1 = 2*(Position[0][0]-Triangle['center'][0][0])*sigmad[0][1]*nud[0][1]+(Position[0][0]-Triangle['center'][0][0])*sigma*nudd[1][1]+(Position[0][0]-Triangle['center'][0][0])*(nu-1)*sigmadd[1][1]
    map_hdd_m1_r0_s0 = 2*(Position[0][1]-Triangle['center'][0][1])*sigmad[0][0]*nud[0][0]+(Position[0][1]-Triangle['center'][0][1])*sigma*nudd[0][0]+(Position[0][1]-Triangle['center'][0][1])*(nu-1)*sigmadd[0][0]
    map_hdd_m1_r0_s1 = sigma*nud[0][0]+(Position[0][1]-Triangle['center'][0][1])*sigmad[0][1]*nud[0][0]+(Position[0][1]-Triangle['center'][0][1])*sigma*nudd[0][1]+(nu-1)*sigmad[0][0]+(Position[0][1]-Triangle['center'][0][1])*sigmad[0][0]*nud[0][1]+(Position[0][1]-Triangle['center'][0][1])*(nu-1)*sigmadd[0][1]
    map_hdd_m1_r1_s0 = sigma*nud[0][0]+(nu-1)*sigmad[0][0]+(Position[0][1]-Triangle['center'][0][1])*sigmad[0][0]*nud[0][1]+(Position[0][1]-Triangle['center'][0][1])*sigma*nudd[0][1]+(Position[0][1]-Triangle['center'][0][1])*sigmad[0][1]*nud[0][0]+(Position[0][1]-Triangle['center'][0][1])*(nu-1)*sigmadd[0][1]
    map_hdd_m1_r1_s1 = 2*sigma*nud[0][1]+2*(nu-1)*sigmad[0][1]+2*(Position[0][1]-Triangle['center'][0][1])*sigmad[0][1]*nud[0][1]+(Position[0][1]-Triangle['center'][0][1])*sigma*nudd[1][1]+(Position[0][1]-Triangle['center'][0][1])*(nu-1)*sigmadd[1][1]
    
    map_hdd = np.array([map_hdd_m0_r0_s0, map_hdd_m0_r0_s1, map_hdd_m0_r1_s0, map_hdd_m0_r1_s1, map_hdd_m1_r0_s0, map_hdd_m1_r0_s1, map_hdd_m1_r1_s0, map_hdd_m1_r1_s1])

    return map_h, map_hd, map_hdd


def polygonDiffeo(Position, PolygonUsed, DiffeoParams):
    """
    Function that computes h(x) (i.e., the position of point x in the model layer), Dh(x) and the derivatives of Dh, when purging a specific polygon in space
    
    Input:
        1) Position: Point to consider - 1x2 numpy.array
        2) PolygonUsed: Description of the polygon - Dictionary
        3) DiffeoParams: Options for the diffeomorphism construction
        
    Output:
        1) map_h: Value of the map
        2) map_hd: Value of the map differential - 1x2 numpy.array
        3) map_hdd: Values of the derivatives of the differential - 8-element numpy array
    """
    # Compute the polygon switch and its gradient
    sigma, sigmad, sigmadd = polygonSwitch(Position, PolygonUsed, DiffeoParams)

    # Compute the polygon deforming factor
    nu, nud, nudd = deformingFactor(Position, PolygonUsed)

    # Find the map and its differential
    map_h = sigma*(PolygonUsed['center']+nu*(Position-PolygonUsed['center'])) + (1-sigma)*Position
    map_hd = (nu-1)*np.dot((Position-PolygonUsed['center']).transpose(),sigmad) + sigma*np.dot((Position-PolygonUsed['center']).transpose(),nud) + (1+sigma*(nu-1))*np.eye(2)

    # Find the derivatives of the jacobian
    map_hdd_m0_r0_s0 = 2*sigma*nud[0][0]+2*(nu-1)*sigmad[0][0]+2*(Position[0][0]-PolygonUsed['center'][0][0])*sigmad[0][0]*nud[0][0]+(Position[0][0]-PolygonUsed['center'][0][0])*sigma*nudd[0][0]+(Position[0][0]-PolygonUsed['center'][0][0])*(nu-1)*sigmadd[0][0]
    map_hdd_m0_r0_s1 = sigma*nud[0][1]+(nu-1)*sigmad[0][1]+(Position[0][0]-PolygonUsed['center'][0][0])*sigmad[0][1]*nud[0][0]+(Position[0][0]-PolygonUsed['center'][0][0])*sigma*nudd[0][1]+(Position[0][0]-PolygonUsed['center'][0][0])*sigmad[0][0]*nud[0][1]+(Position[0][0]-PolygonUsed['center'][0][0])*(nu-1)*sigmadd[0][1]
    map_hdd_m0_r1_s0 = sigma*nud[0][1]+(Position[0][0]-PolygonUsed['center'][0][0])*sigmad[0][0]*nud[0][1]+(Position[0][0]-PolygonUsed['center'][0][0])*sigma*nudd[0][1]+(nu-1)*sigmad[0][1]+(Position[0][0]-PolygonUsed['center'][0][0])*sigmad[0][1]*nud[0][0]+(Position[0][0]-PolygonUsed['center'][0][0])*(nu-1)*sigmadd[0][1]
    map_hdd_m0_r1_s1 = 2*(Position[0][0]-PolygonUsed['center'][0][0])*sigmad[0][1]*nud[0][1]+(Position[0][0]-PolygonUsed['center'][0][0])*sigma*nudd[1][1]+(Position[0][0]-PolygonUsed['center'][0][0])*(nu-1)*sigmadd[1][1]
    map_hdd_m1_r0_s0 = 2*(Position[0][1]-PolygonUsed['center'][0][1])*sigmad[0][0]*nud[0][0]+(Position[0][1]-PolygonUsed['center'][0][1])*sigma*nudd[0][0]+(Position[0][1]-PolygonUsed['center'][0][1])*(nu-1)*sigmadd[0][0]
    map_hdd_m1_r0_s1 = sigma*nud[0][0]+(Position[0][1]-PolygonUsed['center'][0][1])*sigmad[0][1]*nud[0][0]+(Position[0][1]-PolygonUsed['center'][0][1])*sigma*nudd[0][1]+(nu-1)*sigmad[0][0]+(Position[0][1]-PolygonUsed['center'][0][1])*sigmad[0][0]*nud[0][1]+(Position[0][1]-PolygonUsed['center'][0][1])*(nu-1)*sigmadd[0][1]
    map_hdd_m1_r1_s0 = sigma*nud[0][0]+(nu-1)*sigmad[0][0]+(Position[0][1]-PolygonUsed['center'][0][1])*sigmad[0][0]*nud[0][1]+(Position[0][1]-PolygonUsed['center'][0][1])*sigma*nudd[0][1]+(Position[0][1]-PolygonUsed['center'][0][1])*sigmad[0][1]*nud[0][0]+(Position[0][1]-PolygonUsed['center'][0][1])*(nu-1)*sigmadd[0][1]
    map_hdd_m1_r1_s1 = 2*sigma*nud[0][1]+2*(nu-1)*sigmad[0][1]+2*(Position[0][1]-PolygonUsed['center'][0][1])*sigmad[0][1]*nud[0][1]+(Position[0][1]-PolygonUsed['center'][0][1])*sigma*nudd[1][1]+(Position[0][1]-PolygonUsed['center'][0][1])*(nu-1)*sigmadd[1][1]
    
    map_hdd = np.array([map_hdd_m0_r0_s0, map_hdd_m0_r0_s1, map_hdd_m0_r1_s0, map_hdd_m0_r1_s1, map_hdd_m1_r0_s0, map_hdd_m1_r0_s1, map_hdd_m1_r1_s0, map_hdd_m1_r1_s1])

    return map_h, map_hd, map_hdd


def triangleSwitch(Position, Triangle, DiffeoParams):
    """
    Function that computes the overall switch value, its gradient and hessian for a point x outside a triangle
    
    Input:
        1) Position: Point to consider - 1x2 numpy.array
        2) Triangle: Description of the triangle - Dictionary
        3) DiffeoParams: Options for the diffeomorphism construction
        
    Output:
        1) sigma: Value of the switch
        2) sigmad: Value of the switch gradient - 1x2 numpy.array
        3) sigmadd: Value of the switch hessian - 1x2 numpy.array
    """
    # Distinguish whether the triangle to consider is the root or some child
    # Find the separate switch values, gradients and hessians
    sigma_beta, sigma_betad, sigma_betadd = betaSwitchTriangle(Position, Triangle, DiffeoParams)
    sigma_gamma, sigma_gammad, sigma_gammadd = gammaSwitchTriangle(Position, Triangle, DiffeoParams)

    # Find the overall switch value, gradient and hessian
    if (sigma_beta == 1.0) and (sigma_gamma == 0.0):
        sigma = 1.0
        sigmad = np.array([[0.,0.]])
        sigmadd = np.zeros((2,2))
    else:
        nom = sigma_gamma*sigma_beta
        denom = sigma_gamma*sigma_beta + (1-sigma_beta)
        sigma = nom/denom

        nomd = sigma_gamma*sigma_betad + sigma_beta*sigma_gammad
        denomd = sigma_gamma*sigma_betad + sigma_beta*sigma_gammad - sigma_betad
        sigmad = (1/denom)*nomd - (nom/denom**2)*denomd

        nomdd = sigma_gamma*sigma_betadd + np.dot(sigma_betad.transpose(), sigma_gammad) + np.dot(sigma_gammad.transpose(), sigma_betad) + sigma_beta*sigma_gammadd
        denomdd = sigma_gamma*sigma_betadd + np.dot(sigma_betad.transpose(), sigma_gammad) + np.dot(sigma_gammad.transpose(), sigma_betad) + sigma_beta*sigma_gammadd - sigma_betadd
        sigmadd = (1/denom)*nomdd - (1/denom**2)*(np.dot(nomd.transpose(),denomd)+np.dot(denomd.transpose(),nomd)) + 2*(nom/(denom**3))*np.dot(denomd.transpose(),denomd) - (nom/(denom**2))*denomdd
    
    return sigma, sigmad, sigmadd


def polygonSwitch(Position, PolygonUsed, DiffeoParams):
    """
    Function that computes the overall switch value, its gradient and hessian for a point x outside a triangle
    
    Input:
        1) Position: Point to consider - 1x2 numpy.array
        2) PolygonUsed: Description of the polygon - Dictionary
        3) DiffeoParams: Options for the diffeomorphism construction
        
    Output:
        1) sigma: Value of the switch
        2) sigmad: Value of the switch gradient - 1x2 numpy.array
        3) sigmadd: Value of the switch hessian - 1x2 numpy.array
    """
    # Find the separate switch values, gradients and hessians
    sigma_beta, sigma_betad, sigma_betadd = betaSwitchPolygon(Position, PolygonUsed, DiffeoParams)
    sigma_gamma, sigma_gammad, sigma_gammadd = gammaSwitchPolygon(Position, PolygonUsed, DiffeoParams)

    # Find the overall switch value, gradient and hessian
    if (sigma_beta == 1.0) and (sigma_gamma == 0.0):
        sigma = 1.0
        sigmad = np.array([[0.,0.]])
        sigmadd = np.zeros((2,2))
    else:
        nom = sigma_gamma*sigma_beta
        denom = sigma_gamma*sigma_beta + (1-sigma_beta)
        sigma = nom/denom

        nomd = sigma_gamma*sigma_betad + sigma_beta*sigma_gammad
        denomd = sigma_gamma*sigma_betad + sigma_beta*sigma_gammad - sigma_betad
        sigmad = (1/denom)*nomd - (nom/denom**2)*denomd

        nomdd = sigma_gamma*sigma_betadd + np.dot(sigma_betad.transpose(), sigma_gammad) + np.dot(sigma_gammad.transpose(), sigma_betad) + sigma_beta*sigma_gammadd
        denomdd = sigma_gamma*sigma_betadd + np.dot(sigma_betad.transpose(), sigma_gammad) + np.dot(sigma_gammad.transpose(), sigma_betad) + sigma_beta*sigma_gammadd - sigma_betadd
        sigmadd = (1/denom)*nomdd - (1/denom**2)*(np.dot(nomd.transpose(),denomd)+np.dot(denomd.transpose(),nomd)) + 2*(nom/(denom**3))*np.dot(denomd.transpose(),denomd) - (nom/(denom**2))*denomdd
    
    return sigma, sigmad, sigmadd


def deformingFactor(Position, PolygonUsed):
    """
    Function that computes the value, gradient and hessian of the deforming factor for a point x outside a triangle
    
    Input:
        1) Position: Point to consider - 1x2 numpy.array
        2) PolygonUsed: Description of the polygon - Dictionary
        
    Output:
        1) nu: Value of the deforming factor
        2) nud: Value of the deforming factor gradient - 1x2 numpy.array
        3) nudd: Value of the deforming factor hessian - 2x2 numpy.array
    """
    # Distinguish whether the polygon to consider is the root or some child
    if (PolygonUsed['depth'] == 0) and (not PolygonUsed['adj_edge'].any()):
        nu = PolygonUsed['radius']/np.linalg.norm(Position[0]-PolygonUsed['center'][0])
        nud = -(PolygonUsed['radius']/(np.linalg.norm(Position[0]-PolygonUsed['center'][0])**3))*np.array([Position[0]-PolygonUsed['center'][0]])
        nudd = ((3*PolygonUsed['radius'])/(np.linalg.norm(Position[0]-PolygonUsed['center'][0])**5))*np.dot(np.array([Position[0]-PolygonUsed['center'][0]]).transpose(),np.array([Position[0]-PolygonUsed['center'][0]])) - (PolygonUsed['radius']/(np.linalg.norm(Position[0]-PolygonUsed['center'][0])**3))*np.eye(2)
    else:
        # First compute the normal for the adjacency edge
        shared_tangent = np.array(PolygonUsed['adj_edge'][1]-PolygonUsed['adj_edge'][0])/np.linalg.norm(PolygonUsed['adj_edge'][1]-PolygonUsed['adj_edge'][0])
        shared_normal = np.array([-shared_tangent[1],shared_tangent[0]])

        nu = np.dot(PolygonUsed['vertices'][0]-PolygonUsed['center'][0],shared_normal)/np.dot(Position[0]-PolygonUsed['center'][0],shared_normal)
        nud = -((np.dot(PolygonUsed['vertices'][0]-PolygonUsed['center'][0],shared_normal))/(np.dot(Position[0]-PolygonUsed['center'][0],shared_normal))**2)*np.array([shared_normal])
        nudd = 2*((np.dot(PolygonUsed['vertices'][0]-PolygonUsed['center'][0],shared_normal))/(np.dot(Position[0]-PolygonUsed['center'][0],shared_normal))**3)*np.dot(np.array([shared_normal]).transpose(),np.array([shared_normal]))
        
    return nu, nud, nudd


def betaSwitchTriangle(Position, Triangle, DiffeoParams):
    """
    Function that computes the value, gradient and hessian of the beta-switch for a point x outside a triangle
    
    Input:
        1) Position: Point to consider - 1x2 numpy.array
        2) Triangle: Description of the triangle - Dictionary
        3) DiffeoParams: Options for the diffeomorphism construction
        
    Output:
        1) sigma: Value of the switch
        2) sigmad: Value of the switch gradient - 1x2 numpy.array
        3) sigmadd: Value of the switch hessian - 2x2 numpy.array
    """
    # Unwrap parameters
    mu_1 = DiffeoParams['mu_1']
    epsilon = DiffeoParams['epsilon']

    # Compute the value of beta and its gradient and hessian
    beta, betad, betadd = outsideImplicitTriangle(Position, Triangle, DiffeoParams)

    # Compute the value of the switch
    if beta >= epsilon:
        sigma = 0.
        sigmad = np.array([[0.,0.]])
        sigmadd = np.zeros((2,2))
    else:
        sigma = np.exp(-mu_1/(epsilon-beta))/np.exp(-mu_1/(epsilon))
        sigmad = -mu_1*(sigma/((epsilon-beta)**2))*betad
        sigmadd = (mu_1**2*(sigma/((epsilon-beta)**4))-2*mu_1*(sigma/((epsilon-beta)**3)))*np.dot(betad.transpose(),betad) - mu_1*(sigma/((epsilon-beta)**2))*betadd
    
    return sigma, sigmad, sigmadd


def betaSwitchPolygon(Position, PolygonUsed, DiffeoParams):
    """
    Function that computes the value, gradient and hessian of the beta-switch for a point x outside a polygon
    
    Input:
        1) Position: Point to consider - 1x2 numpy.array
        2) PolygonUsed: Description of the polygon - Dictionary
        3) DiffeoParams: Options for the diffeomorphism construction
        
    Output:
        1) sigma: Value of the switch
        2) sigmad: Value of the switch gradient - 1x2 numpy.array
        3) sigmadd: Value of the switch hessian - 2x2 numpy.array
    """
    # Unwrap parameters
    mu_1 = DiffeoParams['mu_1']
    epsilon = DiffeoParams['epsilon']

    # Compute the value of beta and its gradient and hessian
    beta, betad, betadd = outsideImplicitPolygon(Position, PolygonUsed, DiffeoParams)

    # Compute the value of the switch
    if beta >= epsilon:
        sigma = 0.
        sigmad = np.array([[0.,0.]])
        sigmadd = np.zeros((2,2))
    else:
        sigma = np.exp(-mu_1/(epsilon-beta))/np.exp(-mu_1/(epsilon))
        sigmad = -mu_1*(sigma/((epsilon-beta)**2))*betad
        sigmadd = (mu_1**2*(sigma/((epsilon-beta)**4))-2*mu_1*(sigma/((epsilon-beta)**3)))*np.dot(betad.transpose(),betad) - mu_1*(sigma/((epsilon-beta)**2))*betadd
    
    return sigma, sigmad, sigmadd


def gammaSwitchTriangle(Position, Triangle, DiffeoParams):
    """
    Function that computes the value, gradient and hessian of the gamma-switch for a point x outside a triangle
    
    Input:
        1) Position: Point to consider - 1x2 numpy.array
        2) Triangle: Description of the triangle - Dictionary
        3) DiffeoParams: Options for the diffeomorphism construction
        
    Output:
        1) sigma: Value of the switch
        2) sigmad: Value of the switch gradient - 1x2 numpy.array
        3) sigmadd: Value of the switch hessian - 2x2 numpy.array
    """
    # Unwrap parameters
    mu_2 = DiffeoParams['mu_2']

    # Compute the value of gamma and its gradient
    gamma, gammad, gammadd = insideImplicitTriangle(Position, Triangle, DiffeoParams)

    # Compute the value of the switch
    if gamma <= 0.:
        sigma = 0.
        sigmad = np.array([[0.,0.]])
        sigmadd = np.zeros((2,2))
    else:
        # Compute the value of alpha and its gradient and hessian
        nom = gamma
        denom = np.linalg.norm(Position[0]-Triangle['center'][0])
        alpha = nom/denom

        nomd = gammad
        denomd = np.array([Position[0]-Triangle['center'][0]])/np.linalg.norm(Position[0]-Triangle['center'][0])
        alphad = (1/denom)*nomd - (nom/denom**2)*denomd

        nomdd = gammadd
        denomdd = (1/np.linalg.norm(Position[0]-Triangle['center'][0]))*np.eye(2) - (1/np.linalg.norm(Position[0]-Triangle['center'][0])**3)*np.dot(np.array([Position[0]-Triangle['center'][0]]).transpose(),np.array([Position[0]-Triangle['center'][0]]))
        alphadd = (1/denom)*nomdd - (1/denom**2)*(np.dot(nomd.transpose(),denomd)+np.dot(denomd.transpose(),nomd)) + 2*(nom/(denom**3))*np.dot(denomd.transpose(),denomd) - (nom/(denom**2))*denomdd

        sigma = np.exp(-mu_2/alpha)
        sigmad = mu_2*(sigma/(alpha**2))*alphad
        sigmadd = (mu_2**2*(sigma/(alpha**4))-2*mu_2*(sigma/(alpha**3)))*np.dot(alphad.transpose(),alphad) + mu_2*(sigma/(alpha**2))*alphadd
    
    return sigma, sigmad, sigmadd


def gammaSwitchPolygon(Position, PolygonUsed, DiffeoParams):
    """
    Function that computes the value, gradient and hessian of the gamma-switch for a point x outside a polygon
    
    Input:
        1) Position: Point to consider - 1x2 numpy.array
        2) PolygonUsed: Description of the polygon - Dictionary
        3) DiffeoParams: Options for the diffeomorphism construction
        
    Output:
        1) sigma: Value of the switch
        2) sigmad: Value of the switch gradient - 1x2 numpy.array
        3) sigmadd: Value of the switch hessian - 2x2 numpy.array
    """
    # Unwrap parameters
    mu_2 = DiffeoParams['mu_2']

    # Compute the value of gamma and its gradient
    gamma, gammad, gammadd = insideImplicitPolygon(Position, PolygonUsed, DiffeoParams)

    # Compute the value of the switch
    if gamma <= 0.:
        sigma = 0.
        sigmad = np.array([[0.,0.]])
        sigmadd = np.zeros((2,2))
    else:
        # Compute the value of alpha and its gradient and hessian
        nom = gamma
        denom = np.linalg.norm(Position[0]-PolygonUsed['center'][0])
        alpha = nom/denom

        nomd = gammad
        denomd = np.array([Position[0]-PolygonUsed['center'][0]])/np.linalg.norm(Position[0]-PolygonUsed['center'][0])
        alphad = (1/denom)*nomd - (nom/denom**2)*denomd

        nomdd = gammadd
        denomdd = (1/np.linalg.norm(Position[0]-PolygonUsed['center'][0]))*np.eye(2) - (1/np.linalg.norm(Position[0]-PolygonUsed['center'][0])**3)*np.dot(np.array([Position[0]-PolygonUsed['center'][0]]).transpose(),np.array([Position[0]-PolygonUsed['center'][0]]))
        alphadd = (1/denom)*nomdd - (1/denom**2)*(np.dot(nomd.transpose(),denomd)+np.dot(denomd.transpose(),nomd)) + 2*(nom/(denom**3))*np.dot(denomd.transpose(),denomd) - (nom/(denom**2))*denomdd

        sigma = np.exp(-mu_2/alpha)
        sigmad = mu_2*(sigma/(alpha**2))*alphad
        sigmadd = (mu_2**2*(sigma/(alpha**4))-2*mu_2*(sigma/(alpha**3)))*np.dot(alphad.transpose(),alphad) + mu_2*(sigma/(alpha**2))*alphadd
    
    return sigma, sigmad, sigmadd


def outsideImplicitTriangle(Position, Triangle, DiffeoParams):
    """
    Function that computes beta(x) (i.e., the R-function) for a point x outside a triangle, and its gradient and hessian
    
    Input:
        1) Position: Point to consider - 1x2 numpy.array
        2) Triangle: Description of the triangle - Dictionary
        3) DiffeoParams: Options for the diffeomorphism construction
        
    Output:
        1) beta: beta(x)
        2) betad: Gradient of beta(x)
        3) betadd: Hessian of beta(x)
    """
    # Unwrap parameters
    p = DiffeoParams['p']

    # Distinguish between the root triangle and all the other triangles
    if (Triangle['depth'] == 0) and (not Triangle['adj_edge'].any()):
        # Find the hyperplane values
        hyperplane_1 = np.dot(Position[0]-Triangle['vertices'][2],Triangle['r_n'][1])
        hyperplane_2 = np.dot(Position[0]-Triangle['vertices'][2],Triangle['r_n'][2])
        hyperplane_3 = np.dot(Position[0]-Triangle['vertices'][0],Triangle['r_n'][0])

        # Compute the R-function
        hyperplane_12 = hyperplane_1 + hyperplane_2 - (hyperplane_1**p+hyperplane_2**p)**(1/p)
        hyperplane_123 = hyperplane_12 + hyperplane_3 - (hyperplane_12**p+hyperplane_3**p)**(1/p)
        beta = -hyperplane_123

        # Compute the gradients
        hyperplane_1d = np.array([Triangle['r_n'][1]])
        hyperplane_2d = np.array([Triangle['r_n'][2]])
        hyperplane_3d = np.array([Triangle['r_n'][0]])
        hyperplane_12d = (1-((hyperplane_1**(p-1))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p))))*hyperplane_1d + (1-((hyperplane_2**(p-1))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p))))*hyperplane_2d
        hyperplane_123d = (1-((hyperplane_12**(p-1))/((hyperplane_12**p+hyperplane_3**p)**((p-1)/p))))*hyperplane_12d + (1-((hyperplane_3**(p-1))/((hyperplane_12**p+hyperplane_3**p)**((p-1)/p))))*hyperplane_3d
        betad = -hyperplane_123d

        # Compute the hessian
        hyperplane_1dd = np.zeros((2,2))
        hyperplane_2dd = np.zeros((2,2))
        hyperplane_3dd = np.zeros((2,2))
        hyperplane_12dd = (-(p-1)*((hyperplane_1**(p-2))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p)))+(p-1)*((hyperplane_1**(2*p-2))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p+1))))*np.dot(hyperplane_1d.transpose(),hyperplane_1d) + ((p-1)*((hyperplane_1**(p-1)*hyperplane_2**(p-1))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p+1))))*np.dot(hyperplane_1d.transpose(),hyperplane_2d) + (1-((hyperplane_1**(p-1))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p))))*hyperplane_1dd + (-(p-1)*((hyperplane_2**(p-2))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p)))+(p-1)*((hyperplane_2**(2*p-2))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p+1))))*np.dot(hyperplane_2d.transpose(),hyperplane_2d) + ((p-1)*((hyperplane_1**(p-1)*hyperplane_2**(p-1))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p+1))))*np.dot(hyperplane_2d.transpose(),hyperplane_1d) + (1-((hyperplane_2**(p-1))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p))))*hyperplane_2dd
        hyperplane_123dd = (-(p-1)*((hyperplane_12**(p-2))/((hyperplane_12**p+hyperplane_3**p)**((p-1)/p)))+(p-1)*((hyperplane_12**(2*p-2))/((hyperplane_12**p+hyperplane_3**p)**((p-1)/p+1))))*np.dot(hyperplane_12d.transpose(),hyperplane_12d) + ((p-1)*((hyperplane_12**(p-1)*hyperplane_3**(p-1))/((hyperplane_12**p+hyperplane_3**p)**((p-1)/p+1))))*np.dot(hyperplane_12d.transpose(),hyperplane_3d) + (1-((hyperplane_12**(p-1))/((hyperplane_12**p+hyperplane_3**p)**((p-1)/p))))*hyperplane_12dd + (-(p-1)*((hyperplane_3**(p-2))/((hyperplane_12**p+hyperplane_3**p)**((p-1)/p)))+(p-1)*((hyperplane_3**(2*p-2))/((hyperplane_12**p+hyperplane_3**p)**((p-1)/p+1))))*np.dot(hyperplane_3d.transpose(),hyperplane_3d) + ((p-1)*((hyperplane_12**(p-1)*hyperplane_3**(p-1))/((hyperplane_12**p+hyperplane_3**p)**((p-1)/p+1))))*np.dot(hyperplane_3d.transpose(),hyperplane_12d) + (1-((hyperplane_3**(p-1))/((hyperplane_12**p+hyperplane_3**p)**((p-1)/p))))*hyperplane_3dd
        betadd = -hyperplane_123dd
    else:
        # Find the hyperplane values
        hyperplane_1 = np.dot(Position[0]-Triangle['vertices'][2],Triangle['r_n'][1])
        hyperplane_2 = np.dot(Position[0]-Triangle['vertices'][2],Triangle['r_n'][2])
        hyperplane_3 = np.dot(Position[0]-Triangle['center'][0],Triangle['r_center_n'][0])
        hyperplane_4 = np.dot(Position[0]-Triangle['center'][0],Triangle['r_center_n'][1])

        # Compute the R-function
        hyperplane_12 = hyperplane_1 + hyperplane_2 - (hyperplane_1**p+hyperplane_2**p)**(1/p)
        hyperplane_34 = hyperplane_3 + hyperplane_4 - (hyperplane_3**p+hyperplane_4**p)**(1/p)
        hyperplane_1234 = hyperplane_12 + hyperplane_34 - (hyperplane_12**p+hyperplane_34**p)**(1/p)
        beta = -hyperplane_1234

        # Compute the gradients
        hyperplane_1d = np.array([Triangle['r_n'][1]])
        hyperplane_2d = np.array([Triangle['r_n'][2]])
        hyperplane_3d = np.array([Triangle['r_center_n'][0]])
        hyperplane_4d = np.array([Triangle['r_center_n'][1]])
        hyperplane_12d = (1-((hyperplane_1**(p-1))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p))))*hyperplane_1d + (1-((hyperplane_2**(p-1))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p))))*hyperplane_2d
        hyperplane_34d = (1-((hyperplane_3**(p-1))/((hyperplane_3**p+hyperplane_4**p)**((p-1)/p))))*hyperplane_3d + (1-((hyperplane_4**(p-1))/((hyperplane_3**p+hyperplane_4**p)**((p-1)/p))))*hyperplane_4d
        hyperplane_1234d = (1-((hyperplane_12**(p-1))/((hyperplane_12**p+hyperplane_34**p)**((p-1)/p))))*hyperplane_12d + (1-((hyperplane_34**(p-1))/((hyperplane_12**p+hyperplane_34**p)**((p-1)/p))))*hyperplane_34d
        betad = -hyperplane_1234d

        # Compute the hessian
        hyperplane_1dd = np.zeros((2,2))
        hyperplane_2dd = np.zeros((2,2))
        hyperplane_3dd = np.zeros((2,2))
        hyperplane_4dd = np.zeros((2,2))
        hyperplane_12dd = (-(p-1)*((hyperplane_1**(p-2))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p)))+(p-1)*((hyperplane_1**(2*p-2))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p+1))))*np.dot(hyperplane_1d.transpose(),hyperplane_1d) + ((p-1)*((hyperplane_1**(p-1)*hyperplane_2**(p-1))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p+1))))*np.dot(hyperplane_1d.transpose(),hyperplane_2d) + (1-((hyperplane_1**(p-1))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p))))*hyperplane_1dd + (-(p-1)*((hyperplane_2**(p-2))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p)))+(p-1)*((hyperplane_2**(2*p-2))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p+1))))*np.dot(hyperplane_2d.transpose(),hyperplane_2d) + ((p-1)*((hyperplane_1**(p-1)*hyperplane_2**(p-1))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p+1))))*np.dot(hyperplane_2d.transpose(),hyperplane_1d) + (1-((hyperplane_2**(p-1))/((hyperplane_1**p+hyperplane_2**p)**((p-1)/p))))*hyperplane_2dd
        hyperplane_34dd = (-(p-1)*((hyperplane_3**(p-2))/((hyperplane_3**p+hyperplane_4**p)**((p-1)/p)))+(p-1)*((hyperplane_3**(2*p-2))/((hyperplane_3**p+hyperplane_4**p)**((p-1)/p+1))))*np.dot(hyperplane_3d.transpose(),hyperplane_3d) + ((p-1)*((hyperplane_3**(p-1)*hyperplane_4**(p-1))/((hyperplane_3**p+hyperplane_4**p)**((p-1)/p+1))))*np.dot(hyperplane_3d.transpose(),hyperplane_4d) + (1-((hyperplane_3**(p-1))/((hyperplane_3**p+hyperplane_4**p)**((p-1)/p))))*hyperplane_3dd + (-(p-1)*((hyperplane_4**(p-2))/((hyperplane_3**p+hyperplane_4**p)**((p-1)/p)))+(p-1)*((hyperplane_4**(2*p-2))/((hyperplane_3**p+hyperplane_4**p)**((p-1)/p+1))))*np.dot(hyperplane_4d.transpose(),hyperplane_4d) + ((p-1)*((hyperplane_3**(p-1)*hyperplane_4**(p-1))/((hyperplane_3**p+hyperplane_4**p)**((p-1)/p+1))))*np.dot(hyperplane_4d.transpose(),hyperplane_3d) + (1-((hyperplane_4**(p-1))/((hyperplane_3**p+hyperplane_4**p)**((p-1)/p))))*hyperplane_4dd
        hyperplane_1234dd = (-(p-1)*((hyperplane_12**(p-2))/((hyperplane_12**p+hyperplane_34**p)**((p-1)/p)))+(p-1)*((hyperplane_12**(2*p-2))/((hyperplane_12**p+hyperplane_34**p)**((p-1)/p+1))))*np.dot(hyperplane_12d.transpose(),hyperplane_12d) + ((p-1)*((hyperplane_12**(p-1)*hyperplane_34**(p-1))/((hyperplane_12**p+hyperplane_34**p)**((p-1)/p+1))))*np.dot(hyperplane_12d.transpose(),hyperplane_34d) + (1-((hyperplane_12**(p-1))/((hyperplane_12**p+hyperplane_34**p)**((p-1)/p))))*hyperplane_12dd + (-(p-1)*((hyperplane_34**(p-2))/((hyperplane_12**p+hyperplane_34**p)**((p-1)/p)))+(p-1)*((hyperplane_34**(2*p-2))/((hyperplane_12**p+hyperplane_34**p)**((p-1)/p+1))))*np.dot(hyperplane_34d.transpose(),hyperplane_34d) + ((p-1)*((hyperplane_12**(p-1)*hyperplane_34**(p-1))/((hyperplane_12**p+hyperplane_34**p)**((p-1)/p+1))))*np.dot(hyperplane_34d.transpose(),hyperplane_12d) + (1-((hyperplane_34**(p-1))/((hyperplane_12**p+hyperplane_34**p)**((p-1)/p))))*hyperplane_34dd
        betadd = -hyperplane_1234dd

    return beta, betad, betadd


def outsideImplicitPolygon(Position, PolygonUsed, DiffeoParams):
    """
    Function that computes beta(x) (i.e., the R-function) for a point x outside a polygon, and its gradient and hessian
    
    Input:
        1) Position: Point to consider - 1x2 numpy.array
        2) PolygonUsed: Description of the polygon - Dictionary
        3) DiffeoParams: Options for the diffeomorphism construction
        
    Output:
        1) beta: beta(x)
        2) betad: Gradient of beta(x)
        3) betadd: Hessian of beta(x)
    """
    # Unwrap parameters
    p = DiffeoParams['p']

    # Compute hyperplane functions
    hyperplane = np.zeros(PolygonUsed['augmented_vertices'].shape[0])
    for i in range(0, PolygonUsed['augmented_vertices'].shape[0]):
        hyperplane[i] = np.dot(Position[0]-PolygonUsed['augmented_vertices'][i],PolygonUsed['r_n'][i])

    # Compute the R-function and its gradient and hessian
    betadd = (-(p-1)*((hyperplane[0]**(p-2))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p)))+(p-1)*((hyperplane[0]**(2*p-2))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p+1))))*np.dot(np.array([PolygonUsed['r_n'][0]]).transpose(),np.array([PolygonUsed['r_n'][0]])) + ((p-1)*((hyperplane[0]**(p-1)*hyperplane[1]**(p-1))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p+1))))*np.dot(np.array([PolygonUsed['r_n'][0]]).transpose(),np.array([PolygonUsed['r_n'][1]])) + (1-((hyperplane[0]**(p-1))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p))))*np.zeros((2,2)) + (-(p-1)*((hyperplane[1]**(p-2))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p)))+(p-1)*((hyperplane[1]**(2*p-2))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p+1))))*np.dot(np.array([PolygonUsed['r_n'][1]]).transpose(),np.array([PolygonUsed['r_n'][1]])) + ((p-1)*((hyperplane[0]**(p-1)*hyperplane[1]**(p-1))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p+1))))*np.dot(np.array([PolygonUsed['r_n'][1]]).transpose(),np.array([PolygonUsed['r_n'][0]])) + (1-((hyperplane[1]**(p-1))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p))))*np.zeros((2,2))
    betad = (1-((hyperplane[0]**(p-1))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p))))*np.array([PolygonUsed['r_n'][0]]) + (1-((hyperplane[1]**(p-1))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p))))*np.array([PolygonUsed['r_n'][1]])
    beta = hyperplane[0] + hyperplane[1] - (hyperplane[0]**p+hyperplane[1]**p)**(1/p)
    for i in range(2,len(hyperplane)):
        betadd = (-(p-1)*((beta**(p-2))/((beta**p+hyperplane[i]**p)**((p-1)/p)))+(p-1)*((beta**(2*p-2))/((beta**p+hyperplane[i]**p)**((p-1)/p+1))))*np.dot(betad.transpose(),betad) + ((p-1)*((beta**(p-1)*hyperplane[i]**(p-1))/((beta**p+hyperplane[i]**p)**((p-1)/p+1))))*np.dot(betad.transpose(),np.array([PolygonUsed['r_n'][i]])) + (1-((beta**(p-1))/((beta**p+hyperplane[i]**p)**((p-1)/p))))*betadd + (-(p-1)*((hyperplane[i]**(p-2))/((beta**p+hyperplane[i]**p)**((p-1)/p)))+(p-1)*((hyperplane[i]**(2*p-2))/((beta**p+hyperplane[i]**p)**((p-1)/p+1))))*np.dot(np.array([PolygonUsed['r_n'][i]]).transpose(),np.array([PolygonUsed['r_n'][i]])) + ((p-1)*((beta**(p-1)*hyperplane[i]**(p-1))/((beta**p+hyperplane[i]**p)**((p-1)/p+1))))*np.dot(np.array([PolygonUsed['r_n'][i]]).transpose(),betad) + (1-((hyperplane[i]**(p-1))/((beta**p+hyperplane[i]**p)**((p-1)/p))))*np.zeros((2,2))
        betad = (1-((beta**(p-1))/((beta**p+hyperplane[i]**p)**((p-1)/p))))*betad + (1-((hyperplane[i]**(p-1))/((beta**p+hyperplane[i]**p)**((p-1)/p))))*np.array([PolygonUsed['r_n'][i]])
        beta = beta + hyperplane[i] - (beta**p+hyperplane[i]**p)**(1/p)

    return -beta, -betad, -betadd


def insideImplicitTriangle(Position, Triangle, DiffeoParams):
    """
    Function that computes gamma(x) (i.e., the R-function) for a point x inside an enclosing polygon, and its gradient and hessian
    
    Input:
        1) Position: Point to consider - 1x2 numpy.array
        2) Triangle: Description of the triangle - Dictionary
        3) DiffeoParams: Options for the diffeomorphism construction
        
    Output:
        1) gamma: gamma(x)
        2) gammad: Gradient of gamma(x)
        3) gammadd: Hessian of gamma(x)
    """
    # Unwrap parameters
    p = DiffeoParams['p']

    # Compute hyperplane functions
    hyperplane = np.zeros(Triangle['vertices_tilde'].shape[0])
    for i in range(0, Triangle['vertices_tilde'].shape[0]):
        hyperplane[i] = np.dot(Position[0]-Triangle['vertices_tilde'][i],Triangle['r_tilde_n'][i])

    # Compute the R-function and its gradient and hessian
    gammadd = (-(p-1)*((hyperplane[0]**(p-2))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p)))+(p-1)*((hyperplane[0]**(2*p-2))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p+1))))*np.dot(np.array([Triangle['r_tilde_n'][0]]).transpose(),np.array([Triangle['r_tilde_n'][0]])) + ((p-1)*((hyperplane[0]**(p-1)*hyperplane[1]**(p-1))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p+1))))*np.dot(np.array([Triangle['r_tilde_n'][0]]).transpose(),np.array([Triangle['r_tilde_n'][1]])) + (1-((hyperplane[0]**(p-1))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p))))*np.zeros((2,2)) + (-(p-1)*((hyperplane[1]**(p-2))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p)))+(p-1)*((hyperplane[1]**(2*p-2))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p+1))))*np.dot(np.array([Triangle['r_tilde_n'][1]]).transpose(),np.array([Triangle['r_tilde_n'][1]])) + ((p-1)*((hyperplane[0]**(p-1)*hyperplane[1]**(p-1))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p+1))))*np.dot(np.array([Triangle['r_tilde_n'][1]]).transpose(),np.array([Triangle['r_tilde_n'][0]])) + (1-((hyperplane[1]**(p-1))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p))))*np.zeros((2,2))
    gammad = (1-((hyperplane[0]**(p-1))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p))))*np.array([Triangle['r_tilde_n'][0]]) + (1-((hyperplane[1]**(p-1))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p))))*np.array([Triangle['r_tilde_n'][1]])
    gamma = hyperplane[0] + hyperplane[1] - (hyperplane[0]**p+hyperplane[1]**p)**(1/p)
    for i in range(2,len(hyperplane)):
        gammadd = (-(p-1)*((gamma**(p-2))/((gamma**p+hyperplane[i]**p)**((p-1)/p)))+(p-1)*((gamma**(2*p-2))/((gamma**p+hyperplane[i]**p)**((p-1)/p+1))))*np.dot(gammad.transpose(),gammad) + ((p-1)*((gamma**(p-1)*hyperplane[i]**(p-1))/((gamma**p+hyperplane[i]**p)**((p-1)/p+1))))*np.dot(gammad.transpose(),np.array([Triangle['r_tilde_n'][i]])) + (1-((gamma**(p-1))/((gamma**p+hyperplane[i]**p)**((p-1)/p))))*gammadd + (-(p-1)*((hyperplane[i]**(p-2))/((gamma**p+hyperplane[i]**p)**((p-1)/p)))+(p-1)*((hyperplane[i]**(2*p-2))/((gamma**p+hyperplane[i]**p)**((p-1)/p+1))))*np.dot(np.array([Triangle['r_tilde_n'][i]]).transpose(),np.array([Triangle['r_tilde_n'][i]])) + ((p-1)*((gamma**(p-1)*hyperplane[i]**(p-1))/((gamma**p+hyperplane[i]**p)**((p-1)/p+1))))*np.dot(np.array([Triangle['r_tilde_n'][i]]).transpose(),gammad) + (1-((hyperplane[i]**(p-1))/((gamma**p+hyperplane[i]**p)**((p-1)/p))))*np.zeros((2,2))
        gammad = (1-((gamma**(p-1))/((gamma**p+hyperplane[i]**p)**((p-1)/p))))*gammad + (1-((hyperplane[i]**(p-1))/((gamma**p+hyperplane[i]**p)**((p-1)/p))))*np.array([Triangle['r_tilde_n'][i]])
        gamma = gamma + hyperplane[i] - (gamma**p+hyperplane[i]**p)**(1/p)

    return gamma, gammad, gammadd


def insideImplicitPolygon(Position, PolygonUsed, DiffeoParams):
    """
    Function that computes gamma(x) (i.e., the R-function) for a point x inside an enclosing polygon, and its gradient and hessian
    
    Input:
        1) Position: Point to consider - 1x2 numpy.array
        2) PolygonUsed: Description of the polygon - Dictionary
        3) DiffeoParams: Options for the diffeomorphism construction
        
    Output:
        1) gamma: gamma(x)
        2) gammad: Gradient of gamma(x)
        3) gammadd: Hessian of gamma(x)
    """
    # # Unpack the values and compute the hyperplanes
    # PolygonVertices = PolygonUsed['vertices_tilde']
    # PolygonNormals = PolygonUsed['r_tilde_n']
        
    # # Find length of polygon vertices
    # m = PolygonVertices.shape[0]

    # # Span to find values and derivatives
    # hyperplane_values = []
    # hyperplane_values_d = []
    # for k in range(0,m):
    #     hyperplane_values.append(np.dot(Position[0]-PolygonVertices[k], PolygonNormals[k]))
    #     hyperplane_values_d.append(np.array([PolygonNormals[k]]))

    # # Initialize gamma and derivatives
    # gamma = 0.0
    # gamma_d = np.array([[0.0, 0.0]])
    # gamma_dd = np.zeros((2,2))
    # for k in range(0,m):
    #     gamma = gamma + hyperplane_values[k]
    #     gamma_d = gamma_d + hyperplane_values_d[k]

    # # Span through the polygon to find the m by k combinations
    # for k in range(2,m+1):
    #     sorting = list(itertools.combinations(range(0,m),k))
    #     power_k = (-1)**(k+1)

    #     # Span through a specific sort tuple
    #     for sort in sorting:
    #         alpha_k_before_root = 0
    #         alpha_k_before_root_d = np.array([[0.0, 0.0]])
    #         alpha_k_before_root_dd = np.zeros((2,2))

    #         # Span the indices of the sort to find the values, gradients and hessian inside the root - (a_1^2 + ... + a_n^2)
    #         for j in sort:
    #             alpha_k_before_root = alpha_k_before_root + hyperplane_values[j]**2
    #             alpha_k_before_root_d = alpha_k_before_root_d + 2.0*hyperplane_values[j]*hyperplane_values_d[j]
    #             alpha_k_before_root_dd = alpha_k_before_root_dd + 2.0*np.dot(hyperplane_values_d[j].transpose(),hyperplane_values_d[j])
                
    #         # Precompute some terms
    #         alpha_k_before_root_sqrt = np.sqrt(alpha_k_before_root)
                
    #         # Find the final values, gradients and hessians with the root and add it to the current gamma - sqrt(a_1^2 + ... + a_n^2)
    #         gamma = gamma + power_k * alpha_k_before_root_sqrt
    #         gamma_d = gamma_d + power_k * alpha_k_before_root_d/(2.0*alpha_k_before_root_sqrt)
    #         gamma_dd = gamma_dd + power_k*alpha_k_before_root_dd/(2.0*alpha_k_before_root_sqrt) - power_k*np.dot(alpha_k_before_root_d.transpose(),alpha_k_before_root_d)/(4.0*alpha_k_before_root_sqrt**3)

    # return gamma, gamma_d, gamma_dd


    # Unwrap parameters
    p = DiffeoParams['p']

    # Compute hyperplane functions
    hyperplane = np.zeros(PolygonUsed['vertices_tilde'].shape[0])
    for i in range(0, PolygonUsed['vertices_tilde'].shape[0]):
        hyperplane[i] = np.dot(Position[0]-PolygonUsed['vertices_tilde'][i],PolygonUsed['r_tilde_n'][i])

    # Compute the R-function and its gradient and hessian
    gammadd = (-(p-1)*((hyperplane[0]**(p-2))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p)))+(p-1)*((hyperplane[0]**(2*p-2))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p+1))))*np.dot(np.array([PolygonUsed['r_tilde_n'][0]]).transpose(),np.array([PolygonUsed['r_tilde_n'][0]])) + ((p-1)*((hyperplane[0]**(p-1)*hyperplane[1]**(p-1))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p+1))))*np.dot(np.array([PolygonUsed['r_tilde_n'][0]]).transpose(),np.array([PolygonUsed['r_tilde_n'][1]])) + (1-((hyperplane[0]**(p-1))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p))))*np.zeros((2,2)) + (-(p-1)*((hyperplane[1]**(p-2))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p)))+(p-1)*((hyperplane[1]**(2*p-2))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p+1))))*np.dot(np.array([PolygonUsed['r_tilde_n'][1]]).transpose(),np.array([PolygonUsed['r_tilde_n'][1]])) + ((p-1)*((hyperplane[0]**(p-1)*hyperplane[1]**(p-1))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p+1))))*np.dot(np.array([PolygonUsed['r_tilde_n'][1]]).transpose(),np.array([PolygonUsed['r_tilde_n'][0]])) + (1-((hyperplane[1]**(p-1))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p))))*np.zeros((2,2))
    gammad = (1-((hyperplane[0]**(p-1))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p))))*np.array([PolygonUsed['r_tilde_n'][0]]) + (1-((hyperplane[1]**(p-1))/((hyperplane[0]**p+hyperplane[1]**p)**((p-1)/p))))*np.array([PolygonUsed['r_tilde_n'][1]])
    gamma = hyperplane[0] + hyperplane[1] - (hyperplane[0]**p+hyperplane[1]**p)**(1/p)
    for i in range(2,len(hyperplane)):
        gammadd = (-(p-1)*((gamma**(p-2))/((gamma**p+hyperplane[i]**p)**((p-1)/p)))+(p-1)*((gamma**(2*p-2))/((gamma**p+hyperplane[i]**p)**((p-1)/p+1))))*np.dot(gammad.transpose(),gammad) + ((p-1)*((gamma**(p-1)*hyperplane[i]**(p-1))/((gamma**p+hyperplane[i]**p)**((p-1)/p+1))))*np.dot(gammad.transpose(),np.array([PolygonUsed['r_tilde_n'][i]])) + (1-((gamma**(p-1))/((gamma**p+hyperplane[i]**p)**((p-1)/p))))*gammadd + (-(p-1)*((hyperplane[i]**(p-2))/((gamma**p+hyperplane[i]**p)**((p-1)/p)))+(p-1)*((hyperplane[i]**(2*p-2))/((gamma**p+hyperplane[i]**p)**((p-1)/p+1))))*np.dot(np.array([PolygonUsed['r_tilde_n'][i]]).transpose(),np.array([PolygonUsed['r_tilde_n'][i]])) + ((p-1)*((gamma**(p-1)*hyperplane[i]**(p-1))/((gamma**p+hyperplane[i]**p)**((p-1)/p+1))))*np.dot(np.array([PolygonUsed['r_tilde_n'][i]]).transpose(),gammad) + (1-((hyperplane[i]**(p-1))/((gamma**p+hyperplane[i]**p)**((p-1)/p))))*np.zeros((2,2))
        gammad = (1-((gamma**(p-1))/((gamma**p+hyperplane[i]**p)**((p-1)/p))))*gammad + (1-((hyperplane[i]**(p-1))/((gamma**p+hyperplane[i]**p)**((p-1)/p))))*np.array([PolygonUsed['r_tilde_n'][i]])
        gamma = gamma + hyperplane[i] - (gamma**p+hyperplane[i]**p)**(1/p)

    return gamma, gammad, gammadd


def polygonImplicit(Position, DiffeoTree, DiffeoParams):
    """
    Function that computes b(x) (i.e., the implicit function corresponding to a polygon), and its derivative and hessian
    
    Input:
        1) Position: Point to consider - 1x2 numpy.array
        2) DiffeoTree: Tree that contains the diffeomorphism properties for a particular polygon
        3) DiffeoParams: Options for the diffeomorphism construction
        
    Output:
        1) b: Value of the implicit function
        2) bd: Value of the implicit function gradient - 1x2 numpy.array
        3) bdd: Values of the implicit function hessian - 2x2 numpy.array
    """
    # Begin purging process with default values
    if Position.shape == (2,):
        Position = np.array([Position])
    
    # Unwrap parameters
    p = DiffeoParams['p']

    # Compute the R-function and its gradient and Hessian
    b, bd, bdd = outsideImplicitTriangle(Position, DiffeoTree[0], DiffeoParams)
    for i in range(1,len(DiffeoTree)):
        btemp, btempd, btempdd = outsideImplicitTriangle(Position, DiffeoTree[i], DiffeoParams)
        bdd = (-(p-1)*((b**(p-2))/((b**p+btemp**p)**((p-1)/p)))+(p-1)*((b**(2*p-2))/((b**p+btemp**p)**((p-1)/p+1))))*np.dot(bd.transpose(),bd) + ((p-1)*((b**(p-1)*btemp**(p-1))/((b**p+btemp**p)**((p-1)/p+1))))*np.dot(bd.transpose(),btempd) + (1-((b**(p-1))/((b**p+btemp**p)**((p-1)/p))))*bdd + (-(p-1)*((btemp**(p-2))/((b**p+btemp**p)**((p-1)/p)))+(p-1)*((btemp**(2*p-2))/((b**p+btemp**p)**((p-1)/p+1))))*np.dot(btempd.transpose(),btempd) + ((p-1)*((b**(p-1)*btemp**(p-1))/((b**p+btemp**p)**((p-1)/p+1))))*np.dot(btempd.transpose(),bd) + (1-((btemp**(p-1))/((b**p+btemp**p)**((p-1)/p))))*btempdd
        bd = (1-((b**(p-1))/((b**p+btemp**p)**((p-1)/p))))*bd + (1-((btemp**(p-1))/((b**p+btemp**p)**((p-1)/p))))*btempd
        b = b + btemp - (b**p+btemp**p)**(1/p)
        
    return b, bd, bdd


def butterworthLowPass(cutoff, fs, order=5):
    """
    Function that generates the filter coefficients based on a cutoff frequency and a sampling rate
    
    Input:
        1) cutoff: Cutoff frequency
        2) fs: Sampling frequency
        3) order: Order of the filter
        
    Output:
        1) b: Numerator coefficient array
        2) a: Denominator coefficient array
    """
    nyq = 0.5 * fs
    normal_cutoff = cutoff/nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butterworthLowPassFilter(data, cutoff, fs, order=5):
    """
    Function that filters the desired data based on a cutoff frequency and a sampling frequency
    
    Input:
        1) data: Data array
        2) cutoff: Cutoff frequency
        3) fs: Sampling frequency
        4) order: Order of the filter
    
    Output:
        1) y: Filtered version of the input data
    """
    b, a = butterworthLowPass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y
