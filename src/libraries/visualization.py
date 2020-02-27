#!/usr/bin/env python

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
import struct, math, numpy, os, sys, scipy, time, random
import matplotlib.pyplot as plt
import shapely as sp
from matplotlib.animation import FuncAnimation
from shapely.geometry import Point, LineString
from shapely.geometry.polygon import Polygon

# Reactive planner imports
from reactive_planner_lib import LIDARClass, completeLIDAR2D, compensateObstacleLIDAR2D, readLIDAR2D
from reactive_planner_lib import polygonDiffeoTriangulation, polygonDiffeoConvex, diffeoTreeTriangulation, diffeoTreeConvex, triangleDiffeo, polygonDiffeo, polygonImplicit, triangleSwitch, polygonSwitch
from reactive_planner_lib import localfreespaceLIDAR2D

def visualize_diffeoDeterminant_triangulation(Polygons, RobotRadius, PlotBounds, NumPoints, DiffeoParams):
	"""
	Function that visualizes the determinant of the diffeomorphism (based on the ear clipping method) on the plane, given a set of polygons and a robot radius
	
	Input:
		1) Polygons: Vertex Coordinates of input polygons - M-member list of Nx2 numpy.array objects (start and end vertices must be the same)
		2) RobotRadius: Robot radius (m)
		3) PlotBounds: Bounds for the planar plot - 4-member numpy.array ([xmin, xmax, ymin, ymax])
		4) NumPoints: Number of points for the generated grid in x and y - 2-member numpy.array ([x_resolution, y_resolution])
        5) DiffeoParams: Options for the diffeomorphism construction
	
	Test:
		import numpy
		import visualization
		robot_radius = 0.25
		bounds = numpy.array([0, 5, -3, 3])
		num_points = numpy.array([101, 101])
		polygon_list = []
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		polygon_list.append(xy)
		diffeo_params = dict()
		diffeo_params['p'] = 20
		diffeo_params['epsilon'] = 1.5
		diffeo_params['varepsilon'] = 1.5
		diffeo_params['mu_1'] = 1.5
		diffeo_params['mu_2'] = 0.01
		diffeo_params['workspace'] = numpy.array([[-100,-100],[100,-100],[100,100],[-100,100],[-100,-100]])
		visualization.visualize_diffeoDeterminant_triangulation(polygon_list, robot_radius, bounds, num_points, diffeo_params)
	
	Polygon examples to test:
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		xy = numpy.array([[0,5,5,0,0,4,4,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[2.8,6.2,6.2,2.8,2.8,4.8,4.8,2.8,2.8,4.8,4.8,2.8,2.8], [1.4,1.4,10,10,8.6,8.6,7,7,4.4,4.4,2.8,2.8,1.4]]).transpose()
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		xy = numpy.array([[7,9.5,10,9,9,8,9,10,11,10,9,7], [7,8,7,6,7,7,5,5,7,9,9,7]]).transpose()
		xy = numpy.array([[7,7,8,8,7,7,10,10,9,9,10,10,7], [7,6,6,1,1,0,0,1,1,6,6,7,7]]).transpose()
		xy = numpy.array([[0,10,10,0,0,9,9,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[0,0.5,0.5,1.5,1.5,-1,-1,3,3,2,2,0,0], [0,0,-1,-1,1,1,-3,-3,1,1,-2,-2,0]]).transpose()
		xy = numpy.vstack((sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[0],sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[1])).transpose()
	"""
	# Construct list of polygonal objects and enlarge by robot radius
	polygon_list = []
	for i in range(len(Polygons)):
		polygon_list.append(Polygon(Polygons[i]).buffer(RobotRadius, join_style=2))
	
	# Span all the found polygons to check for intersections between the known polygons and keep only the merged polygons
	polygon_list_merged = []
	i = 0
	while (i<len(polygon_list)):
		polygon_list_merged.append(polygon_list[i])

		j = i+1
		while (j<len(polygon_list)):
			if polygon_list_merged[i].intersects(polygon_list[j]):
				polygon_list_merged[i] = polygon_list_merged[i].union(polygon_list[j])
				polygon_list_merged[i] = polygon_list_merged[i].simplify(0.08, preserve_topology=True) # simplify polygon to eliminate strange small corners
				del(polygon_list[j])
			else:
				j = j+1
		polygon_list_merged[i] = sp.geometry.polygon.orient(polygon_list_merged[i], 1.0) # orient polygon to be CCW
		i = i+1
	PolygonList = polygon_list_merged

	# Construct list of diffeo trees for all objects
	DiffeoTreeArray = []
	for i in range(len(polygon_list_merged)):
		coords = numpy.vstack((polygon_list_merged[i].exterior.coords.xy[0],polygon_list_merged[i].exterior.coords.xy[1])).transpose()
		DiffeoTreeArray.append(diffeoTreeTriangulation(coords, DiffeoParams))

	# Generate x and y coordinates
	x_coords = numpy.linspace(PlotBounds[0], PlotBounds[1], NumPoints[0])
	y_coords = numpy.linspace(PlotBounds[2], PlotBounds[3], NumPoints[1])

	# Span all the points
	data_points = numpy.zeros((y_coords.shape[0],x_coords.shape[0]))
	for j in range(y_coords.shape[0]):
		for i in range(x_coords.shape[0]):
			candidate_point = Point(x_coords[i],y_coords[j])

			# Check for inclusion in any of the polygons
			for k in range(len(polygon_list_merged)):
				if polygon_list_merged[k].contains(candidate_point):
					data_points[j][i] = numpy.NAN
					collision = True
					break
				else:
					collision = False
			
			if collision is True:
				continue
			else:
				# Compute the actual diffeomorphism
				PositionTransformed = numpy.array([[x_coords[i],y_coords[j]]])
				PositionTransformedD = numpy.eye(2)
				PositionTransformedDD = numpy.zeros(8)
				for k in range(len(DiffeoTreeArray)):
					TempPositionTransformed, TempPositionTransformedD, TempPositionTransformedDD = polygonDiffeoTriangulation(PositionTransformed, DiffeoTreeArray[k], DiffeoParams)

					res1 = TempPositionTransformedD[0][0]*PositionTransformedDD[0] + TempPositionTransformedD[0][1]*PositionTransformedDD[4] + PositionTransformedD[0][0]*(TempPositionTransformedDD[0]*PositionTransformedD[0][0] + TempPositionTransformedDD[1]*PositionTransformedD[1][0]) + PositionTransformedD[1][0]*(TempPositionTransformedDD[2]*PositionTransformedD[0][0] + TempPositionTransformedDD[3]*PositionTransformedD[1][0])
					res2 = TempPositionTransformedD[0][0]*PositionTransformedDD[1] + TempPositionTransformedD[0][1]*PositionTransformedDD[5] + PositionTransformedD[0][0]*(TempPositionTransformedDD[0]*PositionTransformedD[0][1] + TempPositionTransformedDD[1]*PositionTransformedD[1][1]) + PositionTransformedD[1][0]*(TempPositionTransformedDD[2]*PositionTransformedD[0][1] + TempPositionTransformedDD[3]*PositionTransformedD[1][1])
					res3 = TempPositionTransformedD[0][0]*PositionTransformedDD[2] + TempPositionTransformedD[0][1]*PositionTransformedDD[6] + PositionTransformedD[0][1]*(TempPositionTransformedDD[0]*PositionTransformedD[0][0] + TempPositionTransformedDD[1]*PositionTransformedD[1][0]) + PositionTransformedD[1][1]*(TempPositionTransformedDD[2]*PositionTransformedD[0][0] + TempPositionTransformedDD[3]*PositionTransformedD[1][0])
					res4 = TempPositionTransformedD[0][0]*PositionTransformedDD[3] + TempPositionTransformedD[0][1]*PositionTransformedDD[7] + PositionTransformedD[0][1]*(TempPositionTransformedDD[0]*PositionTransformedD[0][1] + TempPositionTransformedDD[1]*PositionTransformedD[1][1]) + PositionTransformedD[1][1]*(TempPositionTransformedDD[2]*PositionTransformedD[0][1] + TempPositionTransformedDD[3]*PositionTransformedD[1][1])
					res5 = TempPositionTransformedD[1][0]*PositionTransformedDD[0] + TempPositionTransformedD[1][1]*PositionTransformedDD[4] + PositionTransformedD[0][0]*(TempPositionTransformedDD[4]*PositionTransformedD[0][0] + TempPositionTransformedDD[5]*PositionTransformedD[1][0]) + PositionTransformedD[1][0]*(TempPositionTransformedDD[6]*PositionTransformedD[0][0] + TempPositionTransformedDD[7]*PositionTransformedD[1][0])
					res6 = TempPositionTransformedD[1][0]*PositionTransformedDD[1] + TempPositionTransformedD[1][1]*PositionTransformedDD[5] + PositionTransformedD[0][0]*(TempPositionTransformedDD[4]*PositionTransformedD[0][1] + TempPositionTransformedDD[5]*PositionTransformedD[1][1]) + PositionTransformedD[1][0]*(TempPositionTransformedDD[6]*PositionTransformedD[0][1] + TempPositionTransformedDD[7]*PositionTransformedD[1][1])
					res7 = TempPositionTransformedD[1][0]*PositionTransformedDD[2] + TempPositionTransformedD[1][1]*PositionTransformedDD[6] + PositionTransformedD[0][1]*(TempPositionTransformedDD[4]*PositionTransformedD[0][0] + TempPositionTransformedDD[5]*PositionTransformedD[1][0]) + PositionTransformedD[1][1]*(TempPositionTransformedDD[6]*PositionTransformedD[0][0] + TempPositionTransformedDD[7]*PositionTransformedD[1][0])
					res8 = TempPositionTransformedD[1][0]*PositionTransformedDD[3] + TempPositionTransformedD[1][1]*PositionTransformedDD[7] + PositionTransformedD[0][1]*(TempPositionTransformedDD[4]*PositionTransformedD[0][1] + TempPositionTransformedDD[5]*PositionTransformedD[1][1]) + PositionTransformedD[1][1]*(TempPositionTransformedDD[6]*PositionTransformedD[0][1] + TempPositionTransformedDD[7]*PositionTransformedD[1][1])
					PositionTransformedDD[0] = res1
					PositionTransformedDD[1] = res2
					PositionTransformedDD[2] = res3
					PositionTransformedDD[3] = res4
					PositionTransformedDD[4] = res5
					PositionTransformedDD[5] = res6
					PositionTransformedDD[6] = res7
					PositionTransformedDD[7] = res8

					PositionTransformedD = numpy.matmul(TempPositionTransformedD, PositionTransformedD)

					PositionTransformed = TempPositionTransformed
					
				# Add the data point
				data_points[j][i] = numpy.linalg.det(PositionTransformedD)
				print(i+j*y_coords.shape[0])
	
	# Plot the result
	plt.imshow(numpy.log(data_points), vmin=numpy.log(data_points[~numpy.isnan(data_points)]).min(), vmax=numpy.log(data_points[~numpy.isnan(data_points)]).max(), origin='lower', extent=[PlotBounds[0], PlotBounds[1], PlotBounds[2], PlotBounds[3]])
	plt.axis('off')
	plt.colorbar()
	plt.show()

	return


def visualize_diffeoDeterminant_convex(Polygons, RobotRadius, PlotBounds, NumPoints, DiffeoParams):
	"""
	Function that visualizes the determinant of the diffeomorphism (based on the convex decomposition method) on the plane, given a set of polygons and a robot radius
	
	Input:
		1) Polygons: Vertex Coordinates of input polygons - M-member list of Nx2 numpy.array objects (start and end vertices must be the same)
		2) RobotRadius: Robot radius (m)
		3) PlotBounds: Bounds for the planar plot - 4-member numpy.array ([xmin, xmax, ymin, ymax])
		4) NumPoints: Number of points for the generated grid in x and y - 2-member numpy.array ([x_resolution, y_resolution])
        5) DiffeoParams: Options for the diffeomorphism construction
	
	Test:
		import numpy
		import visualization
		robot_radius = 0.25
		bounds = numpy.array([0, 5, -3, 3])
		num_points = numpy.array([101, 101])
		polygon_list = []
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		polygon_list.append(xy)
		diffeo_params = dict()
		diffeo_params['p'] = 20
		diffeo_params['epsilon'] = 1.5
		diffeo_params['varepsilon'] = 1.5
		diffeo_params['mu_1'] = 1.5
		diffeo_params['mu_2'] = 0.01
		diffeo_params['workspace'] = numpy.array([[-100,-100],[100,-100],[100,100],[-100,100],[-100,-100]])
		visualization.visualize_diffeoDeterminant_convex(polygon_list, robot_radius, bounds, num_points, diffeo_params)
	
	Polygon examples to test:
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		xy = numpy.array([[0,5,5,0,0,4,4,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[2.8,6.2,6.2,2.8,2.8,4.8,4.8,2.8,2.8,4.8,4.8,2.8,2.8], [1.4,1.4,10,10,8.6,8.6,7,7,4.4,4.4,2.8,2.8,1.4]]).transpose()
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		xy = numpy.array([[7,9.5,10,9,9,8,9,10,11,10,9,7], [7,8,7,6,7,7,5,5,7,9,9,7]]).transpose()
		xy = numpy.array([[7,7,8,8,7,7,10,10,9,9,10,10,7], [7,6,6,1,1,0,0,1,1,6,6,7,7]]).transpose()
		xy = numpy.array([[0,10,10,0,0,9,9,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[0,0.5,0.5,1.5,1.5,-1,-1,3,3,2,2,0,0], [0,0,-1,-1,1,1,-3,-3,1,1,-2,-2,0]]).transpose()
		xy = numpy.vstack((sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[0],sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[1])).transpose()
	"""
	# Construct list of polygonal objects and enlarge by robot radius
	polygon_list = []
	for i in range(len(Polygons)):
		polygon_list.append(Polygon(Polygons[i]).buffer(RobotRadius, join_style=2))
	
	# Span all the found polygons to check for intersections between the known polygons and keep only the merged polygons
	polygon_list_merged = []
	i = 0
	while (i<len(polygon_list)):
		polygon_list_merged.append(polygon_list[i])

		j = i+1
		while (j<len(polygon_list)):
			if polygon_list_merged[i].intersects(polygon_list[j]):
				polygon_list_merged[i] = polygon_list_merged[i].union(polygon_list[j])
				polygon_list_merged[i] = polygon_list_merged[i].simplify(0.08, preserve_topology=True) # simplify polygon to eliminate strange small corners
				del(polygon_list[j])
			else:
				j = j+1
		polygon_list_merged[i] = sp.geometry.polygon.orient(polygon_list_merged[i], 1.0) # orient polygon to be CCW
		i = i+1
	PolygonList = polygon_list_merged

	# Construct list of diffeo trees for all objects
	DiffeoTreeArray = []
	for i in range(len(polygon_list_merged)):
		coords = numpy.vstack((polygon_list_merged[i].exterior.coords.xy[0],polygon_list_merged[i].exterior.coords.xy[1])).transpose()
		DiffeoTreeArray.append(diffeoTreeConvex(coords, DiffeoParams))

	# Generate x and y coordinates
	x_coords = numpy.linspace(PlotBounds[0], PlotBounds[1], NumPoints[0])
	y_coords = numpy.linspace(PlotBounds[2], PlotBounds[3], NumPoints[1])

	# Span all the points
	data_points = numpy.zeros((y_coords.shape[0],x_coords.shape[0]))
	for j in range(y_coords.shape[0]):
		for i in range(x_coords.shape[0]):
			candidate_point = Point(x_coords[i],y_coords[j])

			# Check for inclusion in any of the polygons
			for k in range(len(polygon_list_merged)):
				if polygon_list_merged[k].contains(candidate_point):
					data_points[j][i] = numpy.NAN
					collision = True
					break
				else:
					collision = False
			
			if collision is True:
				continue
			else:
				# Compute the actual diffeomorphism
				PositionTransformed = numpy.array([[x_coords[i],y_coords[j]]])
				PositionTransformedD = numpy.eye(2)
				PositionTransformedDD = numpy.zeros(8)
				for k in range(len(DiffeoTreeArray)):
					TempPositionTransformed, TempPositionTransformedD, TempPositionTransformedDD = polygonDiffeoConvex(PositionTransformed, DiffeoTreeArray[k], DiffeoParams)

					res1 = TempPositionTransformedD[0][0]*PositionTransformedDD[0] + TempPositionTransformedD[0][1]*PositionTransformedDD[4] + PositionTransformedD[0][0]*(TempPositionTransformedDD[0]*PositionTransformedD[0][0] + TempPositionTransformedDD[1]*PositionTransformedD[1][0]) + PositionTransformedD[1][0]*(TempPositionTransformedDD[2]*PositionTransformedD[0][0] + TempPositionTransformedDD[3]*PositionTransformedD[1][0])
					res2 = TempPositionTransformedD[0][0]*PositionTransformedDD[1] + TempPositionTransformedD[0][1]*PositionTransformedDD[5] + PositionTransformedD[0][0]*(TempPositionTransformedDD[0]*PositionTransformedD[0][1] + TempPositionTransformedDD[1]*PositionTransformedD[1][1]) + PositionTransformedD[1][0]*(TempPositionTransformedDD[2]*PositionTransformedD[0][1] + TempPositionTransformedDD[3]*PositionTransformedD[1][1])
					res3 = TempPositionTransformedD[0][0]*PositionTransformedDD[2] + TempPositionTransformedD[0][1]*PositionTransformedDD[6] + PositionTransformedD[0][1]*(TempPositionTransformedDD[0]*PositionTransformedD[0][0] + TempPositionTransformedDD[1]*PositionTransformedD[1][0]) + PositionTransformedD[1][1]*(TempPositionTransformedDD[2]*PositionTransformedD[0][0] + TempPositionTransformedDD[3]*PositionTransformedD[1][0])
					res4 = TempPositionTransformedD[0][0]*PositionTransformedDD[3] + TempPositionTransformedD[0][1]*PositionTransformedDD[7] + PositionTransformedD[0][1]*(TempPositionTransformedDD[0]*PositionTransformedD[0][1] + TempPositionTransformedDD[1]*PositionTransformedD[1][1]) + PositionTransformedD[1][1]*(TempPositionTransformedDD[2]*PositionTransformedD[0][1] + TempPositionTransformedDD[3]*PositionTransformedD[1][1])
					res5 = TempPositionTransformedD[1][0]*PositionTransformedDD[0] + TempPositionTransformedD[1][1]*PositionTransformedDD[4] + PositionTransformedD[0][0]*(TempPositionTransformedDD[4]*PositionTransformedD[0][0] + TempPositionTransformedDD[5]*PositionTransformedD[1][0]) + PositionTransformedD[1][0]*(TempPositionTransformedDD[6]*PositionTransformedD[0][0] + TempPositionTransformedDD[7]*PositionTransformedD[1][0])
					res6 = TempPositionTransformedD[1][0]*PositionTransformedDD[1] + TempPositionTransformedD[1][1]*PositionTransformedDD[5] + PositionTransformedD[0][0]*(TempPositionTransformedDD[4]*PositionTransformedD[0][1] + TempPositionTransformedDD[5]*PositionTransformedD[1][1]) + PositionTransformedD[1][0]*(TempPositionTransformedDD[6]*PositionTransformedD[0][1] + TempPositionTransformedDD[7]*PositionTransformedD[1][1])
					res7 = TempPositionTransformedD[1][0]*PositionTransformedDD[2] + TempPositionTransformedD[1][1]*PositionTransformedDD[6] + PositionTransformedD[0][1]*(TempPositionTransformedDD[4]*PositionTransformedD[0][0] + TempPositionTransformedDD[5]*PositionTransformedD[1][0]) + PositionTransformedD[1][1]*(TempPositionTransformedDD[6]*PositionTransformedD[0][0] + TempPositionTransformedDD[7]*PositionTransformedD[1][0])
					res8 = TempPositionTransformedD[1][0]*PositionTransformedDD[3] + TempPositionTransformedD[1][1]*PositionTransformedDD[7] + PositionTransformedD[0][1]*(TempPositionTransformedDD[4]*PositionTransformedD[0][1] + TempPositionTransformedDD[5]*PositionTransformedD[1][1]) + PositionTransformedD[1][1]*(TempPositionTransformedDD[6]*PositionTransformedD[0][1] + TempPositionTransformedDD[7]*PositionTransformedD[1][1])
					PositionTransformedDD[0] = res1
					PositionTransformedDD[1] = res2
					PositionTransformedDD[2] = res3
					PositionTransformedDD[3] = res4
					PositionTransformedDD[4] = res5
					PositionTransformedDD[5] = res6
					PositionTransformedDD[6] = res7
					PositionTransformedDD[7] = res8

					PositionTransformedD = numpy.matmul(TempPositionTransformedD, PositionTransformedD)

					PositionTransformed = TempPositionTransformed
					
				# Add the data point
				data_points[j][i] = numpy.linalg.det(PositionTransformedD)
				print(i+j*y_coords.shape[0])
	
	# Plot the result
	plt.imshow(numpy.log(data_points), vmin=numpy.log(data_points[~numpy.isnan(data_points)]).min(), vmax=numpy.log(data_points[~numpy.isnan(data_points)]).max(), origin='lower', extent=[PlotBounds[0], PlotBounds[1], PlotBounds[2], PlotBounds[3]])
	plt.axis('off')
	plt.colorbar()
	plt.show()

	return


def visualize_lyapunov_triangulation(Polygons, RobotRadius, PlotBounds, NumPoints, Goal, DiffeoParams):
	"""
	Function that visualizes the determinant of the diffeomorphism on the plane (based on the ear clipping method), given a set of polygons and a robot radius
	
	Input:
		1) Polygons: Vertex Coordinates of input polygons - M-member list of Nx2 numpy.array objects (start and end vertices must be the same)
		2) RobotRadius: Robot radius (m)
		3) PlotBounds: Bounds for the planar plot - 4-member numpy.array ([xmin, xmax, ymin, ymax])
		4) NumPoints: Number of points for the generated grid in x and y - 2-member numpy.array ([x_resolution, y_resolution])
		5) Goal: The desired navigation goal - 1x2 numpy.array
        6) DiffeoParams: Options for the diffeomorphism construction
	
	Test:
		import numpy
		import visualization
		robot_radius = 0.25
		bounds = numpy.array([0, 5, -3, 3])
		num_points = numpy.array([101, 101])
		goal = numpy.array([[0.0,0.0]])
		polygon_list = []
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		polygon_list.append(xy)
		diffeo_params = dict()
		diffeo_params['p'] = 20
		diffeo_params['epsilon'] = 1.5
		diffeo_params['varepsilon'] = 1.5
		diffeo_params['mu_1'] = 1.5
		diffeo_params['mu_2'] = 0.01
		diffeo_params['workspace'] = numpy.array([[-100,-100],[100,-100],[100,100],[-100,100],[-100,-100]])
		visualization.visualize_lyapunov_triangulation(polygon_list, robot_radius, bounds, num_points, goal, diffeo_params)
	
	Polygon examples to test:
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		xy = numpy.array([[0,5,5,0,0,4,4,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[2.8,6.2,6.2,2.8,2.8,4.8,4.8,2.8,2.8,4.8,4.8,2.8,2.8], [1.4,1.4,10,10,8.6,8.6,7,7,4.4,4.4,2.8,2.8,1.4]]).transpose()
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		xy = numpy.array([[7,9.5,10,9,9,8,9,10,11,10,9,7], [7,8,7,6,7,7,5,5,7,9,9,7]]).transpose()
		xy = numpy.array([[7,7,8,8,7,7,10,10,9,9,10,10,7], [7,6,6,1,1,0,0,1,1,6,6,7,7]]).transpose()
		xy = numpy.array([[0,10,10,0,0,9,9,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[0,0.5,0.5,1.5,1.5,-1,-1,3,3,2,2,0,0], [0,0,-1,-1,1,1,-3,-3,1,1,-2,-2,0]]).transpose()
		xy = numpy.vstack((sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[0],sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[1])).transpose()
	"""
	# Construct list of polygonal objects and enlarge by robot radius
	polygon_list = []
	for i in range(len(Polygons)):
		polygon_list.append(Polygon(Polygons[i]).buffer(RobotRadius, join_style=2))
	
	# Span all the found polygons to check for intersections between the known polygons and keep only the merged polygons
	polygon_list_merged = []
	i = 0
	while (i<len(polygon_list)):
		polygon_list_merged.append(polygon_list[i])

		j = i+1
		while (j<len(polygon_list)):
			if polygon_list_merged[i].intersects(polygon_list[j]):
				polygon_list_merged[i] = polygon_list_merged[i].union(polygon_list[j])
				polygon_list_merged[i] = polygon_list_merged[i].simplify(0.08, preserve_topology=True) # simplify polygon to eliminate strange small corners
				del(polygon_list[j])
			else:
				j = j+1
		polygon_list_merged[i] = sp.geometry.polygon.orient(polygon_list_merged[i], 1.0) # orient polygon to be CCW
		i = i+1
	PolygonList = polygon_list_merged

	# Construct list of diffeo trees for all objects
	DiffeoTreeArray = []
	for i in range(len(polygon_list_merged)):
		coords = numpy.vstack((polygon_list_merged[i].exterior.coords.xy[0],polygon_list_merged[i].exterior.coords.xy[1])).transpose()
		DiffeoTreeArray.append(diffeoTreeTriangulation(coords, DiffeoParams))

	# Generate x and y coordinates
	x_coords = numpy.linspace(PlotBounds[0], PlotBounds[1], NumPoints[0])
	y_coords = numpy.linspace(PlotBounds[2], PlotBounds[3], NumPoints[1])

	# Span all the points
	data_points = numpy.zeros((y_coords.shape[0],x_coords.shape[0]))
	for j in range(y_coords.shape[0]):
		for i in range(x_coords.shape[0]):
			candidate_point = Point(x_coords[i],y_coords[j])

			# Check for inclusion in any of the polygons
			for k in range(len(polygon_list_merged)):
				if polygon_list_merged[k].contains(candidate_point):
					data_points[j][i] = numpy.NAN
					collision = True
					break
				else:
					collision = False
			
			if collision is True:
				continue
			else:
				# Compute the actual diffeomorphism
				PositionTransformed = numpy.array([[x_coords[i],y_coords[j]]])
				PositionTransformedD = numpy.eye(2)
				PositionTransformedDD = numpy.zeros(8)
				GoalTransformed = Goal
				for k in range(len(DiffeoTreeArray)):
					TempPositionTransformed, TempPositionTransformedD, TempPositionTransformedDD = polygonDiffeoTriangulation(PositionTransformed, DiffeoTreeArray[k], DiffeoParams)
					TempGoalTransformed, TempGoalTransformedD, TempGoalTransformedDD = polygonDiffeoTriangulation(GoalTransformed, DiffeoTreeArray[k], DiffeoParams)

					res1 = TempPositionTransformedD[0][0]*PositionTransformedDD[0] + TempPositionTransformedD[0][1]*PositionTransformedDD[4] + PositionTransformedD[0][0]*(TempPositionTransformedDD[0]*PositionTransformedD[0][0] + TempPositionTransformedDD[1]*PositionTransformedD[1][0]) + PositionTransformedD[1][0]*(TempPositionTransformedDD[2]*PositionTransformedD[0][0] + TempPositionTransformedDD[3]*PositionTransformedD[1][0])
					res2 = TempPositionTransformedD[0][0]*PositionTransformedDD[1] + TempPositionTransformedD[0][1]*PositionTransformedDD[5] + PositionTransformedD[0][0]*(TempPositionTransformedDD[0]*PositionTransformedD[0][1] + TempPositionTransformedDD[1]*PositionTransformedD[1][1]) + PositionTransformedD[1][0]*(TempPositionTransformedDD[2]*PositionTransformedD[0][1] + TempPositionTransformedDD[3]*PositionTransformedD[1][1])
					res3 = TempPositionTransformedD[0][0]*PositionTransformedDD[2] + TempPositionTransformedD[0][1]*PositionTransformedDD[6] + PositionTransformedD[0][1]*(TempPositionTransformedDD[0]*PositionTransformedD[0][0] + TempPositionTransformedDD[1]*PositionTransformedD[1][0]) + PositionTransformedD[1][1]*(TempPositionTransformedDD[2]*PositionTransformedD[0][0] + TempPositionTransformedDD[3]*PositionTransformedD[1][0])
					res4 = TempPositionTransformedD[0][0]*PositionTransformedDD[3] + TempPositionTransformedD[0][1]*PositionTransformedDD[7] + PositionTransformedD[0][1]*(TempPositionTransformedDD[0]*PositionTransformedD[0][1] + TempPositionTransformedDD[1]*PositionTransformedD[1][1]) + PositionTransformedD[1][1]*(TempPositionTransformedDD[2]*PositionTransformedD[0][1] + TempPositionTransformedDD[3]*PositionTransformedD[1][1])
					res5 = TempPositionTransformedD[1][0]*PositionTransformedDD[0] + TempPositionTransformedD[1][1]*PositionTransformedDD[4] + PositionTransformedD[0][0]*(TempPositionTransformedDD[4]*PositionTransformedD[0][0] + TempPositionTransformedDD[5]*PositionTransformedD[1][0]) + PositionTransformedD[1][0]*(TempPositionTransformedDD[6]*PositionTransformedD[0][0] + TempPositionTransformedDD[7]*PositionTransformedD[1][0])
					res6 = TempPositionTransformedD[1][0]*PositionTransformedDD[1] + TempPositionTransformedD[1][1]*PositionTransformedDD[5] + PositionTransformedD[0][0]*(TempPositionTransformedDD[4]*PositionTransformedD[0][1] + TempPositionTransformedDD[5]*PositionTransformedD[1][1]) + PositionTransformedD[1][0]*(TempPositionTransformedDD[6]*PositionTransformedD[0][1] + TempPositionTransformedDD[7]*PositionTransformedD[1][1])
					res7 = TempPositionTransformedD[1][0]*PositionTransformedDD[2] + TempPositionTransformedD[1][1]*PositionTransformedDD[6] + PositionTransformedD[0][1]*(TempPositionTransformedDD[4]*PositionTransformedD[0][0] + TempPositionTransformedDD[5]*PositionTransformedD[1][0]) + PositionTransformedD[1][1]*(TempPositionTransformedDD[6]*PositionTransformedD[0][0] + TempPositionTransformedDD[7]*PositionTransformedD[1][0])
					res8 = TempPositionTransformedD[1][0]*PositionTransformedDD[3] + TempPositionTransformedD[1][1]*PositionTransformedDD[7] + PositionTransformedD[0][1]*(TempPositionTransformedDD[4]*PositionTransformedD[0][1] + TempPositionTransformedDD[5]*PositionTransformedD[1][1]) + PositionTransformedD[1][1]*(TempPositionTransformedDD[6]*PositionTransformedD[0][1] + TempPositionTransformedDD[7]*PositionTransformedD[1][1])
					PositionTransformedDD[0] = res1
					PositionTransformedDD[1] = res2
					PositionTransformedDD[2] = res3
					PositionTransformedDD[3] = res4
					PositionTransformedDD[4] = res5
					PositionTransformedDD[5] = res6
					PositionTransformedDD[6] = res7
					PositionTransformedDD[7] = res8

					PositionTransformedD = numpy.matmul(TempPositionTransformedD, PositionTransformedD)

					PositionTransformed = TempPositionTransformed
					GoalTransformed = TempGoalTransformed
					
				# Add the data point
				data_points[j][i] = numpy.linalg.norm(PositionTransformed[0]-GoalTransformed[0])
	
	# Plot the result
	plt.imshow(data_points, vmin=data_points[~numpy.isnan(data_points)].min(), vmax=data_points[~numpy.isnan(data_points)].max(), origin='lower', extent=[PlotBounds[0], PlotBounds[1], PlotBounds[2], PlotBounds[3]])
	plt.colorbar()
	plt.show()

	return


def visualize_diffeoSwitch_triangulation(PolygonVertices, PlotBounds, NumPoints, TriangleNum, DiffeoParams):
	"""
	Function that visualizes the switch function corresponding to a given polygon and a given triangle number in the tree
	
	Input:
		1) PolygonVertices: Vertex Coordinates of input polygon - Nx2 numpy.array (start and end vertices must be the same)
		2) PlotBounds: Bounds for the planar plot - 4-member numpy.array ([xmin, xmax, ymin, ymax])
		3) NumPoints: Number of points for the generated grid in x and y - 2-member numpy.array ([x_resolution, y_resolution])
		4) TriangleNum: Number of triangle for which to visualize the switch function
        5) DiffeoParams: Options for the diffeomorphism construction
	
	Test:
		import numpy
		import visualization
		bounds = numpy.array([0, 5, -3, 3])
		num_points = numpy.array([101, 101])
		triangle_num = 0
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		diffeo_params = dict()
		diffeo_params['p'] = 20
		diffeo_params['epsilon'] = 1.5
		diffeo_params['varepsilon'] = 1.5
		diffeo_params['mu_1'] = 1.5
		diffeo_params['mu_2'] = 0.01
		diffeo_params['workspace'] = numpy.array([[-100,-100],[100,-100],[100,100],[-100,100],[-100,-100]])
		visualization.visualize_diffeoSwitch_triangulation(xy, bounds, num_points, triangle_num, diffeo_params)
	
	Polygon examples to test:
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		xy = numpy.array([[0,5,5,0,0,4,4,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[2.8,6.2,6.2,2.8,2.8,4.8,4.8,2.8,2.8,4.8,4.8,2.8,2.8], [1.4,1.4,10,10,8.6,8.6,7,7,4.4,4.4,2.8,2.8,1.4]]).transpose()
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		xy = numpy.array([[7,9.5,10,9,9,8,9,10,11,10,9,7], [7,8,7,6,7,7,5,5,7,9,9,7]]).transpose()
		xy = numpy.array([[7,7,8,8,7,7,10,10,9,9,10,10,7], [7,6,6,1,1,0,0,1,1,6,6,7,7]]).transpose()
		xy = numpy.array([[0,10,10,0,0,9,9,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[0,0.5,0.5,1.5,1.5,-1,-1,3,3,2,2,0,0], [0,0,-1,-1,1,1,-3,-3,1,1,-2,-2,0]]).transpose()
		xy = numpy.vstack((sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[0],sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[1])).transpose()
	"""
	# Find the tree
	tree = diffeoTreeTriangulation(PolygonVertices, DiffeoParams)

	# Generate x and y coordinates
	x_coords = numpy.linspace(PlotBounds[0], PlotBounds[1], NumPoints[0])
	y_coords = numpy.linspace(PlotBounds[2], PlotBounds[3], NumPoints[1])

	# Span all the points
	data_points = numpy.zeros((y_coords.shape[0],x_coords.shape[0]))
	for j in range(y_coords.shape[0]):
		for i in range(x_coords.shape[0]):
			candidate_point = numpy.array([[x_coords[i],y_coords[j]]])

			# Find the value of the implicit function
			sigma, sigmad, sigmadd = triangleSwitch(candidate_point, tree[TriangleNum], DiffeoParams)

			if sigma == 0. or sigma > 1:
				data_points[j][i] = numpy.nan
				continue
			else:
				# Add the data point
				data_points[j][i] = sigma
				print(i+j*y_coords.shape[0])
	
	# Plot the result
	plt.imshow(data_points, vmin=data_points[~numpy.isnan(data_points)].min(), vmax=1, origin='lower', extent=[PlotBounds[0], PlotBounds[1], PlotBounds[2], PlotBounds[3]])
	plt.axis('off')
	plt.colorbar()
	plt.show()

	return


def visualize_diffeoSwitch_convex(PolygonVertices, PlotBounds, NumPoints, PolygonNum, DiffeoParams):
	"""
	Function that visualizes the switch function corresponding to a given polygon and a given polygon number in the convex decomposition tree
	
	Input:
		1) PolygonVertices: Vertex Coordinates of input polygon - Nx2 numpy.array (start and end vertices must be the same)
		2) PlotBounds: Bounds for the planar plot - 4-member numpy.array ([xmin, xmax, ymin, ymax])
		3) NumPoints: Number of points for the generated grid in x and y - 2-member numpy.array ([x_resolution, y_resolution])
		4) PolygonNum: Number of polygon for which to visualize the switch function
        5) DiffeoParams: Options for the diffeomorphism construction
	
	Test:
		import numpy
		import visualization
		bounds = numpy.array([0, 5, -3, 3])
		num_points = numpy.array([101, 101])
		polygon_num = 0
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		diffeo_params = dict()
		diffeo_params['p'] = 20
		diffeo_params['epsilon'] = 1.5
		diffeo_params['varepsilon'] = 1.5
		diffeo_params['mu_1'] = 1.5
		diffeo_params['mu_2'] = 0.01
		diffeo_params['workspace'] = numpy.array([[-100,-100],[100,-100],[100,100],[-100,100],[-100,-100]])
		visualization.visualize_diffeoSwitch_convex(xy, bounds, num_points, polygon_num, diffeo_params)
	
	Polygon examples to test:
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		xy = numpy.array([[0,5,5,0,0,4,4,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[2.8,6.2,6.2,2.8,2.8,4.8,4.8,2.8,2.8,4.8,4.8,2.8,2.8], [1.4,1.4,10,10,8.6,8.6,7,7,4.4,4.4,2.8,2.8,1.4]]).transpose()
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		xy = numpy.array([[7,9.5,10,9,9,8,9,10,11,10,9,7], [7,8,7,6,7,7,5,5,7,9,9,7]]).transpose()
		xy = numpy.array([[7,7,8,8,7,7,10,10,9,9,10,10,7], [7,6,6,1,1,0,0,1,1,6,6,7,7]]).transpose()
		xy = numpy.array([[0,10,10,0,0,9,9,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[0,0.5,0.5,1.5,1.5,-1,-1,3,3,2,2,0,0], [0,0,-1,-1,1,1,-3,-3,1,1,-2,-2,0]]).transpose()
		xy = numpy.vstack((sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[0],sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[1])).transpose()
	"""
	# Find the tree
	tree = diffeoTreeConvex(PolygonVertices, DiffeoParams)

	# Generate x and y coordinates
	x_coords = numpy.linspace(PlotBounds[0], PlotBounds[1], NumPoints[0])
	y_coords = numpy.linspace(PlotBounds[2], PlotBounds[3], NumPoints[1])

	# Span all the points
	data_points = numpy.zeros((y_coords.shape[0],x_coords.shape[0]))
	for j in range(y_coords.shape[0]):
		for i in range(x_coords.shape[0]):
			candidate_point = numpy.array([[x_coords[i],y_coords[j]]])

			# Find the value of the implicit function
			sigma, sigmad, sigmadd = polygonSwitch(candidate_point, tree[PolygonNum], DiffeoParams)

			if sigma == 0. or sigma > 1:
				data_points[j][i] = numpy.nan
				continue
			else:
				# Add the data point
				data_points[j][i] = sigma
				print(i+j*y_coords.shape[0])
	
	# Plot the result
	plt.imshow(data_points, vmin=data_points[~numpy.isnan(data_points)].min(), vmax=1, origin='lower', extent=[PlotBounds[0], PlotBounds[1], PlotBounds[2], PlotBounds[3]])
	plt.axis('off')
	plt.colorbar()
	plt.show()

	return


def visualize_implicit(PolygonVertices, PlotBounds, NumPoints, DiffeoParams):
	"""
	Function that visualizes the implicit function corresponding to a given polygon
	
	Input:
		1) PolygonVertices: Vertex Coordinates of input polygon - Nx2 numpy.array (start and end vertices must be the same)
		2) PlotBounds: Bounds for the planar plot - 4-member numpy.array ([xmin, xmax, ymin, ymax])
		3) NumPoints: Number of points for the generated grid in x and y - 2-member numpy.array ([x_resolution, y_resolution])
        4) DiffeoParams: Options for the diffeomorphism construction
	
	Test:
		import numpy
		import visualization
		bounds = numpy.array([0, 5, -3, 3])
		num_points = numpy.array([101, 101])
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		diffeo_params = dict()
		diffeo_params['p'] = 20
		diffeo_params['epsilon'] = 1.5
		diffeo_params['varepsilon'] = 1.5
		diffeo_params['mu_1'] = 1.5
		diffeo_params['mu_2'] = 0.01
		diffeo_params['workspace'] = numpy.array([[-100,-100],[100,-100],[100,100],[-100,100],[-100,-100]])
		visualization.visualize_implicit(xy, bounds, num_points, diffeo_params)
	
	Polygon examples to test:
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		xy = numpy.array([[0,5,5,0,0,4,4,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[2.8,6.2,6.2,2.8,2.8,4.8,4.8,2.8,2.8,4.8,4.8,2.8,2.8], [1.4,1.4,10,10,8.6,8.6,7,7,4.4,4.4,2.8,2.8,1.4]]).transpose()
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		xy = numpy.array([[7,9.5,10,9,9,8,9,10,11,10,9,7], [7,8,7,6,7,7,5,5,7,9,9,7]]).transpose()
		xy = numpy.array([[7,7,8,8,7,7,10,10,9,9,10,10,7], [7,6,6,1,1,0,0,1,1,6,6,7,7]]).transpose()
		xy = numpy.array([[0,10,10,0,0,9,9,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[0,0.5,0.5,1.5,1.5,-1,-1,3,3,2,2,0,0], [0,0,-1,-1,1,1,-3,-3,1,1,-2,-2,0]]).transpose()
		xy = numpy.vstack((sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[0],sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[1])).transpose()
	"""
	# Find the tree
	tree = diffeoTreeTriangulation(PolygonVertices, DiffeoParams)

	# Generate x and y coordinates
	x_coords = numpy.linspace(PlotBounds[0], PlotBounds[1], NumPoints[0])
	y_coords = numpy.linspace(PlotBounds[2], PlotBounds[3], NumPoints[1])

	# Span all the points
	data_points = numpy.zeros((y_coords.shape[0],x_coords.shape[0]))
	current_counter = 0
	time_now = time.time()
	for j in range(y_coords.shape[0]):
		for i in range(x_coords.shape[0]):
			candidate_point = numpy.array([[x_coords[i],y_coords[j]]])

			# Find the value of the implicit function
			beta, betad, betadd = polygonImplicit(candidate_point, tree, DiffeoParams)

			if beta < 0:
				data_points[j][i] = numpy.nan
				continue
			else:
				# Add the data point
				data_points[j][i] = beta
				if i+j*y_coords.shape[0] - current_counter >= 100:
					print([i+j*y_coords.shape[0],(time.time()-time_now)/100])
					current_counter = i+j*y_coords.shape[0]
					time_now = time.time()
	
	# Plot the result
	plt.imshow(data_points, vmin=data_points[~numpy.isnan(data_points)].min(), vmax=data_points[~numpy.isnan(data_points)].max(), origin='lower', extent=[PlotBounds[0], PlotBounds[1], PlotBounds[2], PlotBounds[3]])
	plt.axis('off')
	plt.colorbar()
	plt.show()

	return


def visualize_virtualLIDAR(Polygons, RobotState, RobotRadius, DiffeoParams):
	"""
	Function that plots the robot, the obstacles and the LIDAR in the model space
	
	Input:
		1) Polygons: Vertex Coordinates of input polygons - M-member list of Nx2 numpy.array objects (start and end vertices must be the same)
		2) RobotState: State of the robot in the physical space - 3-member numpy.array (x, y, theta)
		3) RobotRadius: Robot radius (m)
        4) DiffeoParams: Options for the diffeomorphism construction
	
	Test:
		import numpy
		import visualization
		robot_state = numpy.array([0, 0, 0])
		robot_radius = 0.25
		polygon_list = []
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		polygon_list.append(xy)
		diffeo_params = dict()
		diffeo_params['p'] = 20
		diffeo_params['epsilon'] = 1.0
		diffeo_params['varepsilon'] = 1.5
		diffeo_params['mu_1'] = 1.0
		diffeo_params['mu_2'] = 0.01
		diffeo_params['workspace'] = numpy.array([[-100,-100],[100,-100],[100,100],[-100,100],[-100,-100]])
		visualization.visualize_virtualLIDAR(polygon_list, robot_state, robot_radius, diffeo_params)
	
	Polygon examples to test:
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		xy = numpy.array([[0,5,5,0,0,4,4,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[2.8,6.2,6.2,2.8,2.8,4.8,4.8,2.8,2.8,4.8,4.8,2.8,2.8], [1.4,1.4,10,10,8.6,8.6,7,7,4.4,4.4,2.8,2.8,1.4]]).transpose()
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		xy = numpy.array([[7,9.5,10,9,9,8,9,10,11,10,9,7], [7,8,7,6,7,7,5,5,7,9,9,7]]).transpose()
		xy = numpy.array([[7,7,8,8,7,7,10,10,9,9,10,10,7], [7,6,6,1,1,0,0,1,1,6,6,7,7]]).transpose()
		xy = numpy.array([[0,10,10,0,0,9,9,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[0,0.5,0.5,1.5,1.5,-1,-1,3,3,2,2,0,0], [0,0,-1,-1,1,1,-3,-3,1,1,-2,-2,0]]).transpose()
		xy = numpy.vstack((sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[0],sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[1])).transpose()
	"""
	# Construct list of polygonal objects and enlarge by robot radius
	polygon_list = []
	for i in range(len(Polygons)):
		polygon_list.append(Polygon(Polygons[i]).buffer(RobotRadius, join_style=2))
	
	# Span all the found polygons to check for intersections between the known polygons and keep only the merged polygons
	polygon_list_merged = []
	i = 0
	while (i<len(polygon_list)):
		polygon_list_merged.append(polygon_list[i])

		j = i+1
		while (j<len(polygon_list)):
			if polygon_list_merged[i].intersects(polygon_list[j]):
				polygon_list_merged[i] = polygon_list_merged[i].union(polygon_list[j])
				polygon_list_merged[i] = polygon_list_merged[i].simplify(0.08, preserve_topology=True) # simplify polygon to eliminate strange small corners
				del(polygon_list[j])
			else:
				j = j+1
		polygon_list_merged[i] = sp.geometry.polygon.orient(polygon_list_merged[i], 1.0) # orient polygon to be CCW
		i = i+1
	PolygonList = polygon_list_merged

	# Register robot state
	RobotPositionX = RobotState[0]
	RobotPositionY = RobotState[1]
	RobotOrientation = RobotState[2]
	RobotPosition = numpy.array([RobotPositionX, RobotPositionY])

	# Create fake LIDAR with range measurements
	NumSample = 101
	MinAngle = -2.35
	MaxAngle = 2.35
	Range = 4
	Infinity = 20
	Resolution = (MaxAngle - MinAngle)/(NumSample-1)
	R = Range*numpy.ones(NumSample)
	LIDAR = LIDARClass(R, Range, Infinity, MinAngle, MaxAngle, Resolution)
	
	# Complete LIDAR readings
	LIDAR = completeLIDAR2D(LIDAR)

	# Set the LIDAR rays that hit known obstacles to the LIDAR range
	for i in range(len(PolygonList)):
		LIDAR = compensateObstacleLIDAR2D(RobotState, PolygonList[i], LIDAR)
	
	# Construct list of diffeo trees for all objects
	DiffeoTreeArray = []
	for i in range(len(polygon_list_merged)):
		coords = numpy.vstack((polygon_list_merged[i].exterior.coords.xy[0],polygon_list_merged[i].exterior.coords.xy[1])).transpose()
		DiffeoTreeArray.append(diffeoTreeTriangulation(coords, DiffeoParams))
	
	# Find list of polygon objects in the model layer based on the known obstacles
	KnownObstaclesModel = []
	for i in range(len(DiffeoTreeArray)):
		theta = numpy.linspace(-numpy.pi, numpy.pi, 15)
		x_coords = DiffeoTreeArray[i][-1]['center'][0][0] + DiffeoTreeArray[i][-1]['radius']*numpy.cos(theta)
		y_coords = DiffeoTreeArray[i][-1]['center'][0][1] + DiffeoTreeArray[i][-1]['radius']*numpy.sin(theta)
		model_disk_coords = numpy.vstack((x_coords,y_coords)).transpose()
		KnownObstaclesModel.append(sp.geometry.polygon.orient(Polygon(model_disk_coords), 1.0))
	
	# Find the diffeomorphism and its jacobian at the robot position, along with the necessary second derivatives
	RobotPositionTransformed = numpy.array([RobotPosition])
	RobotPositionTransformedD = numpy.eye(2)
	RobotPositionTransformedDD = numpy.zeros(8)
	for i in range(len(DiffeoTreeArray)):
		TempPositionTransformed, TempPositionTransformedD, TempPositionTransformedDD = polygonDiffeoTriangulation(RobotPositionTransformed, DiffeoTreeArray[i], DiffeoParams)

		res1 = TempPositionTransformedD[0][0]*RobotPositionTransformedDD[0] + TempPositionTransformedD[0][1]*RobotPositionTransformedDD[4] + RobotPositionTransformedD[0][0]*(TempPositionTransformedDD[0]*RobotPositionTransformedD[0][0] + TempPositionTransformedDD[1]*RobotPositionTransformedD[1][0]) + RobotPositionTransformedD[1][0]*(TempPositionTransformedDD[2]*RobotPositionTransformedD[0][0] + TempPositionTransformedDD[3]*RobotPositionTransformedD[1][0])
		res2 = TempPositionTransformedD[0][0]*RobotPositionTransformedDD[1] + TempPositionTransformedD[0][1]*RobotPositionTransformedDD[5] + RobotPositionTransformedD[0][0]*(TempPositionTransformedDD[0]*RobotPositionTransformedD[0][1] + TempPositionTransformedDD[1]*RobotPositionTransformedD[1][1]) + RobotPositionTransformedD[1][0]*(TempPositionTransformedDD[2]*RobotPositionTransformedD[0][1] + TempPositionTransformedDD[3]*RobotPositionTransformedD[1][1])
		res3 = TempPositionTransformedD[0][0]*RobotPositionTransformedDD[2] + TempPositionTransformedD[0][1]*RobotPositionTransformedDD[6] + RobotPositionTransformedD[0][1]*(TempPositionTransformedDD[0]*RobotPositionTransformedD[0][0] + TempPositionTransformedDD[1]*RobotPositionTransformedD[1][0]) + RobotPositionTransformedD[1][1]*(TempPositionTransformedDD[2]*RobotPositionTransformedD[0][0] + TempPositionTransformedDD[3]*RobotPositionTransformedD[1][0])
		res4 = TempPositionTransformedD[0][0]*RobotPositionTransformedDD[3] + TempPositionTransformedD[0][1]*RobotPositionTransformedDD[7] + RobotPositionTransformedD[0][1]*(TempPositionTransformedDD[0]*RobotPositionTransformedD[0][1] + TempPositionTransformedDD[1]*RobotPositionTransformedD[1][1]) + RobotPositionTransformedD[1][1]*(TempPositionTransformedDD[2]*RobotPositionTransformedD[0][1] + TempPositionTransformedDD[3]*RobotPositionTransformedD[1][1])
		res5 = TempPositionTransformedD[1][0]*RobotPositionTransformedDD[0] + TempPositionTransformedD[1][1]*RobotPositionTransformedDD[4] + RobotPositionTransformedD[0][0]*(TempPositionTransformedDD[4]*RobotPositionTransformedD[0][0] + TempPositionTransformedDD[5]*RobotPositionTransformedD[1][0]) + RobotPositionTransformedD[1][0]*(TempPositionTransformedDD[6]*RobotPositionTransformedD[0][0] + TempPositionTransformedDD[7]*RobotPositionTransformedD[1][0])
		res6 = TempPositionTransformedD[1][0]*RobotPositionTransformedDD[1] + TempPositionTransformedD[1][1]*RobotPositionTransformedDD[5] + RobotPositionTransformedD[0][0]*(TempPositionTransformedDD[4]*RobotPositionTransformedD[0][1] + TempPositionTransformedDD[5]*RobotPositionTransformedD[1][1]) + RobotPositionTransformedD[1][0]*(TempPositionTransformedDD[6]*RobotPositionTransformedD[0][1] + TempPositionTransformedDD[7]*RobotPositionTransformedD[1][1])
		res7 = TempPositionTransformedD[1][0]*RobotPositionTransformedDD[2] + TempPositionTransformedD[1][1]*RobotPositionTransformedDD[6] + RobotPositionTransformedD[0][1]*(TempPositionTransformedDD[4]*RobotPositionTransformedD[0][0] + TempPositionTransformedDD[5]*RobotPositionTransformedD[1][0]) + RobotPositionTransformedD[1][1]*(TempPositionTransformedDD[6]*RobotPositionTransformedD[0][0] + TempPositionTransformedDD[7]*RobotPositionTransformedD[1][0])
		res8 = TempPositionTransformedD[1][0]*RobotPositionTransformedDD[3] + TempPositionTransformedD[1][1]*RobotPositionTransformedDD[7] + RobotPositionTransformedD[0][1]*(TempPositionTransformedDD[4]*RobotPositionTransformedD[0][1] + TempPositionTransformedDD[5]*RobotPositionTransformedD[1][1]) + RobotPositionTransformedD[1][1]*(TempPositionTransformedDD[6]*RobotPositionTransformedD[0][1] + TempPositionTransformedDD[7]*RobotPositionTransformedD[1][1])
		RobotPositionTransformedDD[0] = res1
		RobotPositionTransformedDD[1] = res2
		RobotPositionTransformedDD[2] = res3
		RobotPositionTransformedDD[3] = res4
		RobotPositionTransformedDD[4] = res5
		RobotPositionTransformedDD[5] = res6
		RobotPositionTransformedDD[6] = res7
		RobotPositionTransformedDD[7] = res8

		RobotPositionTransformedD = numpy.matmul(TempPositionTransformedD, RobotPositionTransformedD)

		RobotPositionTransformed = TempPositionTransformed
	
	# Find transformed robot orientation
	RobotOrientationTransformed = numpy.arctan2(RobotPositionTransformedD[1][0]*numpy.cos(RobotOrientation)+RobotPositionTransformedD[1][1]*numpy.sin(RobotOrientation), RobotPositionTransformedD[0][0]*numpy.cos(RobotOrientation)+RobotPositionTransformedD[0][1]*numpy.sin(RobotOrientation))

	# Find transformed robot state
	RobotStateTransformed = numpy.array([RobotPositionTransformed[0][0],RobotPositionTransformed[0][1],RobotOrientationTransformed])

	# Read LIDAR data in the model space to account for the known obstacles
	LIDARmodel = readLIDAR2D(RobotStateTransformed, KnownObstaclesModel, LIDAR.Range-numpy.linalg.norm(RobotPositionTransformed-RobotPosition), LIDAR.MinAngle, LIDAR.MaxAngle, LIDAR.NumSample)

	# Find local freespace; the robot radius can be zero because we have already dilated the obstacles
	LF_model = localfreespaceLIDAR2D(RobotStateTransformed, 0.0, LIDARmodel)

	# Plot LIDAR points
	fig, ax = plt.subplots()
	fig.set_tight_layout(True)
	lidar_plot, = ax.plot(RobotPositionTransformed[0][0] + LIDARmodel.RangeMeasurements*numpy.cos(LIDARmodel.Angle+RobotOrientationTransformed), RobotPositionTransformed[0][1] + LIDARmodel.RangeMeasurements*numpy.sin(LIDARmodel.Angle+RobotOrientationTransformed), '.r')
	ax.set_aspect('equal', 'box')

	# Plot all polygons in the physical space
	for i in range(len(polygon_list_merged)):
		coords = numpy.vstack((polygon_list_merged[i].exterior.coords.xy[0],polygon_list_merged[i].exterior.coords.xy[1])).transpose()
		pgon = plt.Polygon(coords)
		pgon.set_color('c')
		ax.add_patch(pgon)

	# Plot all transformed polygons
	for i in range(len(KnownObstaclesModel)):
		coords = numpy.vstack((KnownObstaclesModel[i].exterior.coords.xy[0],KnownObstaclesModel[i].exterior.coords.xy[1])).transpose()
		pgon = plt.Polygon(coords, alpha=0.3)
		ax.add_patch(pgon)
	
	# Robot polygon points
	bottom_left_point = numpy.array([[-0.25,-0.125]])
	bottom_right_point = numpy.array([[0.25,-0.125]])
	top_right_point = numpy.array([[0.25,0.125]])
	top_left_point = numpy.array([[-0.25,0.125]])
	
	# Plot the robot in the physical space
	RotMat = numpy.array([[numpy.cos(RobotOrientation), -numpy.sin(RobotOrientation)], [numpy.sin(RobotOrientation), numpy.cos(RobotOrientation)]])
	bottom_left_point_physical = numpy.dot(RotMat, bottom_left_point.transpose()).transpose()
	bottom_right_point_physical = numpy.dot(RotMat, bottom_right_point.transpose()).transpose()
	top_right_point_physical = numpy.dot(RotMat, top_right_point.transpose()).transpose()
	top_left_point_physical = numpy.dot(RotMat, top_left_point.transpose()).transpose()
	robot_polygon_physical = numpy.array([bottom_left_point_physical[0], bottom_right_point_physical[0], top_right_point_physical[0], top_left_point_physical[0], bottom_left_point_physical[0]]) + RobotPosition
	pgon = plt.Polygon(robot_polygon_physical, alpha=0.5)
	pgon.set_color('r')
	ax.add_patch(pgon)
	
	# Plot the robot in the model space
	RotMatTransformed = numpy.array([[numpy.cos(RobotOrientationTransformed), -numpy.sin(RobotOrientationTransformed)], [numpy.sin(RobotOrientationTransformed), numpy.cos(RobotOrientationTransformed)]])
	bottom_left_point_transformed = numpy.dot(RotMatTransformed, bottom_left_point.transpose()).transpose()
	bottom_right_point_transformed = numpy.dot(RotMatTransformed, bottom_right_point.transpose()).transpose()
	top_right_point_transformed = numpy.dot(RotMatTransformed, top_right_point.transpose()).transpose()
	top_left_point_transformed = numpy.dot(RotMatTransformed, top_left_point.transpose()).transpose()
	robot_polygon_model = numpy.array([bottom_left_point_transformed[0], bottom_right_point_transformed[0], top_right_point_transformed[0], top_left_point_transformed[0], bottom_left_point_transformed[0]]) + RobotPositionTransformed
	pgon = plt.Polygon(robot_polygon_model)
	ax.add_patch(pgon)

	# Plot the local freespace in the model space
	pgon = plt.Polygon(LF_model, alpha=0.3)
	pgon.set_color('g')
	ax.add_patch(pgon)

	plt.show()

	return


def visualize_tree_triangulation(PolygonVertices, DiffeoParams):
	"""
	Function that plots the generated tree to be used in the diffeomorphism (based on the ear clipping method)
	
	Input:
		1) PolygonVertices: Vertex Coordinates of input polygon - Nx2 numpy.array (start and end vertices must be the same)
        2) DiffeoParams: Options for the diffeomorphism construction
	
	Test:
		import numpy
		import visualization
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		diffeo_params = dict()
		diffeo_params['p'] = 20
		diffeo_params['epsilon'] = 1.0
		diffeo_params['varepsilon'] = 1.5
		diffeo_params['mu_1'] = 1.0
		diffeo_params['mu_2'] = 0.01
		diffeo_params['workspace'] = numpy.array([[-100,-100],[100,-100],[100,100],[-100,100],[-100,-100]])
		visualization.visualize_tree_triangulation(xy, diffeo_params)
	
	Polygon examples to test:
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		xy = numpy.array([[0,5,5,0,0,4,4,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[2.8,6.2,6.2,2.8,2.8,4.8,4.8,2.8,2.8,4.8,4.8,2.8,2.8], [1.4,1.4,10,10,8.6,8.6,7,7,4.4,4.4,2.8,2.8,1.4]]).transpose()
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		xy = numpy.array([[7,9.5,10,9,9,8,9,10,11,10,9,7], [7,8,7,6,7,7,5,5,7,9,9,7]]).transpose()
		xy = numpy.array([[7,7,8,8,7,7,10,10,9,9,10,10,7], [7,6,6,1,1,0,0,1,1,6,6,7,7]]).transpose()
		xy = numpy.array([[0,10,10,0,0,9,9,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[0,0.5,0.5,1.5,1.5,-1,-1,3,3,2,2,0,0], [0,0,-1,-1,1,1,-3,-3,1,1,-2,-2,0]]).transpose()
		xy = numpy.vstack((sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[0],sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[1])).transpose()
	"""
	# Transpose the input
	xy_transpose = PolygonVertices.transpose()

	# Find the dilated polygon
	polygon_dilated = Polygon(PolygonVertices).buffer(DiffeoParams['varepsilon'], join_style=1)
	polygon_dilated_vertices_transpose = numpy.vstack((polygon_dilated.exterior.coords.xy[0], polygon_dilated.exterior.coords.xy[1]))

	# Plot the initial polygon
	plt.plot(xy_transpose[0][:], xy_transpose[1][:], '-b')
	plt.axis('equal')

	# Find the tree
	tree = diffeoTreeTriangulation(PolygonVertices, DiffeoParams)

	# Plot the tree
	for i in range(0,len(tree)):
		triangle = numpy.vstack((tree[i]['vertices'],tree[i]['vertices'][0]))
		triangle = triangle.transpose()
		plt.plot(triangle[0][:], triangle[1][:], '-b')

		polygon_tilde = numpy.vstack((tree[i]['vertices_tilde'],tree[i]['vertices_tilde'][0]))
		polygon_tilde = polygon_tilde.transpose()
		plt.plot(polygon_tilde[0][:], polygon_tilde[1][:])

		plt.plot(tree[i]['center'][0][0], tree[i]['center'][0][1], 'ok')
	
	# Plot the dilated polygon
	plt.plot(polygon_dilated_vertices_transpose[0][:], polygon_dilated_vertices_transpose[1][:], '-m')
	
	plt.show()

	return tree


def visualize_tree_convex(PolygonVertices, DiffeoParams):
	"""
	Function that plots the generated tree to be used in the diffeomorphism (based on convex decomposition)
	
	Input:
		1) PolygonVertices: Vertex Coordinates of input polygon - Nx2 numpy.array (start and end vertices must be the same)
        2) DiffeoParams: Options for the diffeomorphism construction
	
	Test:
		import numpy
		import visualization
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		diffeo_params = dict()
		diffeo_params['p'] = 20
		diffeo_params['epsilon'] = 1.0
		diffeo_params['varepsilon'] = 1.5
		diffeo_params['mu_1'] = 1.0
		diffeo_params['mu_2'] = 0.01
		diffeo_params['workspace'] = numpy.array([[-100,-100],[100,-100],[100,100],[-100,100],[-100,-100]])
		visualization.visualize_tree_convex(xy, diffeo_params)
	
	Polygon examples to test:
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		xy = numpy.array([[0,5,5,0,0,4,4,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[2.8,6.2,6.2,2.8,2.8,4.8,4.8,2.8,2.8,4.8,4.8,2.8,2.8], [1.4,1.4,10,10,8.6,8.6,7,7,4.4,4.4,2.8,2.8,1.4]]).transpose()
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		xy = numpy.array([[7,9.5,10,9,9,8,9,10,11,10,9,7], [7,8,7,6,7,7,5,5,7,9,9,7]]).transpose()
		xy = numpy.array([[7,7,8,8,7,7,10,10,9,9,10,10,7], [7,6,6,1,1,0,0,1,1,6,6,7,7]]).transpose()
		xy = numpy.array([[0,10,10,0,0,9,9,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[0,0.5,0.5,1.5,1.5,-1,-1,3,3,2,2,0,0], [0,0,-1,-1,1,1,-3,-3,1,1,-2,-2,0]]).transpose()
		xy = numpy.vstack((sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[0],sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[1])).transpose()
	"""
	# Transpose the input
	xy_transpose = PolygonVertices.transpose()

	# Find the dilated polygon
	polygon_dilated = Polygon(PolygonVertices).buffer(DiffeoParams['varepsilon'], join_style=1)
	polygon_dilated_vertices_transpose = numpy.vstack((polygon_dilated.exterior.coords.xy[0], polygon_dilated.exterior.coords.xy[1]))

	# Plot the initial polygon
	plt.plot(xy_transpose[0][:], xy_transpose[1][:], '-b')
	plt.axis('equal')

	# Find the tree
	tree = diffeoTreeConvex(PolygonVertices, DiffeoParams)

	# Plot the tree
	for i in range(0,len(tree)):
		polygon = numpy.vstack((tree[i]['augmented_vertices'],tree[i]['augmented_vertices'][0]))
		polygon = polygon.transpose()
		plt.plot(polygon[0][:], polygon[1][:], '-b')

		polygon_tilde = numpy.vstack((tree[i]['vertices_tilde'],tree[i]['vertices_tilde'][0]))
		polygon_tilde = polygon_tilde.transpose()
		plt.plot(polygon_tilde[0][:], polygon_tilde[1][:])

		plt.plot(tree[i]['center'][0][0], tree[i]['center'][0][1], 'ok')
	
	# Plot the dilated polygon
	plt.plot(polygon_dilated_vertices_transpose[0][:], polygon_dilated_vertices_transpose[1][:], '-m')
	
	plt.show()

	return tree


def visualize_map_triangulation(PolygonVertices, DiffeoParams):
	"""
	Function that plots the final sphere constructed from a diffeomorphism of one polygon (based on the ear clipping method)
	
	Input:
		1) PolygonVertices: Vertex Coordinates of input polygon - Nx2 numpy.array (start and end vertices must be the same)
        2) DiffeoParams: Options for the diffeomorphism construction
	
	Test:
		import numpy
		import visualization
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		diffeo_params = dict()
		diffeo_params['p'] = 20
		diffeo_params['epsilon'] = 1.0
		diffeo_params['varepsilon'] = 1.5
		diffeo_params['mu_1'] = 1.0
		diffeo_params['mu_2'] = 0.01
		diffeo_params['workspace'] = numpy.array([[-100,-100],[100,-100],[100,100],[-100,100],[-100,-100]])
		visualization.visualize_map_triangulation(xy, diffeo_params)
	
	Polygon examples to test:
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		xy = numpy.array([[0,5,5,0,0,4,4,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[2.8,6.2,6.2,2.8,2.8,4.8,4.8,2.8,2.8,4.8,4.8,2.8,2.8], [1.4,1.4,10,10,8.6,8.6,7,7,4.4,4.4,2.8,2.8,1.4]]).transpose()
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		xy = numpy.array([[7,9.5,10,9,9,8,9,10,11,10,9,7], [7,8,7,6,7,7,5,5,7,9,9,7]]).transpose()
		xy = numpy.array([[7,7,8,8,7,7,10,10,9,9,10,10,7], [7,6,6,1,1,0,0,1,1,6,6,7,7]]).transpose()
		xy = numpy.array([[0,10,10,0,0,9,9,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[0,0.5,0.5,1.5,1.5,-1,-1,3,3,2,2,0,0], [0,0,-1,-1,1,1,-3,-3,1,1,-2,-2,0]]).transpose()
		xy = numpy.vstack((sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[0],sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[1])).transpose()
	"""
	# Transpose the input
	xy_transpose = PolygonVertices.transpose()

	# Construct the interpolation functions
	t = numpy.linspace(0, PolygonVertices.shape[0]+1, num=PolygonVertices.shape[0], endpoint=True)
	f_x = scipy.interpolate.interp1d(t, xy_transpose[0][:])
	f_y = scipy.interpolate.interp1d(t, xy_transpose[1][:])

	# Find the new points
	t_new = numpy.linspace(0, PolygonVertices.shape[0]+1, num=1001, endpoint=True)
	x_new = f_x(t_new)
	y_new = f_y(t_new)
	xy_new = numpy.vstack((x_new,y_new)).transpose()

	# Plot the initial polygon
	plt.plot(xy_transpose[0][:], xy_transpose[1][:], 'o', x_new, y_new, '-')
	plt.axis('equal')

	# Find the points on the final sphere
	tree = diffeoTreeTriangulation(PolygonVertices, DiffeoParams)
	x_deformed_array = numpy.array([[0,0]])
	for i in range(len(t_new)):
		x_deformed, x_deformedd, x_deformeddd = polygonDiffeoTriangulation(numpy.array([[x_new[i],y_new[i]]]), tree, DiffeoParams)
		x_deformed_array = numpy.vstack((x_deformed_array, x_deformed))
	
	# Plot the sphere
	x_deformed_array = x_deformed_array[1:][:]
	x_deformed_array = x_deformed_array.transpose()
	plt.plot(x_deformed_array[0][:], x_deformed_array[1][:], '')

	plt.show()

	return


def visualize_map_convex(PolygonVertices, DiffeoParams):
	"""
	Function that plots the final sphere constructed from a diffeomorphism of one polygon (based on the convex decomposition method)
	
	Input:
		1) PolygonVertices: Vertex Coordinates of input polygon - Nx2 numpy.array (start and end vertices must be the same)
        2) DiffeoParams: Options for the diffeomorphism construction
	
	Test:
		import numpy
		import visualization
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		diffeo_params = dict()
		diffeo_params['p'] = 20
		diffeo_params['epsilon'] = 1.0
		diffeo_params['varepsilon'] = 1.5
		diffeo_params['mu_1'] = 1.0
		diffeo_params['mu_2'] = 0.01
		diffeo_params['workspace'] = numpy.array([[-100,-100],[100,-100],[100,100],[-100,100],[-100,-100]])
		visualization.visualize_map_convex(xy, diffeo_params)
	
	Polygon examples to test:
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		xy = numpy.array([[0,5,5,0,0,4,4,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[2.8,6.2,6.2,2.8,2.8,4.8,4.8,2.8,2.8,4.8,4.8,2.8,2.8], [1.4,1.4,10,10,8.6,8.6,7,7,4.4,4.4,2.8,2.8,1.4]]).transpose()
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		xy = numpy.array([[7,9.5,10,9,9,8,9,10,11,10,9,7], [7,8,7,6,7,7,5,5,7,9,9,7]]).transpose()
		xy = numpy.array([[7,7,8,8,7,7,10,10,9,9,10,10,7], [7,6,6,1,1,0,0,1,1,6,6,7,7]]).transpose()
		xy = numpy.array([[0,10,10,0,0,9,9,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[0,0.5,0.5,1.5,1.5,-1,-1,3,3,2,2,0,0], [0,0,-1,-1,1,1,-3,-3,1,1,-2,-2,0]]).transpose()
		xy = numpy.vstack((sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[0],sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[1])).transpose()
	"""
	# Transpose the input
	xy_transpose = PolygonVertices.transpose()

	# Construct the interpolation functions
	t = numpy.linspace(0, PolygonVertices.shape[0]+1, num=PolygonVertices.shape[0], endpoint=True)
	f_x = scipy.interpolate.interp1d(t, xy_transpose[0][:])
	f_y = scipy.interpolate.interp1d(t, xy_transpose[1][:])

	# Find the new points
	t_new = numpy.linspace(0, PolygonVertices.shape[0]+1, num=1001, endpoint=True)
	x_new = f_x(t_new)
	y_new = f_y(t_new)
	xy_new = numpy.vstack((x_new,y_new)).transpose()

	# Plot the initial polygon
	plt.plot(xy_transpose[0][:], xy_transpose[1][:], 'o', x_new, y_new, '-')
	plt.axis('equal')

	# Find the points on the final sphere
	tree = diffeoTreeConvex(PolygonVertices, DiffeoParams)
	x_deformed_array = numpy.array([[0,0]])
	for i in range(len(t_new)):
		x_deformed, x_deformedd, x_deformeddd = polygonDiffeoConvex(numpy.array([[x_new[i],y_new[i]]]), tree, DiffeoParams)
		x_deformed_array = numpy.vstack((x_deformed_array, x_deformed))
	
	# Plot the sphere
	x_deformed_array = x_deformed_array[1:][:]
	x_deformed_array = x_deformed_array.transpose()
	plt.plot(x_deformed_array[0][:], x_deformed_array[1][:], '')

	plt.show()

	return


def visualize_purging_triangulation(PolygonVertices, DiffeoParams, FramesPerTriangle, TimeInterval, SaveOption):
	"""
	Function that shows an animation for a purging diffeomorphism of one polygon (based on the ear clipping method)
	
	Input:
		1) PolygonVertices: Vertex Coordinates of input polygon - Nx2 numpy.array (start and end vertices must be the same)
        2) DiffeoParams: Options for the diffeomorphism construction
		3) FramesPerTriangle: Frames per triangle visualization
		4) TimeInterval: Time interval in ms between each frame
		5) SaveOption: True if gif is to be saved, False otherwise
	
	Test:
		import numpy
		import visualization
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		diffeo_params = dict()
		diffeo_params['p'] = 20
		diffeo_params['epsilon'] = 1.0
		diffeo_params['varepsilon'] = 1.5
		diffeo_params['mu_1'] = 1.0
		diffeo_params['mu_2'] = 0.01
		diffeo_params['workspace'] = numpy.array([[-100,-100],[100,-100],[100,100],[-100,100],[-100,-100]])
		visualization.visualize_purging_triangulation(xy, diffeo_params, 30, 50, False)
	
	Polygon examples to test:
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		xy = numpy.array([[0,5,5,0,0,4,4,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[2.8,6.2,6.2,2.8,2.8,4.8,4.8,2.8,2.8,4.8,4.8,2.8,2.8], [1.4,1.4,10,10,8.6,8.6,7,7,4.4,4.4,2.8,2.8,1.4]]).transpose()
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		xy = numpy.array([[7,9.5,10,9,9,8,9,10,11,10,9,7], [7,8,7,6,7,7,5,5,7,9,9,7]]).transpose()
		xy = numpy.array([[7,7,8,8,7,7,10,10,9,9,10,10,7], [7,6,6,1,1,0,0,1,1,6,6,7,7]]).transpose()
		xy = numpy.array([[0,10,10,0,0,9,9,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[0,0.5,0.5,1.5,1.5,-1,-1,3,3,2,2,0,0], [0,0,-1,-1,1,1,-3,-3,1,1,-2,-2,0]]).transpose()
		xy = numpy.vstack((sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[0],sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[1])).transpose()
	"""
	# Find the polygon
	input_polygon = Polygon(PolygonVertices)

	# Transpose the input
	xy_transpose = PolygonVertices.transpose()

	# Construct the interpolation functions
	t = numpy.linspace(0, PolygonVertices.shape[0]+1, num=PolygonVertices.shape[0], endpoint=True)
	f_x = scipy.interpolate.interp1d(t, xy_transpose[0][:])
	f_y = scipy.interpolate.interp1d(t, xy_transpose[1][:])

	# Find the new points
	t_new = numpy.linspace(0, PolygonVertices.shape[0]+1, num=1001, endpoint=True)
	x_new = f_x(t_new)
	y_new = f_y(t_new)
	xy_new = numpy.vstack((x_new,y_new)).transpose()

	# Find the diffeomorphism tree
	tree = diffeoTreeTriangulation(PolygonVertices, DiffeoParams)

	# # Find grid points
	# nrows = 1 + int(numpy.ceil((numpy.max(y_new)+0.5)-(numpy.min(y_new)-0.5))/0.5)
	# ncols = 1 + int(numpy.ceil((numpy.max(x_new)+0.5)-(numpy.min(x_new)-0.5))/0.5)
	# discretization_rows = ((numpy.max(y_new)+0.5)-(numpy.min(y_new)-0.5))/(nrows-1)
	# discretization_cols = ((numpy.max(x_new)+0.5)-(numpy.min(x_new)-0.5))/(ncols-1)
	# num_points = 1000
	# x_grid = []
	# y_grid = []

	# for row in range(0,nrows):
	# 	result_x = numpy.array([])
	# 	result_y = numpy.array([])
	# 	for point_index in range(0,num_points):
	# 		x_point = (numpy.min(x_new)-0.5) + (point_index/num_points)*((numpy.max(x_new)+0.5)-(numpy.min(x_new)-0.5))
	# 		y_point = (numpy.min(y_new)-0.5) + row*discretization_rows
	# 		if not input_polygon.contains(Point(x_point,y_point)):
	# 			result_x = numpy.hstack((result_x,x_point))
	# 			result_y = numpy.hstack((result_y,y_point))
	# 		else:
	# 			x_grid.append(result_x)
	# 			y_grid.append(result_y)
	# 			result_x = numpy.array([])
	# 			result_y = numpy.array([])
	# 	x_grid.append(result_x)
	# 	y_grid.append(result_y)
	
	# for col in range(0,ncols):
	# 	result_x = numpy.array([])
	# 	result_y = numpy.array([])
	# 	for point_index in range(0,num_points):
	# 		x_point = (numpy.min(x_new)-0.5) + col*discretization_cols
	# 		y_point = (numpy.min(y_new)-0.5) + (point_index/num_points)*((numpy.max(y_new)+0.5)-(numpy.min(y_new)-0.5))
	# 		if not input_polygon.contains(Point(x_point,y_point)):
	# 			result_x = numpy.hstack((result_x,x_point))
	# 			result_y = numpy.hstack((result_y,y_point))
	# 		else:
	# 			x_grid.append(result_x)
	# 			y_grid.append(result_y)
	# 			result_x = numpy.array([])
	# 			result_y = numpy.array([])
	# 	x_grid.append(result_x)
	# 	y_grid.append(result_y)
		
	# Initialize plot
	fig, ax = plt.subplots()
	fig.set_tight_layout(True)
	xy_plot, = ax.plot(x_new, y_new, '-', linewidth = 2)
	# xy_grid = [None]*len(x_grid)
	# for i in range(0,len(x_grid)):
	# 	xy_grid[i], = ax.plot(x_grid[i], y_grid[i], '-')
	# 	xy_grid[i].set_color('gray')
	ax.axis([numpy.min(x_new)-0.5,numpy.max(x_new)+0.5,numpy.min(y_new)-0.5,numpy.max(y_new)+0.5])
	ax.set_aspect('equal', 'box')
	ax.set_yticks([])
	ax.set_xticks([])

	# Iterate through the tree for the polygon points
	x_to_animate = numpy.array([x_new])
	y_to_animate = numpy.array([y_new])
	for j in range(len(tree)):
		x_deformed_array = numpy.array([[0,0]])
		for k in range(0,len(t_new)):
			x_deformed, x_deformedd, x_deformeddd = triangleDiffeo(numpy.array([[x_to_animate[-1][k],y_to_animate[-1][k]]]), tree[j], DiffeoParams)
			x_deformed_array = numpy.vstack((x_deformed_array, x_deformed))
		x_deformed_array = x_deformed_array[1:][:].transpose()
		x_to_animate = numpy.vstack((x_to_animate,x_deformed_array[0][:]))
		y_to_animate = numpy.vstack((y_to_animate,x_deformed_array[1][:]))
	x_to_animate = numpy.vstack((x_to_animate,x_deformed_array[0][:]))
	y_to_animate = numpy.vstack((y_to_animate,x_deformed_array[1][:]))

	# # Iterate through the tree for the grid points
	# x_to_animate_grid = [None]*len(x_grid)
	# y_to_animate_grid = [None]*len(y_grid)

	# for i in range(0,len(x_grid)):
	# 	x_to_animate_grid[i] = numpy.array([x_grid[i]])
	# 	y_to_animate_grid[i] = numpy.array([y_grid[i]])
	# 	for j in range(len(tree)):
	# 		x_deformed_array_grid = numpy.array([[0,0]])
	# 		for k in range(0,len(x_grid[i])):
	# 			x_deformed_grid, x_deformedd_gri, x_deformeddd_grid = triangleDiffeo(numpy.array([[x_to_animate_grid[i][-1][k],y_to_animate_grid[i][-1][k]]]), tree[j], DiffeoParams)
	# 			x_deformed_array_grid = numpy.vstack((x_deformed_array_grid, x_deformed_grid))
	# 		x_deformed_array_grid = x_deformed_array_grid[1:][:].transpose()
	# 		x_to_animate_grid[i] = numpy.vstack((x_to_animate_grid[i],x_deformed_array_grid[0][:]))
	# 		y_to_animate_grid[i] = numpy.vstack((y_to_animate_grid[i],x_deformed_array_grid[1][:]))
	# 	x_to_animate_grid[i] = numpy.vstack((x_to_animate_grid[i],x_deformed_array_grid[0][:]))
	# 	y_to_animate_grid[i] = numpy.vstack((y_to_animate_grid[i],x_deformed_array_grid[1][:]))


	def update(i):
		
		# Find the index to consider
		index_to_consider = int(numpy.floor(i/FramesPerTriangle))

		label = 'Purging triangle {0}'.format(index_to_consider+1)

		# Animate data
		xy_plot.set_xdata(((i%FramesPerTriangle)/FramesPerTriangle)*x_to_animate[index_to_consider+1][:] + (1-((i%FramesPerTriangle)/FramesPerTriangle))*x_to_animate[index_to_consider][:])
		xy_plot.set_ydata(((i%FramesPerTriangle)/FramesPerTriangle)*y_to_animate[index_to_consider+1][:] + (1-((i%FramesPerTriangle)/FramesPerTriangle))*y_to_animate[index_to_consider][:])
		# for j in range(0,len(x_grid)):
		# 	xy_grid[j].set_xdata(((i%FramesPerTriangle)/FramesPerTriangle)*x_to_animate_grid[j][index_to_consider+1][:] + (1-((i%FramesPerTriangle)/FramesPerTriangle))*x_to_animate_grid[j][index_to_consider][:])
		# 	xy_grid[j].set_ydata(((i%FramesPerTriangle)/FramesPerTriangle)*y_to_animate_grid[j][index_to_consider+1][:] + (1-((i%FramesPerTriangle)/FramesPerTriangle))*y_to_animate_grid[j][index_to_consider][:])
		ax.set_xlabel(label)
		if SaveOption == True:
			fig.savefig('./../../data/visualizations/figure_' + str(i) + '.pdf', bbox_inches='tight')
		return xy_plot, ax
	
	# Initialize FuncAnimation object
	anim = FuncAnimation(fig, update, frames=numpy.arange(0, 1+FramesPerTriangle*len(tree)), interval=TimeInterval)
	if SaveOption == True:
		anim.save('./../../data/visualizations/diffeomorphism_triangulation.gif', dpi=80, writer='imagemagick')
	else:
		plt.show()
	
	return


def visualize_purging_convex(PolygonVertices, DiffeoParams, FramesPerPolygon, TimeInterval, SaveOption):
	"""
	Function that shows an animation for a purging diffeomorphism of one polygon (based on convex decomposition)
	
	Input:
		1) PolygonVertices: Vertex Coordinates of input polygon - Nx2 numpy.array (start and end vertices must be the same)
        2) DiffeoParams: Options for the diffeomorphism construction
		3) FramesPerPolygon: Frames per triangle visualization
		4) TimeInterval: Time interval in ms between each frame
		5) SaveOption: True if gif is to be saved, False otherwise
	
	Test:
		import numpy
		import visualization
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		diffeo_params = dict()
		diffeo_params['p'] = 20
		diffeo_params['epsilon'] = 1.0
		diffeo_params['varepsilon'] = 1.5
		diffeo_params['mu_1'] = 1.0
		diffeo_params['mu_2'] = 0.01
		diffeo_params['workspace'] = numpy.array([[-100,-100],[100,-100],[100,100],[-100,100],[-100,-100]])
		visualization.visualize_purging_convex(xy, diffeo_params, 30, 50, False)
	
	Polygon examples to test:
		xy = numpy.array([[2.518,1.83,2.043,2.406,2.655,2.518], [0.5048,0.2963,-0.2348,-0.8039,-0.0533,0.5048]]).transpose()
		xy = numpy.array([[0,5,5,0,0,4,4,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[2.8,6.2,6.2,2.8,2.8,4.8,4.8,2.8,2.8,4.8,4.8,2.8,2.8], [1.4,1.4,10,10,8.6,8.6,7,7,4.4,4.4,2.8,2.8,1.4]]).transpose()
		xy = numpy.array([[7,6,4,5.4,5.1,7,8.9,8.6,10,8,7], [9.5,7.6,7.2,5.6,3.5,4.4,3.5,5.6,7.2,7.6,9.5]]).transpose()
		xy = numpy.array([[7,9.5,10,9,9,8,9,10,11,10,9,7], [7,8,7,6,7,7,5,5,7,9,9,7]]).transpose()
		xy = numpy.array([[7,7,8,8,7,7,10,10,9,9,10,10,7], [7,6,6,1,1,0,0,1,1,6,6,7,7]]).transpose()
		xy = numpy.array([[0,10,10,0,0,9,9,0,0], [0,0,5,5,4,4,1,1,0]]).transpose()
		xy = numpy.array([[0,0.5,0.5,1.5,1.5,-1,-1,3,3,2,2,0,0], [0,0,-1,-1,1,1,-3,-3,1,1,-2,-2,0]]).transpose()
		xy = numpy.vstack((sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[0],sp.geometry.polygon.orient(LineString(numpy.array([[0,0],[0,1],[-1,1],[-1,-1],[1,-1],[1,2],[-2,2],[-2,-2],[2,-2],[2,3],[-3,3],[-3,-3],[3,-3],[3,4],[-3,4]])).buffer(0.2).simplify(0.05),1.0).exterior.coords.xy[1])).transpose()
	"""
	# Find the polygon
	input_polygon = Polygon(PolygonVertices)

	# Transpose the input
	xy_transpose = PolygonVertices.transpose()

	# Construct the interpolation functions
	t = numpy.linspace(0, PolygonVertices.shape[0]+1, num=PolygonVertices.shape[0], endpoint=True)
	f_x = scipy.interpolate.interp1d(t, xy_transpose[0][:])
	f_y = scipy.interpolate.interp1d(t, xy_transpose[1][:])

	# Find the new points
	t_new = numpy.linspace(0, PolygonVertices.shape[0]+1, num=1001, endpoint=True)
	x_new = f_x(t_new)
	y_new = f_y(t_new)
	xy_new = numpy.vstack((x_new,y_new)).transpose()

	# Find the diffeomorphism tree
	tree = diffeoTreeConvex(PolygonVertices, DiffeoParams)

	# # Find grid points
	# nrows = 1 + int(numpy.ceil((numpy.max(y_new)+0.5)-(numpy.min(y_new)-0.5))/0.5)
	# ncols = 1 + int(numpy.ceil((numpy.max(x_new)+0.5)-(numpy.min(x_new)-0.5))/0.5)
	# discretization_rows = ((numpy.max(y_new)+0.5)-(numpy.min(y_new)-0.5))/(nrows-1)
	# discretization_cols = ((numpy.max(x_new)+0.5)-(numpy.min(x_new)-0.5))/(ncols-1)
	# num_points = 1000
	# x_grid = []
	# y_grid = []

	# for row in range(0,nrows):
	# 	result_x = numpy.array([])
	# 	result_y = numpy.array([])
	# 	for point_index in range(0,num_points):
	# 		x_point = (numpy.min(x_new)-0.5) + (point_index/num_points)*((numpy.max(x_new)+0.5)-(numpy.min(x_new)-0.5))
	# 		y_point = (numpy.min(y_new)-0.5) + row*discretization_rows
	# 		if not input_polygon.contains(Point(x_point,y_point)):
	# 			result_x = numpy.hstack((result_x,x_point))
	# 			result_y = numpy.hstack((result_y,y_point))
	# 		else:
	# 			x_grid.append(result_x)
	# 			y_grid.append(result_y)
	# 			result_x = numpy.array([])
	# 			result_y = numpy.array([])
	# 	x_grid.append(result_x)
	# 	y_grid.append(result_y)
	
	# for col in range(0,ncols):
	# 	result_x = numpy.array([])
	# 	result_y = numpy.array([])
	# 	for point_index in range(0,num_points):
	# 		x_point = (numpy.min(x_new)-0.5) + col*discretization_cols
	# 		y_point = (numpy.min(y_new)-0.5) + (point_index/num_points)*((numpy.max(y_new)+0.5)-(numpy.min(y_new)-0.5))
	# 		if not input_polygon.contains(Point(x_point,y_point)):
	# 			result_x = numpy.hstack((result_x,x_point))
	# 			result_y = numpy.hstack((result_y,y_point))
	# 		else:
	# 			x_grid.append(result_x)
	# 			y_grid.append(result_y)
	# 			result_x = numpy.array([])
	# 			result_y = numpy.array([])
	# 	x_grid.append(result_x)
	# 	y_grid.append(result_y)
		
	# Initialize plot
	fig, ax = plt.subplots()
	fig.set_tight_layout(True)
	xy_plot, = ax.plot(x_new, y_new, '-', linewidth = 2)
	# xy_grid = [None]*len(x_grid)
	# for i in range(0,len(x_grid)):
	# 	xy_grid[i], = ax.plot(x_grid[i], y_grid[i], '-')
	# 	xy_grid[i].set_color('gray')
	ax.axis([numpy.min(x_new)-0.5,numpy.max(x_new)+0.5,numpy.min(y_new)-0.5,numpy.max(y_new)+0.5])
	ax.set_aspect('equal', 'box')
	ax.set_yticks([])
	ax.set_xticks([])

	# Iterate through the tree for the polygon points
	x_to_animate = numpy.array([x_new])
	y_to_animate = numpy.array([y_new])
	for j in range(len(tree)):
		x_deformed_array = numpy.array([[0,0]])
		for k in range(0,len(t_new)):
			x_deformed, x_deformedd, x_deformeddd = polygonDiffeo(numpy.array([[x_to_animate[-1][k],y_to_animate[-1][k]]]), tree[j], DiffeoParams)
			x_deformed_array = numpy.vstack((x_deformed_array, x_deformed))
		x_deformed_array = x_deformed_array[1:][:].transpose()
		x_to_animate = numpy.vstack((x_to_animate,x_deformed_array[0][:]))
		y_to_animate = numpy.vstack((y_to_animate,x_deformed_array[1][:]))
	x_to_animate = numpy.vstack((x_to_animate,x_deformed_array[0][:]))
	y_to_animate = numpy.vstack((y_to_animate,x_deformed_array[1][:]))

	# # Iterate through the tree for the grid points
	# x_to_animate_grid = [None]*len(x_grid)
	# y_to_animate_grid = [None]*len(y_grid)

	# for i in range(0,len(x_grid)):
	# 	x_to_animate_grid[i] = numpy.array([x_grid[i]])
	# 	y_to_animate_grid[i] = numpy.array([y_grid[i]])
	# 	for j in range(len(tree)):
	# 		x_deformed_array_grid = numpy.array([[0,0]])
	# 		for k in range(0,len(x_grid[i])):
	# 			x_deformed_grid, x_deformedd_gri, x_deformeddd_grid = triangleDiffeo(numpy.array([[x_to_animate_grid[i][-1][k],y_to_animate_grid[i][-1][k]]]), tree[j], DiffeoParams)
	# 			x_deformed_array_grid = numpy.vstack((x_deformed_array_grid, x_deformed_grid))
	# 		x_deformed_array_grid = x_deformed_array_grid[1:][:].transpose()
	# 		x_to_animate_grid[i] = numpy.vstack((x_to_animate_grid[i],x_deformed_array_grid[0][:]))
	# 		y_to_animate_grid[i] = numpy.vstack((y_to_animate_grid[i],x_deformed_array_grid[1][:]))
	# 	x_to_animate_grid[i] = numpy.vstack((x_to_animate_grid[i],x_deformed_array_grid[0][:]))
	# 	y_to_animate_grid[i] = numpy.vstack((y_to_animate_grid[i],x_deformed_array_grid[1][:]))


	def update(i):
		
		# Find the index to consider
		index_to_consider = int(numpy.floor(i/FramesPerPolygon))

		label = 'Purging polygon {0}'.format(index_to_consider+1)

		# Animate data
		xy_plot.set_xdata(((i%FramesPerPolygon)/FramesPerPolygon)*x_to_animate[index_to_consider+1][:] + (1-((i%FramesPerPolygon)/FramesPerPolygon))*x_to_animate[index_to_consider][:])
		xy_plot.set_ydata(((i%FramesPerPolygon)/FramesPerPolygon)*y_to_animate[index_to_consider+1][:] + (1-((i%FramesPerPolygon)/FramesPerPolygon))*y_to_animate[index_to_consider][:])
		# for j in range(0,len(x_grid)):
		# 	xy_grid[j].set_xdata(((i%FramesPerPolygon)/FramesPerPolygon)*x_to_animate_grid[j][index_to_consider+1][:] + (1-((i%FramesPerPolygon)/FramesPerPolygon))*x_to_animate_grid[j][index_to_consider][:])
		# 	xy_grid[j].set_ydata(((i%FramesPerPolygon)/FramesPerPolygon)*y_to_animate_grid[j][index_to_consider+1][:] + (1-((i%FramesPerPolygon)/FramesPerPolygon))*y_to_animate_grid[j][index_to_consider][:])
		ax.set_xlabel(label)
		if SaveOption == True:
			fig.savefig('./../../data/visualizations/figure_' + str(i) + '.pdf', bbox_inches='tight')
		return xy_plot, ax
	
	# Initialize FuncAnimation object
	anim = FuncAnimation(fig, update, frames=numpy.arange(0, 1+FramesPerPolygon*len(tree)), interval=TimeInterval)
	if SaveOption == True:
		anim.save('./../../data/visualizations/diffeomorphism_convex.gif', dpi=80, writer='imagemagick')
	else:
		plt.show()
	
	return