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

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/foreach.hpp>

// CGAL imports
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>

// Import ear clipping method
#include <earcut.hpp>

// ROS imports
#include <ros/ros.h>

// Other imports
#include <cmath>
#include <vector>
#include <iostream>
#include <functional>
#include <algorithm>

// Define namespaces
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// Define various geometries
using point = bg::model::point<double, 2, bg::cs::cartesian>;
using polygon = bg::model::polygon<point, false, true>;
using line = bg::model::linestring<point>;
using multi_point = bg::model::multi_point<point>;
using multi_polygon = bg::model::multi_polygon<polygon>;

// Define CGAL types
typedef CGAL::Exact_predicates_inexact_constructions_kernel         CGAL_K;
typedef CGAL::Partition_traits_2<CGAL_K>                            CGAL_Traits;
typedef CGAL_Traits::Point_2                                        CGAL_Point_2;
typedef CGAL_Traits::Polygon_2                                      CGAL_Polygon_2;
typedef CGAL_Polygon_2::Vertex_iterator                             CGAL_Vertex_iterator;
typedef std::list<CGAL_Polygon_2>                                   CGAL_Polygon_list;
typedef CGAL::Creator_uniform_2<int, CGAL_Point_2>                  CGAL_Creator;
typedef CGAL::Random_points_in_square_2<CGAL_Point_2, CGAL_Creator> CGAL_Point_generator;


std::vector<std::vector<double>> MatrixMatrixMultiplication(std::vector<std::vector<double>> Matrix1, std::vector<std::vector<double>> Matrix2);


std::vector<double> MatrixVectorMultiplication(std::vector<std::vector<double>> Matrix, std::vector<double> Vector);


std::vector<std::vector<double>> VectorOuterProduct(std::vector<double> Vector1, std::vector<double> Vector2);


double MatrixDeterminant(std::vector<std::vector<double>> Matrix);


std::vector<point> StdToBoostPoint(std::vector<std::vector<double>> input);


std::vector<std::vector<double>> BoostPointToStd(std::vector<point> input);


std::vector<point> BoostPolyToBoostPoint(polygon input);


std::vector<point> BoostLineToBoostPoint(line input);


polygon BoostPointToBoostPoly(std::vector<point> input);


line BoostPointToBoostLine(std::vector<point> input);


class TriangleClass {
    /** 
     * Class that describes the properties of a single triangle in the triangulation
     * 
     * Properties:
     *  1) radius: Radius of the final sphere to be constructed (0 unless the triangle is a root triangle to be deformed to disk)
     *  2) r_center_t: The tangents from vertices 0 and 1 to the center
     *  3) r_center_n: The normals corresponding to r_center_t
     *  4) vertices: The vertices of the triangle - 3-element vector of points
     *  5) vertices_tilde: The vertices of the polygonal collar that encompasses the triangle - M-element vector of points in CCW order starting from the center in the parent
     *  6) r_t: The unit tangents for the triangle to be deformed - 3-element vector of points in CCW order (the 1st row is the shared tangent between the parent and the child)
     *  7) r_n: The unit normals for the triangle to be deformed corresponding to r_t
     *  8) r_tilde_t: The unit tangents for the polygonal collar in CCW order
     *  9) r_tilde_n: The unit normals corresponding to r_tilde_t
     *  10) center: The center in the parent used for the purging transformation, or the center of the root used for the final transformation
     *  11) adj_edge: Adjacency edge to the parent
     *  12) depth: Depth in the transformation
     *  13) predecessor: Predecessor in the purging transformation
     *  14) index: Index of the tree
     */
    public:
        double get_radius() const {return this->radius;}
        std::vector<point> get_r_center_t() const {return this->r_center_t;}
        std::vector<point> get_r_center_n() const {return this->r_center_n;}
        std::vector<point> get_vertices() const {return this->vertices;}
        std::vector<point> get_vertices_tilde() const {return this->vertices_tilde;}
        std::vector<point> get_r_t() const {return this->r_t;}
        std::vector<point> get_r_n() const {return this->r_n;}
        std::vector<point> get_r_tilde_t() const {return this->r_tilde_t;}
        std::vector<point> get_r_tilde_n() const {return this->r_tilde_n;}
        point get_center() const {return this->center;}
        std::vector<point> get_adj_edge() const {return this->adj_edge;}
        uint16_t get_depth() const {return this->depth;}
        int get_predecessor() const {return this->predecessor;}
        uint16_t get_index() const {return this->index;}

        void set_radius(double radius_in) {this->radius = radius_in;}
        void set_r_center_t(std::vector<point> r_center_t_in) {this->r_center_t = r_center_t_in;}
        void set_r_center_n(std::vector<point> r_center_n_in) {this->r_center_n = r_center_n_in;}
        void set_vertices(std::vector<point> vertices_in) {this->vertices = vertices_in;}
        void set_vertices_tilde(std::vector<point> vertices_tilde_in) {this->vertices_tilde = vertices_tilde_in;}
        void set_r_t(std::vector<point> r_t_in) {this->r_t = r_t_in;}
        void set_r_n(std::vector<point> r_n_in) {this->r_n = r_n_in;}
        void set_r_tilde_t(std::vector<point> r_tilde_t_in) {this->r_tilde_t = r_tilde_t_in;}
        void set_r_tilde_n(std::vector<point> r_tilde_n_in) {this->r_tilde_n = r_tilde_n_in;}
        void set_center(point center_in) {this->center = center_in;}
        void set_adj_edge(std::vector<point> adj_edge_in) {this->adj_edge = adj_edge_in;}
        void set_depth(uint16_t depth_in) {this->depth = depth_in;}
        void set_predecessor(int predecessor_in) {this->predecessor = predecessor_in;}
        void set_index(uint16_t index_in) {this->index = index_in;}
    
    private:
        double radius = 0.0;
        std::vector<point> r_center_t;
        std::vector<point> r_center_n;
        std::vector<point> vertices;
        std::vector<point> vertices_tilde;
        std::vector<point> r_t;
        std::vector<point> r_n;
        std::vector<point> r_tilde_t;
        std::vector<point> r_tilde_n;
        point center;
        std::vector<point> adj_edge;
        uint16_t depth;
        int predecessor;
        uint16_t index;
};

class PolygonClass {
    /** 
     * Class that describes the properties of a single polygon in the convex decomposition tree
     * 
     * Properties:
     *  1) radius: Radius of the final sphere to be constructed (0 unless the polygon is a root polygon to be deformed to disk)
     *  2) r_center_t: The tangents from vertices 0 and 1 to the center
     *  3) r_center_n: The normals corresponding to r_center_t
     *  4) vertices: The vertices of the polygon without the center of transformation - N-element vector of points
     *  5) augmented_vertices: The vertices of the interior polygon that includes the center of transformation (2nd element of the vector) - (N+1)-element vector of points
     *  6) vertices_tilde: The vertices of the polygonal collar that encompasses the polygon - M-element vector of points in CCW order
     *  7) r_t: The unit tangents for the polygon to be deformed corresponding to augmented_vertices - (N+1)-element vector of points in CCW order (the 1st row is essentially the first element of r_center_t)
     *  8) r_n: The unit normals for the polygon to be deformed corresponding to r_t
     *  9) r_tilde_t: The unit tangents for the polygonal collar in CCW order
     *  10) r_tilde_n: The unit normals corresponding to r_tilde_t
     *  11) center: The center in the parent used for the purging transformation, or the center of the root used for the final transformation
     *  12) adj_edge: Adjacency edge to the parent (empty if the polygon is the root polygon)
     *  13) depth: Depth in the transformation
     *  14) predecessor: Predecessor in the purging transformation
     *  15) index: Index of the node in the tree
     */
    public:
        double get_radius() const {return this->radius;}
        std::vector<point> get_r_center_t() const {return this->r_center_t;}
        std::vector<point> get_r_center_n() const {return this->r_center_n;}
        std::vector<point> get_vertices() const {return this->vertices;}
        std::vector<point> get_augmented_vertices() const {return this->augmented_vertices;}
        std::vector<point> get_vertices_tilde() const {return this->vertices_tilde;}
        std::vector<point> get_r_t() const {return this->r_t;}
        std::vector<point> get_r_n() const {return this->r_n;}
        std::vector<point> get_r_tilde_t() const {return this->r_tilde_t;}
        std::vector<point> get_r_tilde_n() const {return this->r_tilde_n;}
        point get_center() const {return this->center;}
        std::vector<point> get_adj_edge() const {return this->adj_edge;}
        uint16_t get_depth() const {return this->depth;}
        int get_predecessor() const {return this->predecessor;}
        uint16_t get_index() const {return this->index;}

        void set_radius(double radius_in) {this->radius = radius_in;}
        void set_r_center_t(std::vector<point> r_center_t_in) {this->r_center_t = r_center_t_in;}
        void set_r_center_n(std::vector<point> r_center_n_in) {this->r_center_n = r_center_n_in;}
        void set_vertices(std::vector<point> vertices_in) {this->vertices = vertices_in;}
        void set_augmented_vertices(std::vector<point> augmented_vertices_in) {this->augmented_vertices = augmented_vertices_in;}
        void set_vertices_tilde(std::vector<point> vertices_tilde_in) {this->vertices_tilde = vertices_tilde_in;}
        void set_r_t(std::vector<point> r_t_in) {this->r_t = r_t_in;}
        void set_r_n(std::vector<point> r_n_in) {this->r_n = r_n_in;}
        void set_r_tilde_t(std::vector<point> r_tilde_t_in) {this->r_tilde_t = r_tilde_t_in;}
        void set_r_tilde_n(std::vector<point> r_tilde_n_in) {this->r_tilde_n = r_tilde_n_in;}
        void set_center(point center_in) {this->center = center_in;}
        void set_adj_edge(std::vector<point> adj_edge_in) {this->adj_edge = adj_edge_in;}
        void set_depth(uint16_t depth_in) {this->depth = depth_in;}
        void set_predecessor(int predecessor_in) {this->predecessor = predecessor_in;}
        void set_index(uint16_t index_in) {this->index = index_in;}
    
    private:
        double radius = 0.0;
        std::vector<point> r_center_t;
        std::vector<point> r_center_n;
        std::vector<point> vertices;
        std::vector<point> augmented_vertices;
        std::vector<point> vertices_tilde;
        std::vector<point> r_t;
        std::vector<point> r_n;
        std::vector<point> r_tilde_t;
        std::vector<point> r_tilde_n;
        point center;
        std::vector<point> adj_edge;
        uint16_t depth;
        int predecessor;
        uint16_t index;
};


struct ProjectionResultStruct {
    /** 
     * Struct that contains the output of a projection result (i.e., point and distance)
     * 
     * Properties:
     *  1) projected_point: Result of the projection as a point on the projected set
     *  2) dist: Distance to the projected set
     */
    public:
        point projected_point;
        double dist;
};


polygon cvxpolyxhplane(polygon xy, point m, point n);


line polyxline(polygon xy, point m, point n);


point polyxray(polygon xy, point b, point v);


ProjectionResultStruct polydist(polygon xy, point p);


ProjectionResultStruct linedist(line xy, point p);


void polytriangulation(std::vector<std::vector<double>> xy, std::vector<std::vector<double>> workspace, bool touching_boundary, std::vector<TriangleClass> *tree);

void polyconvexdecomposition(std::vector<std::vector<double>> xy, std::vector<std::vector<double>> workspace, bool touching_boundary, std::vector<PolygonClass> *tree);