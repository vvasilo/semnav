// MIT License (modified)

// Copyright (c) 2019 The Trustees of the University of Pennsylvania
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

#include <polygeom_lib.h>


std::vector<std::vector<double>> MatrixMatrixMultiplication(std::vector<std::vector<double>> Matrix1, std::vector<std::vector<double>> Matrix2) {
    std::vector<std::vector<double>> MatrixOut = {{0.0, 0.0}, {0.0, 0.0}};
    for (size_t i = 0; i < Matrix1.size(); i++) {
        for (size_t j = 0; j < Matrix2.size(); j++) {
            for (size_t k = 0; k < Matrix1.size(); k++) {
                MatrixOut[i][j] += Matrix1[i][k] * Matrix2[k][j];
            }
        }
    }
    return MatrixOut;
}


std::vector<double> MatrixVectorMultiplication(std::vector<std::vector<double>> Matrix, std::vector<double> Vector) {
    std::vector<double> VectorOut = {0.0, 0.0};
    for (size_t i = 0; i < Matrix.size(); i++) {
        for (size_t j = 0; j < Vector.size(); j++) {
            VectorOut[i] += Matrix[i][j] * Vector[j];
        }
    }
    return VectorOut;
}


std::vector<std::vector<double>> VectorOuterProduct(std::vector<double> Vector1, std::vector<double> Vector2) {
    std::vector<std::vector<double>> MatrixOut = {{0.0, 0.0}, {0.0, 0.0}};
    for (size_t i = 0; i < Vector1.size(); i++) {
        for (size_t j = 0; j < Vector2.size(); j++) {
            MatrixOut[i][j] += Vector1[i] * Vector2[j];
        }
    }
    return MatrixOut;
}


double MatrixDeterminant(std::vector<std::vector<double>> Matrix) {
    return (Matrix[0][0]*Matrix[1][1]-Matrix[0][1]*Matrix[1][0]);
}


std::vector<point> StdToBoostPoint(std::vector<std::vector<double>> input) {
    /**
     * Function that takes as input an array of points in std::vector<std::vector<double>> format and converts them to std::vector<point>
     * 
     * Input:
     *  1) input: Array of points
     * 
     * Output:
     *  2) output: Output of points in Boost format
     */
    std::vector<point> output;
    for (size_t i = 0; i < input.size(); i++) {
        output.push_back(point(input[i][0], input[i][1]));
    }

    return output;
}


std::vector<std::vector<double>> BoostPointToStd(std::vector<point> input) {
    /**
     * Function that takes as input an array of points in std::vector<point> format and converts them to std::vector<std::vector<double>>
     * 
     * Input:
     *  1) input: Array of points
     * 
     * Output:
     *  2) output: Output of points in std format
     */
    std::vector<std::vector<double>> output;
    for (size_t i = 0; i < input.size(); i++) {
        output.push_back({input[i].get<0>(), input[i].get<1>()});
    }

    return output;
}


std::vector<point> BoostPolyToBoostPoint(polygon input) {
    /**
     * Function that takes as input a polygon in Boost format and converts it to std::vector<point>
     * 
     * Input:
     *  1) input: Polygon in Boost format
     * 
     * Output:
     *  2) output: Vector of points in Boost format
     */
    std::vector<point> output;
    for (auto it = boost::begin(bg::exterior_ring(input)); it != boost::end(bg::exterior_ring(input)); ++it) {
        point new_point = point(bg::get<0>(*it), bg::get<1>(*it));
        output.push_back(new_point);
    }

    return output;
}


std::vector<point> BoostLineToBoostPoint(line input) {
    /**
     * Function that takes as input a line in Boost format and converts it to std::vector<point>
     * 
     * Input:
     *  1) input: Line in Boost format
     * 
     * Output:
     *  2) output: Vector of points in Boost format
     */
    std::vector<point> output;
    for (auto it = input.begin(); it != input.end(); ++it) {
        point new_point = point(bg::get<0>(*it), bg::get<1>(*it));
        output.push_back(new_point);
    }

    return output;
}


polygon BoostPointToBoostPoly(std::vector<point> input) {
    /**
     * Function that takes as input a std::vector<point> series of points and converts them to Boost polygon format
     * 
     * Input:
     *  1) input: Vector of points in Boost format
     * 
     * Output:
     *  2) output: Polygon in Boost format
     */
    polygon output;
    for (size_t i = 0; i < input.size(); i++) {
        output.outer().push_back(input[i]);
    }
    if (!bg::is_valid(output)) {
        bg::correct(output);
    }

    return output;
}


line BoostPointToBoostLine(std::vector<point> input) {
    /**
     * Function that takes as input a std::vector<point> series of points and converts them to Boost line format
     * 
     * Input:
     *  1) input: Vector of points in Boost format
     * 
     * Output:
     *  2) output: Line in Boost format
     */
    line output;
    for (size_t i = 0; i < input.size(); i++) {
        output.push_back(input[i]);
    }

    return output;
}


polygon cvxpolyxhplane(polygon xy, point m, point n) {
    /** 
     * Function that computes the the intersection of a polygon, with vertex coordinates xy, and a halfplane, defined by a boundary point m and the inward normal n.
     * 
     * Input:
     *  1) xy: Polygon
     *  2) m: Boundary point of the halfplane
     *  2) n: Inward normal of the halfplane
     * 
     * Output:
     *  1) xyNew: Intersection of the polygon with the halfplane
     */

    // Create a dummy origin
    point origin(0.0, 0.0);
    
    // Check if the input polygon is empty
    polygon polyout;
    if (bg::is_empty(xy)) {
        return polyout;
    }

    // Get list of points for the polygon and erase the last element
    std::vector<std::vector<double>> VertexList = BoostPointToStd(BoostPolyToBoostPoint(xy));
    VertexList.pop_back();

    // Compute distance of polygon vertices to the halfspace boundary
    n = point(n.get<0>()/bg::distance(n,origin), n.get<1>()/bg::distance(n,origin)); // Normalize once again
    std::vector<double> dist2hplane(VertexList.size(), 0.0);
    for (size_t i = 0; i < VertexList.size(); i++) {
        dist2hplane[i] = (VertexList[i][0]-m.get<0>())*n.get<0>() + (VertexList[i][1]-m.get<1>())*n.get<1>();
    }

    std::vector<point> VertexListNew;
    size_t numVertex = VertexList.size();
    for (size_t ck = 0; ck < numVertex; ck++) {
        size_t cn = (ck+1)%numVertex;

        if ((dist2hplane[ck]*dist2hplane[cn]) < 0.0) {
            // Compute the point on the boundary and include it into the new vertex list
            double w = ((m.get<0>()-VertexList[cn][0])*n.get<0>() + (m.get<1>()-VertexList[cn][1])*n.get<1>())/((VertexList[ck][0]-VertexList[cn][0])*n.get<0>() + (VertexList[ck][1]-VertexList[cn][1])*n.get<1>());
            point b(w*VertexList[ck][0]+(1-w)*VertexList[cn][0], w*VertexList[ck][1]+(1-w)*VertexList[cn][1]);
            VertexListNew.push_back(b);
        }

        if (dist2hplane[cn] >= 0.0) {
            // Include the next vertex since it is included in the halfspace
            VertexListNew.push_back(point(VertexList[cn][0], VertexList[cn][1]));
        }
    }

    // Finally, if an intersection was found, push back the first element to close the polygon
    if (!VertexListNew.empty()) {
        VertexListNew.push_back(VertexListNew[0]);
    }

    polygon xyNew = BoostPointToBoostPoly(VertexListNew);

    return xyNew;
}


line polyxline(polygon xy, point m, point n) {
    /** 
     * Function that computes the the intersection of a polygon, with vertex coordinates xy, and a line, defined by a point m on the line and the normal vector n.
     * 
     * Input:
     *  1) xy: Polygon
     *  2) m: Point on the line
     *  2) n: Line normal vector
     * 
     * Output:
     *  1) output: Line intersection
     */

    // Create a dummy origin
    point origin(0.0, 0.0);
    
    // Check if v is trivial
    line output;
    if (bg::distance(n,origin) == 0.0) {
        return output;
    }

    // Get list of points for the polygon and erase the last element
    std::vector<std::vector<double>> VertexList = BoostPointToStd(BoostPolyToBoostPoint(xy));
    VertexList.pop_back();

    // Normalize the input vector
    n = {n.get<0>()/bg::distance(n,origin), n.get<1>()/bg::distance(n,origin)};

    // Find distance of all vertices to line
    std::vector<double> dist2line(VertexList.size(), 0.0);
    for (size_t i = 0; i < VertexList.size(); i++) {
        dist2line[i] = (VertexList[i][0]-m.get<0>())*n.get<0>() + (VertexList[i][1]-m.get<1>())*n.get<1>();
    }

    std::vector<point> VertexListNew;
    size_t numVertex = VertexList.size();
    for (size_t ck = 0; ck < numVertex; ck++) {
        if (dist2line[ck] == 0.0) {
            VertexListNew.push_back(point(VertexList[ck][0], VertexList[ck][1]));
        } else {
            size_t cn = (ck+1)%numVertex;
            if (dist2line[ck]*dist2line[cn] < 0.0) {
                double a = -dist2line[cn]/(dist2line[ck]-dist2line[cn]);
                point b(a*VertexList[ck][0]+(1-a)*VertexList[cn][0], a*VertexList[ck][1]+(1-a)*VertexList[cn][1]);
                VertexListNew.push_back(b);
            }
        }
    }

    // Populate line
    output = BoostPointToBoostLine(VertexListNew);
    return output;
}


point polyxray(polygon xy, point b, point v) {
    /** 
     * Function that computes the the intersection of a polygon, with vertex coordinates xy, and a ray, defined by a base point b and the direction vector v.
     * 
     * Input:
     *  1) xy: Polygon
     *  2) b: Base vector
     *  2) v: Direction vector
     * 
     * Output:
     *  1) output: Point intersection
     */

    // Create a dummy origin
    point origin(0.0, 0.0);
    
    // Check if v is trivial
    point output;
    if (bg::distance(v,origin) == 0.0) {
        return output;
    }

    // Get list of points for the polygon and erase the last element
    std::vector<std::vector<double>> VertexList = BoostPointToStd(BoostPolyToBoostPoint(xy));
    VertexList.pop_back();

    // Normalize the input vector
    v = {v.get<0>()/bg::distance(v,origin), v.get<1>()/bg::distance(v,origin)};

    // Get normal vector
    point vn(-v.get<1>(),v.get<0>());

    // Get the intersection with the line itself
    line c = polyxline(xy, b, vn);

    std::vector<point> cvector = BoostLineToBoostPoint(c);
    if (cvector.size() < 2) {
	    return b;
    }
    std::vector<double> a = {(cvector[0].get<0>()-b.get<0>())*v.get<0>() + (cvector[0].get<1>()-b.get<1>())*v.get<1>(), 
                             (cvector[1].get<0>()-b.get<0>())*v.get<0>() + (cvector[1].get<1>()-b.get<1>())*v.get<1>()};
    
    if ((a[0] > 0.0) && (a[1] <= 0.0)) {
        output = point(cvector[0].get<0>(), cvector[0].get<1>());
    } else if ((a[1] > 0.0) && (a[0] <= 0.0)) {
        output = point(cvector[1].get<0>(), cvector[1].get<1>());
    } else {
        output = point(cvector[0].get<0>(), cvector[0].get<1>());
    }

    return output;
}


ProjectionResultStruct polydist(polygon xy, point p) {
    /**
     * Function that computes the distance between a point p and a polygon xy and returns the closest point on the polygon boundary
     * 
     * Input:
     *  1) xy: Polygon
     *  2) p: Point
     * 
     * Output:
     *  1) output: Output that contains the point of minimum distance and the corresponding distance
     */

    // Create a dummy origin
    point origin(0.0, 0.0);

    // Distance to empty set is infinity
    ProjectionResultStruct output;
    if (bg::is_empty(xy)) {
        output.dist = 100000000.0;
        output.projected_point = point(0.0, 0.0);
    }

    // Convert point to std
    std::vector<double> PointCoord = {p.get<0>(), p.get<1>()};

    // Get list of points for the polygon and erase the last element
    std::vector<std::vector<double>> VertexList = BoostPointToStd(BoostPolyToBoostPoint(xy));
    VertexList.pop_back();

    // Construct a new list with all the vertices rolled forward by 1
    size_t numVertex = VertexList.size();
    std::vector<std::vector<double>> VertexListRolled(numVertex, {0.0, 0.0});
    std::vector<std::vector<double>> dxy(numVertex, {0.0, 0.0});
    std::vector<double> diff_norm(numVertex, 0.0);
    for (size_t i = 0; i < numVertex; i++) {
        size_t j = (i+1)%numVertex;
        VertexListRolled[j] = {VertexList[j][0], VertexList[j][1]};
        dxy[j] = {VertexList[j][0]-VertexList[i][0], VertexList[j][1]-VertexList[i][1]};
        diff_norm[j] = bg::distance(point(dxy[j][0],dxy[j][1]), origin);
        if (diff_norm[j] == 0.0) {
            diff_norm[j] = 1.0;
        }
    }

    // Iterate through the edges to find the desired point and distance
    std::vector<double> w(numVertex, 0.0);
    std::vector<double> dtemp(numVertex, 0.0);
    std::vector<std::vector<double>> ctemp(numVertex, {0.0, 0.0});
    for (size_t i = 0; i < numVertex; i++) {
        double w_temp = (PointCoord[0]-VertexList[i][0])*(dxy[i][0]/pow(diff_norm[i],2)) + (PointCoord[1]-VertexList[i][1])*(dxy[i][1]/pow(diff_norm[i],2));
        w[i] = std::max(std::min(w_temp, 1.0), 0.0);
        ctemp[i] = {(1-w[i])*VertexList[i][0] + w[i]*VertexListRolled[i][0], (1-w[i])*VertexList[i][1] + w[i]*VertexListRolled[i][1]};
        dtemp[i] = bg::distance(p, point(ctemp[i][0], ctemp[i][1]));
    }

    // Find the minimum distance and extract the corresponding point
    size_t dist_argmin = std::distance(dtemp.begin(), std::min_element(dtemp.begin(), dtemp.end()));
    double dist = dtemp[dist_argmin];
    std::vector<double> projected_point = ctemp[dist_argmin];

    // Populate and return the output
    output.dist = dist;
    output.projected_point = point(projected_point[0], projected_point[1]);
    return output;
}


ProjectionResultStruct linedist(line xy, point p) {
    /**
     * Function that computes the distance between a point p and a line xy and returns the closest point on the line
     * 
     * Input:
     *  1) xy: Line
     *  2) p: Point
     * 
     * Output:
     *  1) output: Output that contains the point of minimum distance and the corresponding distance
     */

    // Create a dummy origin
    point origin(0.0, 0.0);

    // Distance to empty set is infinity
    ProjectionResultStruct output;
    if (bg::is_empty(xy)) {
        output.dist = 100000000.0;
        output.projected_point = point(0.0, 0.0);
    }

    // Convert point to std
    std::vector<double> PointCoord = {p.get<0>(), p.get<1>()};

    // Get list of points for the polygon and erase the last element
    std::vector<std::vector<double>> VertexList = BoostPointToStd(BoostLineToBoostPoint(xy));

    // Construct a new list with all the vertices rolled forward by 1
    size_t numEdge = VertexList.size()-1;
    std::vector<std::vector<double>> dxy(numEdge, {0.0, 0.0});
    std::vector<double> diff_norm(numEdge, 0.0);
    for (size_t i = 0; i < numEdge; i++) {
        dxy[i] = {VertexList[i+1][0]-VertexList[i][0], VertexList[i+1][1]-VertexList[i][1]};
        diff_norm[i] = bg::distance(point(dxy[i][0],dxy[i][1]), origin);
        if (diff_norm[i] == 0.0) {
            diff_norm[i] = 1.0;
        }
    }

    // Iterate through the edges to find the desired point and distance
    std::vector<double> w(numEdge, 0.0);
    std::vector<double> dtemp(numEdge, 0.0);
    std::vector<std::vector<double>> ctemp(numEdge, {0.0, 0.0});
    for (size_t i = 0; i < numEdge; i++) {
        double w_temp = (PointCoord[0]-VertexList[i][0])*(dxy[i][0]/pow(diff_norm[i],2)) + (PointCoord[1]-VertexList[i][1])*(dxy[i][1]/pow(diff_norm[i],2));
        w[i] = std::max(std::min(w_temp, 1.0), 0.0);
        ctemp[i] = {(1-w[i])*VertexList[i][0] + w[i]*VertexList[i+1][0], (1-w[i])*VertexList[i][1] + w[i]*VertexList[i+1][1]};
        dtemp[i] = bg::distance(p, point(ctemp[i][0], ctemp[i][1]));
    }

    // Find the minimum distance and extract the corresponding point
    size_t dist_argmin = std::distance(dtemp.begin(), std::min_element(dtemp.begin(), dtemp.end()));
    double dist = dtemp[dist_argmin];
    std::vector<double> projected_point = ctemp[dist_argmin];

    // Populate and return the output
    output.dist = dist;
    output.projected_point = point(projected_point[0], projected_point[1]);
    return output;
}


void polytriangulation(std::vector<std::vector<double>> xy, std::vector<std::vector<double>> workspace, bool touching_boundary, std::vector<TriangleClass> *tree) {
    /**
     * Compute the triangulation of the input polygon and its dual (adjacency) graph.
     * 
     * Input:
     *  1) xy: Vertex Coordinates of input polygon - start and end vertices must be the same
     *  2) workspace: Convex boundary of the workspace - start and end vertices must be the same
     *  3) touching_boundary: Flag that is True if the polygon is touching the boundary of the workspace and False otherwise
     *  4) tree: Array of dictionaries with triangles and generated adjacency graph
     * 
     */

    // Eliminate first element
    xy.pop_back();

    // Convert to array notation
    std::vector<std::array<double,2>> xy_array;
    for (size_t i = 0; i < xy.size(); i++) {
        xy_array.push_back({xy[i][0], xy[i][1]});
    }

    // Find triangulation
    std::vector<std::vector<std::array<double,2>>> polygoninput;
    polygoninput.push_back(xy_array);
    std::vector<uint16_t> indices = mapbox::earcut<uint16_t>(polygoninput);

    //for (std::vector<uint16_t>::const_iterator i = indices.begin(); i != indices.end(); ++i)
    //std::cout << *i << ' ';

    // Make the triangles CCW
    //for (std::vector<uint16_t>::iterator it = indices.begin(); it != indices.end(); it = it+3) {
    //    std::iter_swap(it+1, it+2);
    //}

    // Populate an initial list of triangles
    std::vector<std::vector<std::vector<double>>> triangle_list;
    for (size_t i = 0; i < indices.size(); i = i+3) {
        triangle_list.push_back({{xy[indices[i]][0], xy[indices[i]][1]}, {xy[indices[i+1]][0], xy[indices[i+1]][1]}, {xy[indices[i+2]][0], xy[indices[i+2]][1]}});
    }

    // Sort triangles - Area if not touching boundary, Min distance to boundary if touching boundary
    if (!touching_boundary) {
        // Create comparator for area
        struct area_comparator {
            bool operator() (const std::vector<std::vector<double>> a, const std::vector<std::vector<double>> b) {
                return (bg::area(BoostPointToBoostPoly(StdToBoostPoint(a))) > bg::area(BoostPointToBoostPoly(StdToBoostPoint(b))));
            }
        };

        // Create comparator for distance
        class distance_comparator {
            point robotposition_;
            public:
                distance_comparator(point robotposition) : robotposition_(robotposition) {}
                bool operator() (const std::vector<std::vector<double>> a, const std::vector<std::vector<double>> b) {
                    return (bg::distance(robotposition_, BoostPointToBoostPoly(StdToBoostPoint(a))) > bg::distance(robotposition_, BoostPointToBoostPoly(StdToBoostPoint(b))));
                }
        };

        // Sort in order of descending area
        std::sort(triangle_list.begin(), triangle_list.end(), area_comparator());
    } else {
        // Construct line object for boundary
        line workspaceLine;
        for (size_t i = 0; i < workspace.size(); i++) {
            workspaceLine.push_back(point(workspace[i][0], workspace[i][1]));
        }
        // Create comparator for distance to the boundary
        class boundary_distance_comparator {
            line workspaceLine_;
            public:
                boundary_distance_comparator(line workspaceLine) : workspaceLine_(workspaceLine) {}
                bool operator() (const std::vector<std::vector<double>> a, const std::vector<std::vector<double>> b) {
                    return (bg::distance(BoostPointToBoostPoly(StdToBoostPoint(a)), workspaceLine_) < bg::distance(BoostPointToBoostPoly(StdToBoostPoint(b)), workspaceLine_));
                }
        };

        // Sort in order of ascending distance to the boundary
        std::sort(triangle_list.begin(), triangle_list.end(), boundary_distance_comparator(workspaceLine));
    }

    // Construct the first node of the tree that will act as the root
    TriangleClass root_triangle;
    root_triangle.set_vertices(StdToBoostPoint(triangle_list[0]));
    root_triangle.set_predecessor(-1);
    root_triangle.set_depth(0);
    root_triangle.set_index(0);
    uint16_t tree_index = 0;
    tree->push_back(root_triangle);

    // Initialize search
    triangle_list.erase(triangle_list.begin());
    std::vector<TriangleClass> stack;
    stack.push_back((*tree)[0]);

    // Build the tree by expanding nodes until the stack is empty
    // (The stack will be empty when the leaf nodes consist of only one edge)
    while (stack.size() != 0) {
        // Pop the first element of the stack and delete it from the stack
        TriangleClass expanded_node = stack[0];
        stack.erase(stack.begin());
        size_t i = 0;
        while (i < triangle_list.size()) {
            // Construct two edge arrays: one for the parent and one for the candidate child
            // Orient the parent CCW as desired and the child CW to check for collisions
            std::vector<point> expanded_node_vertices = expanded_node.get_vertices();
            std::vector<std::vector<point>> polygon1_edges, polygon2_edges;
            polygon1_edges.push_back({point(expanded_node_vertices[0].get<0>(), expanded_node_vertices[0].get<1>()), 
                                      point(expanded_node_vertices[1].get<0>(), expanded_node_vertices[1].get<1>())});
            polygon1_edges.push_back({point(expanded_node_vertices[1].get<0>(), expanded_node_vertices[1].get<1>()), 
                                      point(expanded_node_vertices[2].get<0>(), expanded_node_vertices[2].get<1>())});
            polygon1_edges.push_back({point(expanded_node_vertices[2].get<0>(), expanded_node_vertices[2].get<1>()), 
                                      point(expanded_node_vertices[0].get<0>(), expanded_node_vertices[0].get<1>())});
            polygon2_edges.push_back({point(triangle_list[i][0][0], triangle_list[i][0][1]), 
                                      point(triangle_list[i][2][0], triangle_list[i][2][1])});
            polygon2_edges.push_back({point(triangle_list[i][2][0], triangle_list[i][2][1]), 
                                      point(triangle_list[i][1][0], triangle_list[i][1][1])});
            polygon2_edges.push_back({point(triangle_list[i][1][0], triangle_list[i][1][1]), 
                                      point(triangle_list[i][0][0], triangle_list[i][0][1])});
            bool triangles_touch = false;
            size_t adj_edge_index = -1;
            for (size_t polygon1_edge_index = 0; polygon1_edge_index < polygon1_edges.size(); polygon1_edge_index++) {
                for (size_t polygon2_edge_index = 0; polygon2_edge_index < polygon2_edges.size(); polygon2_edge_index++) {
                    line line1 = BoostPointToBoostLine(polygon1_edges[polygon1_edge_index]);
                    line line2 = BoostPointToBoostLine(polygon2_edges[polygon2_edge_index]);
                    if (bg::equals(line1,line2)) {
                        triangles_touch = true;
                        adj_edge_index = polygon2_edge_index;

                        // Do something to complete the iterations
                        polygon1_edge_index = polygon1_edges.size();
                        break;
                    }
                }
            }

            // Check if the triangles touch, otherwise continue
            if (!triangles_touch) {
                i++;
                continue;
            } else {
                // Add the child to the tree
                tree_index++;
                TriangleClass new_triangle;
                new_triangle.set_predecessor(expanded_node.get_index());
                new_triangle.set_depth(expanded_node.get_depth()+1);
                new_triangle.set_index(tree_index);
                new_triangle.set_adj_edge(polygon2_edges[adj_edge_index]);

                // Find the 3rd point of the child triangle (that does not belong to the shared edge) and arrange the vertices so that this is the 3rd vertex
                size_t third_vertex_index;
                switch(adj_edge_index) {
                    case 0: 
                        third_vertex_index = 1;
                        break;
                    case 1:
                        third_vertex_index = 0;
                        break;
                    case 2:
                        third_vertex_index = 2;
                        break;
                }
                new_triangle.set_vertices({polygon2_edges[adj_edge_index][1], polygon2_edges[adj_edge_index][0], point(triangle_list[i][third_vertex_index][0], triangle_list[i][third_vertex_index][1])});

                // Delete the child from the input stack and append it to the stack to be expanded
                triangle_list.erase(triangle_list.begin()+i);
                stack.push_back(new_triangle);

                // Add the triangle to the list
                tree->push_back(new_triangle);
            }
        }
    }

    // Create comparator for depth
    struct depth_comparator {
        bool operator() (const TriangleClass a, const TriangleClass b) {
            return (a.get_depth() > b.get_depth());
        }
    };

    // As a final step, sort the tree as a stack in order of descending depth
    std::sort(tree->begin(), tree->end(), depth_comparator());

    // Make sure to change the node and predecessor indices to indicate the index change
    std::vector<size_t> indices_new(tree->size(), 0);
    std::vector<size_t> indices_old(tree->size(), 0);
    for (size_t i = 0; i < tree->size(); i++) {
        indices_new[i] = i;
        indices_old[i] = (*tree)[i].get_index();
    }

    // Update indices and predecessors appropriately
    for (size_t i = 0; i < tree->size()-1; i++) {
        std::vector<size_t>::iterator it = std::find(indices_old.begin(), indices_old.end(), (*tree)[i].get_predecessor());
        (*tree)[i].set_predecessor(indices_new[std::distance(indices_old.begin(),it)]);
        (*tree)[i].set_index(i);
    }
    (*tree)[tree->size()-1].set_index(tree->size()-1);

    return;
}


void polyconvexdecomposition(std::vector<std::vector<double>> xy, std::vector<std::vector<double>> workspace, bool touching_boundary, std::vector<PolygonClass> *tree) {
    /**
     * Compute the convex decomposition of the input polygon and its dual (adjacency) graph.
     * 
     * Input:
     *  1) xy: Vertex Coordinates of input polygon - start and end vertices must be the same
     *  2) workspace: Convex boundary of the workspace - start and end vertices must be the same
     *  3) touching_boundary: Flag that is True if the polygon is touching the boundary of the workspace and False otherwise
     *  4) tree: Array of dictionaries with polygons and generated adjacency graph
     * 
     */

    // Eliminate first element
    xy.pop_back();

    // Convert to CGAL array notation
    CGAL_Polygon_2 cgal_polygon;
    CGAL_Polygon_list partition_polys;
    CGAL_Traits partition_traits;
    CGAL::set_pretty_mode(std::cout);
    for (size_t i = 0; i < xy.size(); i++) {
        cgal_polygon.push_back(CGAL_Point_2(xy[i][0], xy[i][1]));
    }

    // Find convex decomposition
    std::vector<std::vector<std::vector<double>>> polygon_list;
    CGAL::optimal_convex_partition_2(cgal_polygon.vertices_begin(),
                                     cgal_polygon.vertices_end(), 
                                     std::back_inserter(partition_polys),
                                     partition_traits);
    assert(CGAL::convex_partition_is_valid_2(cgal_polygon.vertices_begin(),
                                             cgal_polygon.vertices_end(),
                                             partition_polys.begin(),
                                             partition_polys.end(),
                                             partition_traits));
    
    for (size_t k = 0; k < partition_polys.size(); k++) {
        auto it = std::next(partition_polys.begin(), k);
        CGAL_Polygon_2 next_polygon = *it;
        std::vector<point> next_polygon_vertices = {};
        for (size_t l = 0; l < next_polygon.size(); l++) {
            next_polygon_vertices.push_back(point(next_polygon.vertex(l).x(), next_polygon.vertex(l).y()));
        }
        polygon_list.push_back(BoostPointToStd(next_polygon_vertices));
    }

    // Sort polygons - Area if not touching boundary, Min distance to boundary if touching boundary
    if (!touching_boundary) {
        // Create comparator for area
        struct area_comparator {
            bool operator() (const std::vector<std::vector<double>> a, const std::vector<std::vector<double>> b) {
                return (bg::area(BoostPointToBoostPoly(StdToBoostPoint(a))) > bg::area(BoostPointToBoostPoly(StdToBoostPoint(b))));
            }
        };

        // Sort in order of descending area
        std::sort(polygon_list.begin(), polygon_list.end(), area_comparator());
    } else {
        // Construct line object for boundary
        line workspaceLine;
        for (size_t i = 0; i < workspace.size(); i++) {
            workspaceLine.push_back(point(workspace[i][0], workspace[i][1]));
        }
        // Create comparator for distance to the boundary
        class boundary_distance_comparator {
            line workspaceLine_;
            public:
                boundary_distance_comparator(line workspaceLine) : workspaceLine_(workspaceLine) {}
                bool operator() (const std::vector<std::vector<double>> a, const std::vector<std::vector<double>> b) {
                    double size_a = static_cast<double>(a.size());
                    double size_b = static_cast<double>(b.size());
                    double sum_a = 0.0;
                    double sum_b = 0.0;
                    double boundary_a = 0.0;
                    double boundary_b = 0.0;
                    for (size_t i = 0; i < a.size(); i++) {
                        sum_a = sum_a + bg::distance(point(a[i][0],a[i][1]), workspaceLine_);
                        if (bg::distance(point(a[i][0],a[i][1]), workspaceLine_) < 1e-2) {
                            boundary_a = boundary_a + 1.0;
                        }
                    }
                    for (size_t i = 0; i < b.size(); i++) {
                        sum_b = sum_b + bg::distance(point(b[i][0],b[i][1]), workspaceLine_);
                        if (bg::distance(point(b[i][0],b[i][1]), workspaceLine_) < 1e-2) {
                            boundary_b = boundary_b + 1.0;
                        }
                    }
                    double normalized_a = sum_a/size_a;
                    double normalized_b = sum_b/size_b;
                    return (boundary_a/normalized_a > boundary_b/normalized_b);
                }
        };

        // Sort in order of ascending distance to the boundary
        std::sort(polygon_list.begin(), polygon_list.end(), boundary_distance_comparator(workspaceLine));
    }

    // Construct the first node of the tree that will act as the root
    PolygonClass root_polygon;
    root_polygon.set_vertices(StdToBoostPoint(polygon_list[0]));
    root_polygon.set_predecessor(-1);
    root_polygon.set_depth(0);
    root_polygon.set_index(0);
    uint16_t tree_index = 0;
    tree->push_back(root_polygon);

    // Initialize search
    polygon_list.erase(polygon_list.begin());
    std::vector<PolygonClass> stack;
    stack.push_back((*tree)[0]);

    // Build the tree by expanding nodes until the stack is empty
    // (The stack will be empty when the leaf nodes consist of only one edge)
    while (stack.size() != 0) {
        // Pop the first element of the stack and delete it from the stack
        PolygonClass expanded_node = stack[0];
        stack.erase(stack.begin());

        // Find edges of expanded node - CW
        std::vector<point> expanded_node_vertices = expanded_node.get_vertices();
        std::vector<std::vector<point>> polygon1_edges;
        for (size_t j = 0; j < expanded_node_vertices.size(); j++) {
            polygon1_edges.push_back({point(expanded_node_vertices[(j+1)%expanded_node_vertices.size()].get<0>(), expanded_node_vertices[(j+1)%expanded_node_vertices.size()].get<1>()), 
                                      point(expanded_node_vertices[j%expanded_node_vertices.size()].get<0>(), expanded_node_vertices[j%expanded_node_vertices.size()].get<1>())});
        }

        size_t i = 0;
        while (i < polygon_list.size()) {
            // Find edges of candidate child - CCW
            std::vector<std::vector<point>> polygon2_edges;
            for (size_t j = 0; j < polygon_list[i].size(); j++) {
                polygon2_edges.push_back({point(polygon_list[i][j%polygon_list[i].size()][0], polygon_list[i][j%polygon_list[i].size()][1]),
                                          point(polygon_list[i][(j+1)%polygon_list[i].size()][0], polygon_list[i][(j+1)%polygon_list[i].size()][1])});
            }
            bool polygons_touch = false;
            size_t adj_edge_index = -1;
            for (size_t polygon1_edge_index = 0; polygon1_edge_index < polygon1_edges.size(); polygon1_edge_index++) {
                for (size_t polygon2_edge_index = 0; polygon2_edge_index < polygon2_edges.size(); polygon2_edge_index++) {
                    line line1 = BoostPointToBoostLine(polygon1_edges[polygon1_edge_index]);
                    line line2 = BoostPointToBoostLine(polygon2_edges[polygon2_edge_index]);
                    if (bg::equals(line1,line2)) {
                        polygons_touch = true;
                        adj_edge_index = polygon2_edge_index;

                        // Do something to complete the iterations
                        polygon1_edge_index = polygon1_edges.size();
                        break;
                    }
                }
            }

            // Check if the polygons touch, otherwise continue
            if (!polygons_touch) {
                i++;
                continue;
            } else {
                // Add the child to the tree
                PolygonClass new_polygon;
                new_polygon.set_predecessor(expanded_node.get_index());
                new_polygon.set_depth(expanded_node.get_depth()+1);
                new_polygon.set_index(++tree_index);
                new_polygon.set_adj_edge(polygon2_edges[adj_edge_index]);

                // Find the vertices of the added polygon
                std::vector<point> new_polygon_vertices;
                for (size_t j = 0; j < polygon_list[i].size(); j++) {
                    new_polygon_vertices.push_back(point(polygon_list[i][(adj_edge_index+j)%polygon_list[i].size()][0], polygon_list[i][(adj_edge_index+j)%polygon_list[i].size()][1]));
                }
                new_polygon.set_vertices(new_polygon_vertices);
                
                // As a final preprocessing step, check whether the edges before and after the adj_edge are parallel with adj_edge
                // If that's the case, cut the triangle corresponding to that edge as an extra polygon
                double dist_last_first = bg::distance(new_polygon_vertices[0], new_polygon_vertices.back());
                double dist_second_third = bg::distance(new_polygon_vertices[2], new_polygon_vertices[1]);
                double dist_adj_edge = bg::distance(new_polygon_vertices[1], new_polygon_vertices[0]);
                std::vector<double> tangent_before = {(new_polygon_vertices[0].get<0>()-new_polygon_vertices.back().get<0>())/dist_last_first, 
                                                      (new_polygon_vertices[0].get<1>()-new_polygon_vertices.back().get<1>())/dist_last_first};
                std::vector<double> tangent_after = {(new_polygon_vertices[2].get<0>()-new_polygon_vertices[1].get<0>())/dist_second_third, 
                                                     (new_polygon_vertices[2].get<1>()-new_polygon_vertices[1].get<1>())/dist_second_third};
                std::vector<double> tangent_adj_edge = {(new_polygon_vertices[1].get<0>()-new_polygon_vertices[0].get<0>())/dist_adj_edge, 
                                                        (new_polygon_vertices[1].get<1>()-new_polygon_vertices[0].get<1>())/dist_adj_edge};
                std::vector<double> normal_adj_edge = {-tangent_adj_edge[1], tangent_adj_edge[0]};
                if (abs(tangent_before[0]*normal_adj_edge[0]+tangent_before[1]*normal_adj_edge[1]) < 0.001) {
                    // Add triangle
                    PolygonClass polygon_before;
                    polygon_before.set_predecessor(new_polygon.get_index());
                    polygon_before.set_depth(new_polygon.get_depth()+1);
                    polygon_before.set_index(++tree_index);
                    std::vector<point> polygon_before_adj_edge = {new_polygon_vertices[0], new_polygon_vertices[new_polygon_vertices.size()-2]};
                    std::vector<point> polygon_before_vertices = {new_polygon_vertices[0], new_polygon_vertices[new_polygon_vertices.size()-2], new_polygon_vertices[new_polygon_vertices.size()-1]};
                    polygon_before.set_adj_edge(polygon_before_adj_edge);
                    polygon_before.set_vertices(polygon_before_vertices);

                    // Delete the last vertex from the original polygon
                    new_polygon_vertices.pop_back();
                    new_polygon.set_vertices(new_polygon_vertices);

                    // Add the new triangle to the tree
                    tree->push_back(polygon_before);
                    stack.push_back(polygon_before);
                }

                if (abs(tangent_after[0]*normal_adj_edge[0]+tangent_after[1]*normal_adj_edge[1]) < 0.001) {
                    // Add triangle
                    PolygonClass polygon_after;
                    polygon_after.set_predecessor(new_polygon.get_index());
                    polygon_after.set_depth(new_polygon.get_depth()+1);
                    polygon_after.set_index(++tree_index);
                    std::vector<point> polygon_after_adj_edge = {new_polygon_vertices[3], new_polygon_vertices[1]};
                    std::vector<point> polygon_after_vertices = {new_polygon_vertices[3], new_polygon_vertices[1], new_polygon_vertices[2]};
                    polygon_after.set_adj_edge(polygon_after_adj_edge);
                    polygon_after.set_vertices(polygon_after_vertices);

                    // Delete the second vertex from the original polygon
                    new_polygon_vertices.erase(new_polygon_vertices.begin()+2);
                    new_polygon.set_vertices(new_polygon_vertices);

                    // Add the new triangle to the tree and stack
                    tree->push_back(polygon_after);
                    stack.push_back(polygon_after);
                }

                // Delete the child from the input stack
                polygon_list.erase(polygon_list.begin()+i);

                // Add the polygon to the tree and stack
                tree->push_back(new_polygon);
                stack.push_back(new_polygon);
            }
        }
    }

    // Create comparator for depth
    struct depth_comparator {
        bool operator() (const PolygonClass a, const PolygonClass b) {
            return (a.get_depth() > b.get_depth());
        }
    };

    // As a final step, sort the tree as a stack in order of descending depth
    std::sort(tree->begin(), tree->end(), depth_comparator());

    // Make sure to change the node and predecessor indices to indicate the index change
    std::vector<size_t> indices_new(tree->size(), 0);
    std::vector<size_t> indices_old(tree->size(), 0);
    for (size_t i = 0; i < tree->size(); i++) {
        indices_new[i] = i;
        indices_old[i] = (*tree)[i].get_index();
    }

    // Update indices and predecessors appropriately
    for (size_t i = 0; i < tree->size()-1; i++) {
        std::vector<size_t>::iterator it = std::find(indices_old.begin(), indices_old.end(), (*tree)[i].get_predecessor());
        (*tree)[i].set_predecessor(indices_new[std::distance(indices_old.begin(),it)]);
        (*tree)[i].set_index(i);
    }
    (*tree)[tree->size()-1].set_index(tree->size()-1);

    return;
}
