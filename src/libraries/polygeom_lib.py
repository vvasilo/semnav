#!/usr/bin/env python

"""
MIT License (modified)

Copyright (c) 2020 The Trustees of the University of Pennsylvania
Authors:
Omur Arslan <omur@seas.upenn.edu>
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

from __future__ import division
import numpy as np
import matplotlib.path as mpath
import shapely as sp
import tripy, time, math
from shapely.geometry import Polygon
from shapely.ops import cascaded_union
from operator import itemgetter
from itertools import groupby


def polysignarea(xy):
    """
    polysignarea(xy) determines the signed area of a non-self-intersecting 
    polygon with vertices xy
    
    Input:
        xy   : Vertex coordinated of a non-self-intersecting polygon
               (Nx2 numpy.array)   
    Output:
        area : Signed area of the polygon
    Usage:
        import numpy as np
        from cvxpolygeom import polysignarea 
        xy = np.array([[0,0],[0,1],[1,0]])
        area = polysignarea(xy)
    """
    xy = xy.reshape(-1,2) # Convert the input data into a 2D array 
    numVertex = xy.shape[0] # Number of vertices
    area = 0.0
    for ck in range(0,numVertex):
        cn = (ck + 1) % numVertex
        area = area + np.cross(xy[ck],xy[cn])
    area = 0.5*area

    return area    


def polyarea(xy):
    """
    polyarea(xy) determines the area of a non-self-intersecting polygon 
    with vertices xy
    Input:
        xy   : Vertex coordinated of a non-self-intersecting polygon
               (Nx2 numpy.array)   
    Output:
        area : Area of the polygon
    Usage:
        import numpy as np
        from cvxpolygeom import polyarea 
        xy = np.array([[0,0],[0,1],[1,0]])
        area = polyarea(xy)
    """

    return abs(polysignarea(xy))


def ispolycw(xy):
    """
    ispolycw(xy) determines if the vertices, xy, of a non-self-intersecting polygon 
    are in clockwise order. Its computation is based on the signed are of the polygon. 
    
    Input:
        xy : Vertex coordinated of a non-self-intersecting polygon
             (Nx2 numpy.array)   
    Output:
        cw : a boolean variable which is True if the input polygon is in clockwise order 
             (Boolean [True/False])
    Usage:
        import numpy as np
        from cvxpolygeom import ispolycw 
        xy = np.array([[0,0],[0,1],[1,0]])
        cw = ispolycw(xy)
    """
    return (polysignarea(xy) <= 0)

def inpolygon(xy, p):
    """
    inpolygon(xy, p) determines if a given set of points, p, are contained in 
    a polygon, with vertex set xy.
    Input:  
        xy : Vertex coordinates of a polygon
              (Nx2 numpy.array)
        p  : Coordinates of a set of points
              (Mx2 numpy.array)
    Output: 
        I  : a boolean array indicating which points are contained in the polygon 
    Usage:
        import numpy as np
        from cvxpolygeom import inpolygon
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches 
        n = 10
        v = 8
        p = 2 * np.random.rand(n,2) - 1
        th = np.linspace(0, 2*np.pi, v)
        xy = np.array([np.cos(th), np.sin(th)]).T
        I = inpolygon(xy, p)
        pIn = p[I]
        pOut = p[~I]
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.add_patch(mpatches.Polygon(xy, closed=True, alpha=0.5))
        plt.plot(pIn[:,0], pIn[:,1], 'g*', pOut[:,0], pOut[:,1], 'ro')
        ax.axis('equal') 
        fig.show()   
    """
    
    # Convert input data into 2D arrays
    xy = xy.reshape(-1,2)
    p = p.reshape(-1,2)
    
    # Create a path decribing the polygon and check if each point is contained in the polygon
    polypath = mpath.Path(xy)
    I = polypath.contains_points(p)
    
    return I    


def polydist(xy, p):
    """
    polydist(xy, p) computes the distance between a set of points, p, and 
    a polygon, xy, and return the closest points on the polygon boundary.   
    Here, distance is defined as the minimum distance between an input 
    point and any point on the polygon boundary.

    Input:  
        xy : Vertex coordinates of a polygon
              (Nx2 numpy.array)
        p  : Coordinates of a set of points
              (Mx2 numpy.array)
    Output: 
        D  : Distance between points and the polygon 
        C  : Coordinates of the closest points on the polygon to the input points
    Usage:
        import numpy as np
        from cvxpolygeom import polydist 
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        n = 2
        v = 7
        p = 2 * np.random.rand(n,2) - 1
        th = np.linspace(0, 2*np.pi, v)
        xy = np.array([np.cos(th), np.sin(th)]).T
        D, C = polydist(xy, p)
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.add_patch(mpatches.Polygon(xy, closed=True, alpha=0.5))
        plt.plot(p[:,0],p[:,1],'ro')
        plt.plot(C[:,0],C[:,1], 'r*')   
        LX = np.array([p[:,0],C[:,0]])
        LY = np.array([p[:,1],C[:,1]])  
        plt.plot(LX, LY, 'r-')
        ax.axis('equal') 
        fig.show()   
    """
    # Convert input data into 2D arrays
    xy = xy.reshape(-1,2)
    p = p.reshape(-1,2)
     
    # Distance to empty set is infinity
    if (xy.shape[0] == 0):
        D = np.zeros(p.shape[0])
        D.fill(np.inf)
        C = np.zeros(p.shape)
        C.fill(np.inf) 
        return D,C
    
    orientsign = 1 - 2 * ispolycw(xy) # orientation of the polygon
    numPoint = p.shape[0] # number of points
    # Relative coordinates of polygon rims
    xyPre = np.roll(xy,1, axis=0)
    dxy = xyPre - xy
    dxyNorm = np.power(np.linalg.norm(dxy,axis=1)[:,np.newaxis],2)
    dxyNorm[(dxyNorm==0)] = 1

    # Compute distances and closest points on the polygon boundary  
    D = np.zeros(numPoint)
    C = np.zeros([numPoint,2])
    for k in range(numPoint):
        w = np.sum((p[k] - xy)*dxy,axis=1)[:,np.newaxis]/dxyNorm
        w = np.fmax(np.fmin(w,1),0)
        ctemp = (1-w)*xy + w*xyPre
        dtemp = np.linalg.norm(p[k] - ctemp, axis=1)
        iMin = dtemp.argmin()
        D[k] = dtemp[iMin]
        C[k] = ctemp[iMin]  
    
    return D,C  


def cvxpolyxhplane(xy, m, n):
    """
    cvxpolyxhplane(xy, m, n) computes the intersection of a polygon, 
    with vertex coordinates xy, and a halfplane, defined by a boundary
    point m and the inward normal n.
    
    Input:
        xy   : Vertex coordinates of a polygon 
               (Nx2 numpy.array)
        m, n : A boundary point and the inward normal of a halfplane 
               (1x2 numpy.array)
    Output:
        xyNew : Vertex coordinates of the intersection
    Usage:
        import numpy as np
        from cvxpolygeom import cvxpolyxhplane 
        xy = np.array([[0,0],[0,1],[1,1],[1,0]])
        m = np.array([0.25,0.25])
        n = np.array([1,1])
        xyNew = cvxpolyxhplane(xy,m,n)
    """

    # Check if the input polygon is empty
    if (xy.size == 0):  
        return xy

    # Compute distance of polygon vertices to the halfspace boundary    
    n = n/np.linalg.norm(n) # Normalize this vector once again
    dist2hplane = np.dot(xy - m, n)
    
    xyNew = []
    numVertex =xy.shape[0]
    for ck in range(0,numVertex):
        cn = (ck + 1) % numVertex
        
        if ((dist2hplane[ck]*dist2hplane[cn]) < 0):
        # Compute the point on the boundary and include it into the new vertex list
            w = np.dot(m - xy[cn],n)/np.dot(xy[ck] - xy[cn],n)
            b = w*xy[ck] + (1 - w)*xy[cn]
            xyNew.append(b.tolist())

        if (dist2hplane[cn] >= 0): 
            #Include the next vertex since it is included in the halfspace 
            xyNew.append(xy[cn].tolist())
    xyNew = np.array(xyNew)
    xyNew = xyNew.reshape(-1,2)
    
    return xyNew

def polyxline(xy, m, n):
    """
    polyxline(xy,m,n) computes the intersection of the boundary of a polygon, 
    with vertex coordinates xy, and a line, defined by a point m on the line 
    and its normal n.      

    Input:
        xy    : Vertex coordinates of a polygon
        m, n  : Line parameters,  m is a point on the line and n is its normal 
    Output:
        xyNew : Vertex coordinates of the intersecting points of the polygon
                boundary and line
    Usage:
      import numpy as np
      import matplotlib.pyplot as plt
      import matplotlib.patches as mpatches
      from cvxpolygeom import polyxline
      xy = np.array([[0,0],[0,1],[1,1],[1,0]])
      m = np.array([0.5, 0.5])
      n = np.array([np.cos(np.pi/6), np.sin(np.pi/6)])
      xyNew = polyxline(xy, m, n)
      fig = plt.figure()
      ax = fig.add_subplot(1,1,1)
      ax.add_patch(mpatches.Polygon(xy, closed=True, facecolor='#FFFF00', edgecolor='b'))
      plt.plot(xyNew[:,0], xyNew[:,1], 'ro')
      ax.axis('equal')
      fig.show()
    """
    if (np.linalg.norm(n) == 0):
        xyNew = []
    else:
        xy = xy.reshape(-1,2) # Convert the input polygon into 2D array
        n = n/np.linalg.norm(n) # Normalize the line normal vector 
        dist2line = np.dot(xy-m,n) # Signed perpendicular distance to the line
        numVertex = xy.shape[0] # Number of vertices
        xyNew = []
        for ck in range(numVertex):
            if (dist2line[ck] == 0):
                xyNew.append(xy[ck])
            else:
                cn = (ck+1) % numVertex # Next vertex index
                if ((dist2line[ck]*dist2line[cn]) < 0):
                    a = -dist2line[cn]/(dist2line[ck] - dist2line[cn])
                    xyNew.append(a*xy[ck]+(1-a)*xy[cn])
    xyNew = np.array(xyNew) 
    return xyNew           


def cvxpolyerode(xy, r):
    """
    Erosion (Contraction) of a convex polygon, with vertices xy, by a closed 
    circle of radius r

    Input: 
        xy    : Vertex Coordinates of a convex polygon
                (Nx2 numpy.array)
        r     : Radius of the erosion (contraction) disk
    Output:
        xyNew : Coordinates of contracted convex polygon
    Usage:    
        import numpy as np
        from cvxpolygeom import cvxpolyerode 
        xy = np.array([[0,0],[0,1],[1,1],[1,0]])
        r = 0.25
        xyNew = cvxpolyerode(xy,r)
    """

    numVertex = xy.shape[0] # Number of vertices
    if (numVertex < 3):
    # The input polygon is trivial
        if (r > 0):
            xyNew = np.array([[]])
        else: 
            xyNew = xy
    else:
        # For nontrivial convex polygons: 
        
        # Determine if the input polygon is in clockwise or counter-clockwise order
        orientsign = 1 - 2 * ispolycw(xy)
    
        # Compute the contraction of input polygon
        xyNew = xy
        for ck in range(0,numVertex):
            cp = (ck - 1) % numVertex # Previous vertex index        
            ek = xy[ck] - xy[cp] # Edge vector        
            if (np.linalg.norm(ek) > 0):
                # Erode any nontrivial edge
                nk = np.array([-ek[1], ek[0]]) 
                nk = orientsign * nk/np.linalg.norm(nk) # inward pointing edge normal
                mk = xy[ck] + r * nk # A point on the edge
                xyNew = cvxpolyxhplane(xyNew, mk, nk) 

    return xyNew  

def cvxpolyintersect(xy1,xy2):
    """
    Compute the intersection of two convex polygons, with vertices xy1 and xy2

    Input:
        xy1     : Vertex Coordinates of a convex polygon
                  (Nx2 numpy.array)
        xy2     : Vertex Coordinates of a convex polygon
                  (Nx2 numpy.array)
    Output:
        xyNew   : Coordinates of polygon intersection (polygon, line, point or empty)
    Usage:    
        import numpy as np
        from shapely.geometry import Polygon
        from cvxpolygeom import cvxpolyintersect 
        xy1 = np.array([[0,0],[0,1],[1,1],[1,0]])
        xy2 = np.array([[0.5,0],[0,1],[1,1],[1,0]])
        xyNew = cvxpolyintersect(xy1,xy2)
    """
    numVertex1 = xy1.shape[0] # Number of vertices
    numVertex2 = xy2.shape[0] # Number of vertices
    polygon1 = Polygon(xy1)
    polygon2 = Polygon(xy2)
    res = polygon1.intersection(polygon2)
    if isinstance(res,sp.geometry.polygon.Polygon) == True:
        vertlist = list(res.exterior.coords)
        numVertex = len(vertlist)-1
        xyNew = np.zeros((numVertex,2))
        for ctr in range(0,numVertex):
            xyNew[ctr,0] = vertlist[ctr][0]
            xyNew[ctr,1] = vertlist[ctr][1]
    elif isinstance(res,sp.geometry.linestring.LineString) == True:
        vertlist = list(res.coords)
        xyNew = np.zeros((2,2))
        xyNew[0,0] = vertlist[0][0]
        xyNew[0,1] = vertlist[0][1]
        xyNew[1,0] = vertlist[1][0]
        xyNew[1,1] = vertlist[1][1]
    elif isinstance(res,sp.geometry.point.Point) == True:
        vertlist = list(res.coords)
        xyNew = np.zeros((1,2))
        xyNew[0,0] = vertlist[0][0]
        xyNew[0,1] = vertlist[0][1]
    else:
        xyNew = np.array([[]])
        
    return xyNew

def polyxray(x,b,v):
    """
    Compute the intersection of the boundary of a polygon,
    with vertex coordinates x and y, and a ray, defined by
    a base point b on the line and a direction vector v.

    Input:
        x       : Vertex Coordinates of a convex polygon
                  (Nx2 numpy.array)
        b, v    : Ray parameters starting at the base point b and in direction v 
                  (1x2 numpy.array)
    Output:
        c       : Vertex coordinates of the intersecting points of the polygon boundary and ray
    Usage:
        import numpy as np
        from shapely.geometry import Polygon
        from cvxpolygeom import polyxray, polyxline
        x = np.array([[0,0],[0,1],[1,1],[1,0]])
        b = np.array([0.5, 0.5])
        v = np.array([np.cos(np.pi/6), np.sin(np.pi/6)])
        c = polyxray(x,b,v)
    """
    if (np.linalg.norm(v) == 0):
        c = []
    else:
        x = x.reshape(-1,2) # Convert the input polygon into 2D array
        v = v/np.linalg.norm(v) # Normalize the input vector
        vn = np.array([-v[1],v[0]])
        c = polyxline(x,b,vn)
        a = np.array([(c[0][0]-b[0])*v[0]+(c[0][1]-b[1])*v[1],(c[1][0]-b[0])*v[0]+(c[1][1]-b[1])*v[1]])
        if a[0]>0 and a[1]<=0: # Choose only the point in the positive direction
            c = np.array([c[0][0],c[0][1]])
        elif a[1]>0 and a[0]<=0:
            c = np.array([c[1][0],c[1][1]])
        else:
            c = np.array([c[0][0],c[0][1]])
    return c

def polydilate(xy,offset):
    """
    Compute the dilation of a polygon by a fixed offset.

    Input:
        xy      : Vertex Coordinates of a polygon
                  (Nx2 numpy.array)
        offset  : Offset of dilation
    Output:
        xyNew   : Coordinates of dilated polygon
                  (Nx2 numpy.array)
    """

    # Construct a polygon based on the input coordinate vertices
    polygon_input = Polygon(xy)

    # Find the output after the offset
    polygon_output = polygon_input.buffer(offset, join_style=2)
    polygon_output = sp.geometry.polygon.orient(polygon_output, 1.0) # orient polygon to be CCW

    # Compute the coordinates
    xyNew = np.array(polygon_output.exterior.coords.xy).transpose()

    return xyNew

def polytriangulation(xy,workspace,touching_boundary):
    """
    Compute the triangulation of the input polygon and its dual (adjacency) graph.

    Input:
        xy                  : Vertex Coordinates of input polygon - start and end vertices must be the same
                              (Nx2 numpy.array)
        workspace           : Convex boundary of the workspace - start and end vertices must be the same
                              (Nx2 numpy.array)
        touching_boundary   : Flag that is True if the polygon is touching the boundary of the workspace and False otherwise
    Output:
        tree    : Array of dictionaries with triangles and generated adjacency graph
                  Each dictionary contains:
                    1) 'vertices': vertices of each triangle (arranged CCW - start and end vertices are NOT the same) - for each of the children, the vertices of the adjacency edge are the first two vertices
                    2) 'predecessor': the index of the triangle predecessor in the adjacency tree (-1 for the root)
                    3) 'depth': the depth of the (triangle) node in the adjacency tree (0 for the root)
                    4) 'index': the index of the triangle in the tree (its serial number)
                    5) 'adj_edge': the edge the triangle shares with its predecessor (CCW oriented with respect to the triangle)
    """
    # Construct a polygon based on the input coordinate vertices 
    polygon_in = Polygon(np.array(xy))

    # Find polygon vertices (except the last one because of tripy)
    polygon_vertices = np.vstack((polygon_in.exterior.coords.xy[0][0:-1], polygon_in.exterior.coords.xy[1][0:-1])).transpose()

    # Find triangulation
    triangles = np.array(tripy.earclip(polygon_vertices))

    # Sorting argument for triangles - Area if not touching boundary, Min distance to boundary if touching boundary
    sorting = np.zeros(triangles.shape[0])
    if not touching_boundary:
        for i in range(triangles.shape[0]):
            sorting[i] = -polyarea(triangles[i])
    elif touching_boundary:
        for i in range(triangles.shape[0]):
            D, C = polydist(workspace[0:-1],triangles[i])
            sorting[i] = np.min(D)
    
    # Sort the triangles
    inds = (sorting).argsort()
    triangles = triangles[inds]

    # Construct the first node of the tree that will act as the root
    input_triangles = triangles
    tree = [dict() for x in range(input_triangles.shape[0])]
    tree[0]['vertices'] = triangles[0]
    tree[0]['predecessor'] = -1
    tree[0]['depth'] = 0
    tree[0]['index'] = 0
    tree[0]['adj_edge'] = np.array([])
    tree_index = 0

    # Initialize search
    input_triangles = np.delete(input_triangles,0,axis=0)
    stack = [tree[0]]

    # Build the tree by expanding nodes until the stack is empty
    # (The stack will be empty when the leaf nodes consist of only one edge)
    while len(stack) is not 0:
        # Pop the first element of the stack and delete it from the stack
        expanded_node = stack[0]
        del(stack[0])
        i = 0
        while i<input_triangles.shape[0]:
            # Construct two edge arrays: one for the parent and one for the candidate child
            # Orient the parent CCW as desired and the child CW to check for collisions
            polygon1_edges = np.array([np.vstack((expanded_node['vertices'][0], expanded_node['vertices'][1])), np.vstack((expanded_node['vertices'][1], expanded_node['vertices'][2])), np.vstack((expanded_node['vertices'][2], expanded_node['vertices'][0]))])
            polygon2_edges = np.array([np.vstack((input_triangles[i][0], input_triangles[i][2])), np.vstack((input_triangles[i][2], input_triangles[i][1])), np.vstack((input_triangles[i][1], input_triangles[i][0]))])
            triangles_touch = False
            for polygon1_edge_index in range(polygon1_edges.shape[0]):
                for polygon2_edge_index in range(polygon2_edges.shape[0]):
                    if (np.abs(polygon1_edges[polygon1_edge_index]-polygon2_edges[polygon2_edge_index])<1e-5).all():
                        triangles_touch = True
                        adj_edge_index = polygon2_edge_index

            # Check if the triangles touch, otherwise continue
            if not triangles_touch:
                i = i+1
                continue
            else:
                # Add the child to the tree
                tree_index = tree_index+1
                tree[tree_index]['predecessor'] = expanded_node['index']
                tree[tree_index]['depth'] = tree[tree[tree_index]['predecessor']]['depth']+1
                tree[tree_index]['index'] = tree_index
                tree[tree_index]['adj_edge'] = polygon2_edges[adj_edge_index]

                # Find the 3rd point of the child triangle (that does not belong to the shared edge) and arrange the vertices so that this is the 3rd point
                nrows, ncols = input_triangles[i].shape
                dtype = {'names':['f{}'.format(j) for j in range(ncols)], 'formats':ncols * [input_triangles[i].dtype]}
                set_diff = np.setdiff1d(input_triangles[i].view(dtype), np.array(tree[tree_index]['adj_edge'].transpose()).view(dtype))
                third_vertex = set_diff.view(input_triangles[i].dtype).reshape(-1, ncols)
                tree[tree_index]['vertices'] = np.array([tree[tree_index]['adj_edge'][1], tree[tree_index]['adj_edge'][0], third_vertex[0]]) # change the direction of adj_edge to make the child CCW again 

                # Delete the child from the input
                input_triangles = np.delete(input_triangles,i,axis=0)

                # Add the child to the stack to be expanded
                stack.append(tree[tree_index])

    # As a final step, sort the tree as a stack, in order of descending depth
    tree = sorted(tree, key=itemgetter('depth'), reverse=True)

    # Make sure to change the node and predecessor indices since the indices have now changed
    indices_new = np.zeros(len(tree))
    indices_old = np.zeros(len(tree))
    for i in range(0,len(tree)):
        indices_new[i] = i
        indices_old[i] = tree[i]['index']

    for i in range(0,len(tree)-1):
        tree[i]['predecessor'] = int(indices_new[indices_old==tree[i]['predecessor']][0])
        tree[i]['index'] = i
    tree[len(tree)-1]['index'] = len(tree)-1
    
    return tree

def polyconvexdecomposition(xy,workspace,touching_boundary):
    """
    Compute the convex decomposition of the input polygon and its dual (adjacency) graph.

    Input:
        xy                  : Vertex Coordinates of input polygon - start and end vertices must be the same
                              (Nx2 numpy.array)
        workspace           : Convex boundary of the workspace - start and end vertices must be the same
                              (Nx2 numpy.array)
        touching_boundary   : Flag that is True if the polygon is touching the boundary of the workspace and False otherwise
    Output:
        tree    : Array of dictionaries with polygons and generated adjacency graph
                  Each dictionary contains:
                    1) 'vertices': vertices of each polygon (arranged CCW - start and end vertices are NOT the same) - for each of the children, the vertices of the adjacency edge are the first two vertices
                    2) 'predecessor': the index of the polygon predecessor in the adjacency tree (-1 for the root)
                    3) 'depth': the depth of the (polygon) node in the adjacency tree (0 for the root)
                    4) 'index': the index of the polygon in the tree (its serial number)
                    5) 'adj_edge': the edge the polygon shares with its predecessor (CCW oriented with respect to the polygon)
    """
    # Construct a polygon based on the input coordinate vertices 
    polygon_in = Polygon(np.array(xy))

    # Find polygon vertices
    polygon_vertices = np.vstack((polygon_in.exterior.coords.xy[0][0:-1], polygon_in.exterior.coords.xy[1][0:-1])).transpose()

    # Find convex decomposition and delete problematic areas
    polygons = np.array([np.array(xi) for xi in polycvxdecomp(polygon_vertices.tolist())])
    idx = 0
    while idx <= len(polygons)-1:
        if Polygon(polygons[idx]).area < 0.01:
            polygons = np.delete(polygons, idx, axis=0)
        else:
            idx = idx+1

    # Sorting argument for polygons - Area if not touching boundary, Min distance to boundary if touching boundary
    sorting = np.zeros(polygons.shape[0])
    if not touching_boundary:
        for i in range(polygons.shape[0]):
            sorting[i] = -polyarea(polygons[i])
    elif touching_boundary:
        for i in range(polygons.shape[0]):
            if (Polygon(workspace).exterior).intersection(Polygon(polygons[i])).geom_type == 'LineString':
                sorting[i] = 0.0
            else:
                sorting[i] = 1.0
    
    # Sort the polygons
    inds = (sorting).argsort()
    polygons = polygons[inds]

    # Construct the first node of the tree that will act as the root
    input_polygons = polygons
    tree = [dict() for x in range(input_polygons.shape[0])]
    tree[0]['vertices'] = polygons[0]
    tree[0]['predecessor'] = -1
    tree[0]['depth'] = 0
    tree[0]['index'] = 0
    tree[0]['adj_edge'] = np.array([])
    tree_index = 0

    # Initialize search
    input_polygons = np.delete(input_polygons,0,axis=0)
    stack = [tree[0]]

    # Build the tree by expanding nodes until the stack is empty
    while len(stack) is not 0:
        # Pop the first element of the stack and delete it from the stack
        expanded_node = stack[0]
        del(stack[0])

        # Find edges of expanded node - CW
        polygon1_edges = []
        for j in range(0,expanded_node['vertices'].shape[0]):
            polygon1_edges.append(np.array([expanded_node['vertices'][(j+1)%expanded_node['vertices'].shape[0]], expanded_node['vertices'][j%expanded_node['vertices'].shape[0]]]))
        polygon1_edges = np.array(polygon1_edges)

        i = 0
        while i<input_polygons.shape[0]:
            # Find edges of candidate child - CCW
            polygon2_edges = []
            for j in range(0,input_polygons[i].shape[0]):
                polygon2_edges.append(np.array([input_polygons[i][j%input_polygons[i].shape[0]], input_polygons[i][(j+1)%input_polygons[i].shape[0]]]))
            polygon2_edges = np.array(polygon2_edges)

            polygons_touch = False
            for polygon1_edge_index in range(polygon1_edges.shape[0]):
                for polygon2_edge_index in range(polygon2_edges.shape[0]):
                    if (np.abs(polygon1_edges[polygon1_edge_index]-polygon2_edges[polygon2_edge_index])<1e-5).all():
                        polygons_touch = True
                        adj_edge_index = polygon2_edge_index

            # Check if the polygons touch, otherwise continue
            if not polygons_touch:
                i = i+1
                continue
            else:
                # Add the child to the tree with the adjacency edge being first
                tree_index = tree_index+1
                tree[tree_index]['predecessor'] = expanded_node['index']
                tree[tree_index]['depth'] = tree[tree[tree_index]['predecessor']]['depth']+1
                tree[tree_index]['index'] = tree_index
                tree[tree_index]['adj_edge'] = polygon2_edges[adj_edge_index]
                tree[tree_index]['vertices'] = np.roll(input_polygons[i],-adj_edge_index,axis=0)

                # As a final preprocessing step, check whether the edges before and after the adj_edge are parallel with adj_edge
                # If that's the case, cut the triangle corresponding to that edge as an extra polygon
                tangent_before = np.array(tree[tree_index]['vertices'][0]-tree[tree_index]['vertices'][-1])/np.linalg.norm(tree[tree_index]['vertices'][0]-tree[tree_index]['vertices'][-1])
                tangent_after = np.array(tree[tree_index]['vertices'][2]-tree[tree_index]['vertices'][1])/np.linalg.norm(tree[tree_index]['vertices'][2]-tree[tree_index]['vertices'][1])
                tangent_adj_edge = np.array(tree[tree_index]['adj_edge'][1]-tree[tree_index]['adj_edge'][0])/np.linalg.norm(tree[tree_index]['adj_edge'][1]-tree[tree_index]['adj_edge'][0])
                normal_adj_edge = np.array([-tangent_adj_edge[1],tangent_adj_edge[0]])
                if np.abs(np.dot(tangent_before,normal_adj_edge)) < 0.001:
                    # Add triangle
                    tree_before = dict()
                    tree_before['predecessor'] = tree_index
                    tree_before['depth'] = tree[tree_index]['depth']+1
                    tree_before['index'] = len(tree)
                    tree_before['adj_edge'] = np.vstack((tree[tree_index]['vertices'][0],tree[tree_index]['vertices'][-2]))
                    tree_before['vertices'] = np.vstack((tree[tree_index]['vertices'][0],tree[tree_index]['vertices'][-2],tree[tree_index]['vertices'][-1]))

                    # Delete the last vertex from the original polygon
                    tree[tree_index]['vertices'] = np.delete(tree[tree_index]['vertices'], -1, axis=0)

                    # Add the new triangle to the tree and stack
                    tree.append(tree_before)
                    stack.append(tree_before)
                
                if np.abs(np.dot(tangent_after,normal_adj_edge)) < 0.001:
                    # Add triangle
                    tree_after = dict()
                    tree_after['predecessor'] = tree_index
                    tree_after['depth'] = tree[tree_index]['depth']+1
                    tree_after['index'] = len(tree)
                    tree_after['adj_edge'] = np.vstack((tree[tree_index]['vertices'][3],tree[tree_index]['vertices'][1]))
                    tree_after['vertices'] = np.vstack((tree[tree_index]['vertices'][3],tree[tree_index]['vertices'][1],tree[tree_index]['vertices'][2]))

                    # Delete the third vertex from the original polygon
                    tree[tree_index]['vertices'] = np.delete(tree[tree_index]['vertices'], 2, axis=0)

                    # Add the new triangle to the tree and stack
                    tree.append(tree_after)
                    stack.append(tree_after)

                # Delete the child from the input
                input_polygons = np.delete(input_polygons,i,axis=0)

                # Add the child to the stack to be expanded
                stack.append(tree[tree_index])

    # As a final step, sort the tree as a stack, in order of descending depth
    tree = sorted(tree, key=itemgetter('depth'), reverse=True)

    # Make sure to change the node and predecessor indices since the indices have now changed
    indices_new = np.zeros(len(tree))
    indices_old = np.zeros(len(tree))
    for i in range(0,len(tree)):
        indices_new[i] = i
        indices_old[i] = tree[i]['index']

    for i in range(0,len(tree)-1):
        tree[i]['predecessor'] = int(indices_new[indices_old==tree[i]['predecessor']][0])
        tree[i]['index'] = i
    tree[len(tree)-1]['index'] = len(tree)-1
    
    return tree

def polyintersect(xy1,xy2):
    """
    Checks if polygon xy1 intersects polygon xy2

    Input:
        xy1     : Vertex Coordinates of a polygon
                  (Nx2 numpy.array)
        xy2     : Vertex Coordinates of a polygon
                  (Nx2 numpy.array)
    Output:
        outcome : True if xy1 intersects xy2, False otherwise
    """

    # Construct polygon objects based on the input vertices
    polygon1 = Polygon(xy1)
    polygon2 = Polygon(xy2)

    outcome = polygon1.intersects(polygon2)

    return outcome

def polyunion(xy1,xy2):
    """
    Computes the union of two polygons xy1 and xy2

    Input:
        xy1 : Vertex Coordinates of a polygon
              (Nx2 numpy.array)
        xy2 : Vertex Coordinates of a polygon
              (Nx2 numpy.array)
    Output:
        xy  : Coordinates of the output polygon
              (Nx2 numpy.array)
    """

    # Construct polygon objects based on the input vertices
    polygon1 = Polygon(xy1)
    polygon2 = Polygon(xy2)
    polygons = [polygon1,polygon2]

    # Find the union and orient appropriately
    output = cascaded_union(polygons)
    output = sp.geometry.polygon.orient(output, 1.0) # orient polygon to be CCW

    # Find the actual vertices
    xy = np.array(output.exterior.coords.xy).transpose()

    return xy

def lineint(l1, l2, precision=0):
    """Compute the intersection between two lines.

    Input:
        l1 : first line
        l2 : second line
        precision : precision to check if lines are parallel (default 0)

    Output:
        The intersection point
    """
    i = [0, 0] # point
    a1 = l1[1][1] - l1[0][1]
    b1 = l1[0][0] - l1[1][0]
    c1 = a1 * l1[0][0] + b1 * l1[0][1]
    a2 = l2[1][1] - l2[0][1]
    b2 = l2[0][0] - l2[1][0]
    c2 = a2 * l2[0][0] + b2 * l2[0][1]
    det = a1 * b2 - a2 * b1
    if not scalar_eq(det, 0, precision): # lines are not parallel
        i[0] = (b2 * c1 - b1 * c2) / det
        i[1] = (a1 * c2 - a2 * c1) / det
    return i

def linesegmentsintersect(p1, p2, q1, q2):
    """Checks if two line segments intersect.

    Input:
        p1 : The start vertex of the first line segment.
        p2 : The end vertex of the first line segment.
        q1 : The start vertex of the second line segment.
        q2 : The end vertex of the second line segment.

    Output:
        True if the two line segments intersect
    """
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    da = q2[0] - q1[0]
    db = q2[1] - q1[1]

    # segments are parallel
    if (da*dy - db*dx) == 0:
        return False

    s = (dx * (q1[1] - p1[1]) + dy * (p1[0] - q1[0])) / (da * dy - db * dx)
    t = (da * (p1[1] - q1[1]) + db * (q1[0] - p1[0])) / (db * dx - da * dy)

    return s >= 0 and s <= 1 and t >= 0 and t <= 1

def trianglearea(a, b, c):
    """Calculates the area of a triangle spanned by three points.
    Note that the area will be negative if the points are not given in counter-clockwise order.

    Input:
        a : First point
        b : Second point
        c : Third point

    Output:
        Area of triangle
    """
    return ((b[0] - a[0])*(c[1] - a[1]))-((c[0] - a[0])*(b[1] - a[1]))

def isleft(a, b, c):
    return trianglearea(a, b, c) > 0

def islefton(a, b, c):
    return trianglearea(a, b, c) >= 0

def isright(a, b, c):
    return trianglearea(a, b, c) < 0

def isrighton(a, b, c):
    return trianglearea(a, b, c) <= 0

def collinear(a, b, c, thresholdAngle=0):
    """Checks if three points are collinear.

    Input:
        a : First point
        b : Second point
        c : Third point
        thresholdAngle : threshold to consider if points are collinear, in radians (default 0)

    Output:
        True if points are collinear
    """
    if thresholdAngle == 0:
        return trianglearea(a, b, c) == 0
    else:
        ab = [None] * 2
        bc = [None] * 2

        ab[0] = b[0]-a[0]
        ab[1] = b[1]-a[1]
        bc[0] = c[0]-b[0]
        bc[1] = c[1]-b[1]

        dot = ab[0]*bc[0] + ab[1]*bc[1]
        magA = math.sqrt(ab[0]*ab[0] + ab[1]*ab[1])
        magB = math.sqrt(bc[0]*bc[0] + bc[1]*bc[1])
        angle = math.acos(dot/(magA*magB))
        return angle < thresholdAngle

def sqdist(a, b):
    dx = b[0] - a[0]
    dy = b[1] - a[1]
    return dx * dx + dy * dy

def polyat(polygon, i):
    """Gets a vertex at position i on the polygon.
    It does not matter if i is out of bounds.

    Input:
        polygon : The polygon
        i : Position desired on the polygon

    Output:
        Vertex at position i
    """
    s = len(polygon)
    return polygon[i % s]

def polyclear(polygon):
    """Clears the polygon data

    Input:
        polygon : The polygon
    """
    del polygon[:]

def polyappend(polygon, poly, start, end):
    """Grabs points at indicies `start` to `end` from `poly`
    and appends them to `polygon`

    Input:
        polygon : The destination polygon
        poly : The source polygon
        start : Starting source index
        end : Ending source index (not included in the slice)
    """
    for i in range(start, end):
        polygon.append(poly[i])

def polymakeccw(polygon):
    """Makes sure that the polygon vertices are ordered counter-clockwise.

    Input:
        polygon : The polygon
    """
    br = 0
    v = polygon

    # find bottom right point
    for i in range(1, len(polygon)):
        if v[i][1] < v[br][1] or (v[i][1] == v[br][1] and v[i][0] > v[br][0]):
            br = i

    # reverse poly if clockwise
    if not isleft(polyat(polygon, br - 1), polyat(polygon, br), polyat(polygon, br + 1)):
        polyreverse(polygon)

def polyreverse(polygon):
    """Reverses the vertices in the polygon.

    Input:
        polygon : The polygon
    """
    polygon.reverse()

def polyisreflex(polygon, i):
    """Checks if a point in the polygon is a reflex point.

    Input:
        polygon : The polygon
        i : index of point to check
    
    Output:
        True is point is a reflex point
    """
    return isright(polyat(polygon, i - 1), polyat(polygon, i), polyat(polygon, i + 1))

def polycansee(polygon, a, b):
    """Checks if two vertices in the polygon can see each other.

    Input:
        polygon : The polygon
        a : Vertex 1
        b : Vertex 2

    Output:
        True if vertices can see each other
    """

    l1 = [None]*2
    l2 = [None]*2

    if islefton(polyat(polygon, a + 1), polyat(polygon, a), polyat(polygon, b)) and isrighton(polyat(polygon, a - 1), polyat(polygon, a), polyat(polygon, b)):
        return False

    dist = sqdist(polyat(polygon, a), polyat(polygon, b))
    for i in range(0, len(polygon)): # for each edge
        if (i + 1) % len(polygon) == a or i == a: # ignore incident edges
            continue

        if islefton(polyat(polygon, a), polyat(polygon, b), polyat(polygon, i + 1)) and isrighton(polyat(polygon, a), polyat(polygon, b), polyat(polygon, i)): # if diag intersects an edge
            l1[0] = polyat(polygon, a)
            l1[1] = polyat(polygon, b)
            l2[0] = polyat(polygon, i)
            l2[1] = polyat(polygon, i + 1)
            p = lineint(l1, l2)
            if sqdist(polyat(polygon, a), p) < dist: # if edge is blocking visibility to b
                return False

    return True

def polycopy(polygon, i, j, targetPoly=None):
    """Copies the polygon from vertex i to vertex j to targetPoly.

    Input:
        polygon : The source polygon
        i : start vertex
        j : end vertex (inclusive)
        targetPoly -- Optional target polygon

    Output:
        The resulting copy.
    """
    p = targetPoly or []
    polyclear(p)
    if i < j:
        # Insert all vertices from i to j
        for k in range(i, j+1):
            p.append(polygon[k])

    else:
        # Insert vertices 0 to j
        for k in range(0, j+1):
            p.append(polygon[k])

        # Insert vertices i to end
        for k in range(i, len(polygon)):
            p.append(polygon[k])

    return p

def polygetcutedges(polygon):
    """Decomposes the polygon into convex pieces.
    Note that this algorithm has complexity O(N^4) and will be very slow for polygons with many vertices.

    Input:
        polygon : The polygon

    Output:
        A list of edges [[p1,p2],[p2,p3],...] that cut the polygon.
    """
    mins = []
    tmp1 = []
    tmp2 = []
    tmpPoly = []
    nDiags = float('inf')

    for i in range(0, len(polygon)):
        if polyisreflex(polygon, i):
            for j in range(0, len(polygon)):
                if polycansee(polygon, i, j):
                    tmp1 = polygetcutedges(polycopy(polygon, i, j, tmpPoly))
                    tmp2 = polygetcutedges(polycopy(polygon, j, i, tmpPoly))

                    for k in range(0, len(tmp2)):
                        tmp1.append(tmp2[k])

                    if len(tmp1) < nDiags:
                        mins = tmp1
                        nDiags = len(tmp1)
                        mins.append([polyat(polygon, i), polyat(polygon, j)])

    return mins

def polydecomp(polygon):
    """Decomposes the polygon into one or more convex sub-polygons.

    Input:
        polygon : The polygon

    Output:
        An array or polygon objects.
    """
    edges = polygetcutedges(polygon)
    if len(edges) > 0:
        return polyslice(polygon, edges)
    else:
        return [polygon]

def polyslice(polygon, cutEdges):
    """Slices the polygon given one or more cut edges. If given one, this function will return two polygons (false on failure). If many, an array of polygons.
    
    Input:
        polygon : The polygon
        cutEdges : A list of edges to cut on, as returned by getCutEdges()
    
    Output:
        An array of polygon objects.
    """
    if len(cutEdges) == 0:
        return [polygon]

    if isinstance(cutEdges, list) and len(cutEdges) != 0 and isinstance(cutEdges[0], list) and len(cutEdges[0]) == 2 and isinstance(cutEdges[0][0], list):

        polys = [polygon]

        for i in range(0, len(cutEdges)):
            cutEdge = cutEdges[i]
            # Cut all polys
            for j in range(0, len(polys)):
                poly = polys[j]
                result = polyslice(poly, cutEdge)
                if result:
                    # Found poly! Cut and quit
                    del polys[j:j+1]
                    polys.extend((result[0], result[1]))
                    break

        return polys
    else:

        # Was given one edge
        cutEdge = cutEdges
        i = polygon.index(cutEdge[0])
        j = polygon.index(cutEdge[1])

        if i != -1 and j != -1:
            return [polycopy(polygon, i, j),
                    polycopy(polygon, j, i)]
        else:
            return False

def polyissimple(polygon):
    """Checks that the line segments of this polygon do not intersect each other.
    
    Input:
        polygon : The polygon
    
    Output:
        True is polygon is simple (not self-intersecting)
    
    Todo:
        Should it check all segments with all others?
    """
    path = polygon
    # Check
    for i in range(0,len(path)-1):
        for j in range(0, i-1):
            if linesegmentsintersect(path[i], path[i+1], path[j], path[j+1]):
                return False

    # Check the segment between the last and the first point to all others
    for i in range(1,len(path)-2):
        if linesegmentsintersect(path[0], path[len(path)-1], path[i], path[i+1]):
            return False

    return True

def getintersection(p1, p2, q1, q2, delta=0):
    """Gets the intersection point 
    
    Input:
        p1 : The start vertex of the first line segment.
        p2 : The end vertex of the first line segment.
        q1 : The start vertex of the second line segment.
        q2 : The end vertex of the second line segment.
        delta : Optional precision to check if lines are parallel (default 0)
    
    Output:
        The intersection point.
    """
    a1 = p2[1] - p1[1]
    b1 = p1[0] - p2[0]
    c1 = (a1 * p1[0]) + (b1 * p1[1])
    a2 = q2[1] - q1[1]
    b2 = q1[0] - q2[0]
    c2 = (a2 * q1[0]) + (b2 * q1[1])
    det = (a1 * b2) - (a2 * b1)

    if not scalar_eq(det, 0, delta):
        return [((b2 * c1) - (b1 * c2)) / det, ((a1 * c2) - (a2 * c1)) / det]
    else:
        return [0, 0]

def polycvxdecomp(polygon, result=None, reflexVertices=None, steinerPoints=None, delta=25, maxlevel=10000, level=0):
    """Quickly decompose the Polygon into convex sub-polygons. Algorithm based on Mark Bayazit's polygon decomposition.
    
    Input:
        polygon : The polygon to decompose
        result : Stores result of decomposed polygon, passed recursively
        reflexVertices : 
        steinerPoints :
        delta : Currently unused
        maxlevel : The maximum allowed level of recursion
        level : The current level of recursion
    
    Output:
        List of decomposed convex polygons
    """
    if result is None:
        result = []
    reflexVertices = reflexVertices or []
    steinerPoints = steinerPoints or []

    upperInt = [0, 0]
    lowerInt = [0, 0]
    p = [0, 0]         # Points
    upperDist = 0
    lowerDist = 0
    d = 0
    closestDist = 0 # scalars
    upperIndex = 0
    lowerIndex = 0
    closestIndex = 0 # integers
    lowerPoly = []
    upperPoly = [] # polygons
    poly = polygon
    v = polygon

    if len(v) < 3:
        return result

    level += 1
    if level > maxlevel:
        print("quickDecomp: max level ("+str(maxlevel)+") reached.")
        return result

    for i in range(0, len(polygon)):
        if polyisreflex(poly, i):
            reflexVertices.append(poly[i])
            upperDist = float('inf')
            lowerDist = float('inf')

            for j in range(0, len(polygon)):
                if isleft(polyat(poly, i - 1), polyat(poly, i), polyat(poly, j)) and isrighton(polyat(poly, i - 1), polyat(poly, i), polyat(poly, j - 1)): # if line intersects with an edge
                    p = getintersection(polyat(poly, i - 1), polyat(poly, i), polyat(poly, j), polyat(poly, j - 1)) # find the point of intersection
                    if isright(polyat(poly, i + 1), polyat(poly, i), p): # make sure it's inside the poly
                        d = sqdist(poly[i], p)
                        if d < lowerDist: # keep only the closest intersection
                            lowerDist = d
                            lowerInt = p
                            lowerIndex = j

                if isleft(polyat(poly, i + 1), polyat(poly, i), polyat(poly, j + 1)) and isrighton(polyat(poly, i + 1), polyat(poly, i), polyat(poly, j)):
                    p = getintersection(polyat(poly, i + 1), polyat(poly, i), polyat(poly, j), polyat(poly, j + 1))
                    if isleft(polyat(poly, i - 1), polyat(poly, i), p):
                        d = sqdist(poly[i], p)
                        if d < upperDist:
                            upperDist = d
                            upperInt = p
                            upperIndex = j

            # if there are no vertices to connect to, choose a point in the middle
            if lowerIndex == (upperIndex + 1) % len(polygon):
                #print("Case 1: Vertex("+str(i)+"), lowerIndex("+str(lowerIndex)+"), upperIndex("+str(upperIndex)+"), poly.size("+str(len(polygon))+")")
                p[0] = (lowerInt[0] + upperInt[0]) / 2
                p[1] = (lowerInt[1] + upperInt[1]) / 2
                steinerPoints.append(p)

                if i < upperIndex:
                    #lowerPoly.insert(lowerPoly.end(), poly.begin() + i, poly.begin() + upperIndex + 1)
                    polyappend(lowerPoly, poly, i, upperIndex+1)
                    lowerPoly.append(p)
                    upperPoly.append(p)
                    if lowerIndex != 0:
                        #upperPoly.insert(upperPoly.end(), poly.begin() + lowerIndex, poly.end())
                        polyappend(upperPoly, poly, lowerIndex, len(poly))

                    #upperPoly.insert(upperPoly.end(), poly.begin(), poly.begin() + i + 1)
                    polyappend(upperPoly, poly, 0, i+1)
                else:
                    if i != 0:
                        #lowerPoly.insert(lowerPoly.end(), poly.begin() + i, poly.end())
                        polyappend(lowerPoly, poly, i, len(poly))

                    #lowerPoly.insert(lowerPoly.end(), poly.begin(), poly.begin() + upperIndex + 1)
                    polyappend(lowerPoly, poly, 0, upperIndex+1)
                    lowerPoly.append(p)
                    upperPoly.append(p)
                    #upperPoly.insert(upperPoly.end(), poly.begin() + lowerIndex, poly.begin() + i + 1)
                    polyappend(upperPoly, poly, lowerIndex, i+1)

            else:
                # connect to the closest point within the triangle
                #print("Case 2: Vertex("+str(i)+"), closestIndex("+str(closestIndex)+"), poly.size("+str(len(polygon))+")\n")

                if lowerIndex > upperIndex:
                    upperIndex += len(polygon)

                closestDist = float('inf')

                if upperIndex < lowerIndex:
                    return result

                for j in range(lowerIndex, upperIndex+1):
                    if islefton(polyat(poly, i - 1), polyat(poly, i), polyat(poly, j)) and isrighton(polyat(poly, i + 1), polyat(poly, i), polyat(poly, j)):
                        d = sqdist(polyat(poly, i), polyat(poly, j))
                        if d < closestDist:
                            closestDist = d
                            closestIndex = j % len(polygon)

                if i < closestIndex:
                    polyappend(lowerPoly, poly, i, closestIndex+1)
                    if closestIndex != 0:
                        polyappend(upperPoly, poly, closestIndex, len(v))

                    polyappend(upperPoly, poly, 0, i+1)
                else:
                    if i != 0:
                        polyappend(lowerPoly, poly, i, len(v))

                    polyappend(lowerPoly, poly, 0, closestIndex+1)
                    polyappend(upperPoly, poly, closestIndex, i+1)

            # solve smallest poly first
            if len(lowerPoly) < len(upperPoly):
                polycvxdecomp(lowerPoly, result, reflexVertices, steinerPoints, delta, maxlevel, level)
                polycvxdecomp(upperPoly, result, reflexVertices, steinerPoints, delta, maxlevel, level)
            else:
                polycvxdecomp(upperPoly, result, reflexVertices, steinerPoints, delta, maxlevel, level)
                polycvxdecomp(lowerPoly, result, reflexVertices, steinerPoints, delta, maxlevel, level)

            return result

    result.append(polygon)

    return result

def polyremovecollinear(polygon, precision=0):
    """Remove collinear points in the polygon.
    
    Input:
        polygon : The polygon
        precision : The threshold angle to use when determining whether two edges are collinear. (default is 0)
    
    Output:
        The number of points removed
    """
    num = 0
    i = len(polygon) - 1
    while len(polygon) > 3 and i >= 0:
    #(var i=polygon.length-1; polygon.length>3 && i>=0; --i){
        if collinear(polyat(polygon, i - 1), polyat(polygon, i), polyat(polygon, i+1), precision):
            # Remove the middle point
            del polygon[i % len(polygon):(i % len(polygon))+1]
            num += 1
        i -= 1
    return num

def scalar_eq(a, b, precision=0):
    """Check if two scalars are equal.
    
    Input:
        a : first scalar
        b : second scalar
        precision : precision to check equality
    
    Output:
        True if scalars are equal
    """
    return abs(a - b) <= precision