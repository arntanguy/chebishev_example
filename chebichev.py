#!/bin/env python3
import numpy as np
import pypoman
from scipy.spatial import ConvexHull

# generate vertices for a 3D cube of side length l
cw = 40
ch = 10
ss = 4
cwh = cw/2 # large half side (top of pyramid)
chh = ch/2
shs = ss/2 # small half side (bottom of pyramid)

vertex_list_cube = [
    [ -cwh,-cwh,-chh ],
    [ cwh,-cwh,-chh ],
    [ cwh,-cwh, chh ],
    [ -cwh,-cwh, chh ],
    [ -cwh, cwh,-chh ],
    [ cwh, cwh,-chh ],
    [ cwh, cwh, chh ],
    [ -cwh, cwh, chh ]
    ]

# generate vertices for a 3D parallelogram
# top square is of half-side lhs
# bottom square is of half-side shs
# height is of half-side hL
L = 30 
hL = L/2
vertex_list_parallelogram = [
    [ -shs,-shs,-hL ],
    [ shs,-shs,-hL ],
    [ shs, shs,-hL ],
    [ -shs, shs,-hL ],
    [ -cwh,-cwh, hL ],
    [ cwh,-cwh, hL ],
    [ cwh, cwh, hL ],
    [ -cwh, cwh, hL ]
]

# shift down by half of the top cube
vertex_list_shifted_parallelogram = [[vertex[0], vertex[1], vertex[2]-chh-hL] for vertex in vertex_list_parallelogram ]

vertex_list_concat = vertex_list_cube + vertex_list_shifted_parallelogram
print(vertex_list_parallelogram)
print(vertex_list_shifted_parallelogram)

# compute convex hull
hull = ConvexHull(vertex_list_concat)
# vertices are guaranteed to be in counter clocwise order
print(hull.vertices)
vertex_list = [vertex_list_concat[vertex] for vertex in hull.vertices]
vertex_list = vertex_list_concat

print(vertex_list)

vertices = map(
    np.array,
    vertex_list
)
type(vertices)

A, b = pypoman.compute_polytope_halfspaces(vertices)
vertices_computed = pypoman.compute_polytope_vertices(A, b)
print(vertex_list)
print(vertices_computed)
result = pypoman.compute_chebyshev_center(A, b)
center = result[0:3]
radius = result[3]

print("Center is :", center)
print("Radius is :", radius)


def set_aspect_equal_3d(ax):
    """Fix equal aspect bug for 3D plots."""

    xlim = ax.get_xlim3d()
    ylim = ax.get_ylim3d()
    zlim = ax.get_zlim3d()

    from numpy import mean
    xmean = mean(xlim)
    ymean = mean(ylim)
    zmean = mean(zlim)

    plot_radius = max([abs(lim - mean_)
                       for lims, mean_ in ((xlim, xmean),
                                           (ylim, ymean),
                                           (zlim, zmean))
                       for lim in lims])

    ax.set_xlim3d([xmean - plot_radius, xmean + plot_radius])
    ax.set_ylim3d([ymean - plot_radius, ymean + plot_radius])
    ax.set_zlim3d([zmean - plot_radius, zmean + plot_radius])

# Plot a sphere of radius "radius" centered at center
# Plot all 3D vertices as points
def plot(vertices, center, radius):

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for vertex in vertices:
        ax.scatter(vertex[0], vertex[1], vertex[2])

    u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:20j]
    x = center[0] + radius*np.cos(u)*np.sin(v)
    y = center[1] + radius*np.sin(u)*np.sin(v)
    z = center[2] + radius*np.cos(v)
    ax.plot_wireframe(x, y, z)

    ax.scatter(center[0], center[1], center[2], color='r')

    # plot a vector from center of length radius
    ax.quiver(center[0], center[1], center[2], radius, 0, 0, color='r')

    set_aspect_equal_3d(ax)


    plt.show()

plot(vertex_list, center, radius) 
