#!/bin/env python3
import numpy as np
import pypoman
from scipy.spatial import ConvexHull

cw = 40 # width of top square
ch = 10 # height of top cube
ss = 4 # width/length of bottom square
# Half side lengths
cwh = cw/2
chh = ch/2
shs = ss/2

# Cube vertices
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

# vertices for a 3D parallelogram
# top square is of half-side cwh
# bottom square is of half-side shs
# height is of half-side hL
L = 30  # height of parallelogram
hL = L/2 # half height
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

# shift down to match the bottom surface of the cube
vertex_list_shifted_parallelogram = [[vertex[0], vertex[1], vertex[2]-chh-hL] for vertex in vertex_list_parallelogram ]
vertex_list_concat = vertex_list_cube + vertex_list_shifted_parallelogram

# compute convex hull
hull = ConvexHull(vertex_list_concat)
# convex hull vertices are guaranteed to be in counter clocwise order
vertex_list = [vertex_list_concat[vertex] for vertex in hull.vertices]

A, b = pypoman.compute_polytope_halfspaces(np.array(vertex_list))

vertices_computed = pypoman.compute_polytope_vertices(A, b)
result = pypoman.compute_chebyshev_center(A, b)
center = result[0:3]
radius = result[3]

print("Top square side is {} / Bottom square side is {}".format(cwh, shs))
print("Top cube height = {} / Parallelogram height = {} / Total height = {}".format(ch, L, L + ch))
print("Chebyshev Center is :", center)
print("Chebyshev Radius is :", radius)


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
def plot(hull, vertices, center, radius):
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

    for s in hull.simplices:
        s = np.append(s, s[0])  # Here we cycle back to the first coordinate
        v = np.array(vertices)
        ax.plot(v[s, 0], v[s, 1], v[s, 2], "b-")

    plt.title("Chebyshev Center: " + str(np.array(center).round(3)) + " Radius: " + str(round(radius,3)))
    set_aspect_equal_3d(ax)
    plt.show()

plot(hull, vertex_list, center, radius)
