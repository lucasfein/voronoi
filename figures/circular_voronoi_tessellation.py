#!/usr/bin/python3
import matplotlib.pyplot as plt
from voronoi import voronoi

generators = voronoi.Generators()
generators.randomize(num=10, weight_range=(25.0, 30.0), min_distance=0.5)

delaunay_triangulation = voronoi.DelaunayTriangulation()
delaunay_triangulation.compute_from_generators(generators)

voronoi_tessellation = voronoi.CircularVoronoiTessellation()
voronoi_tessellation.compute_from_generators(generators)

delaunay_triangulation.plot(include_generators=False,
                            edge_color="gray",
                            edge_style="dashed")

voronoi_tessellation.plot()
plt.show()