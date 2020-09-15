#!/usr/bin/python3
import matplotlib.pyplot as plt
from voronoi import voronoi

generators = voronoi.Generators()
generators.randomize(num=25, weight_range=(5.0, 5.0))

delaunay_triangulation = voronoi.DelaunayTriangulation()
delaunay_triangulation.compute_from_generators(generators)

voronoi_tessellation = voronoi.VoronoiTessellation()
voronoi_tessellation.compute_from_generators(generators)

delaunay_triangulation.plot(include_generators=False,
                            edge_color="gray",
                            edge_style="dashed")
voronoi_tessellation.plot()

plt.show()
