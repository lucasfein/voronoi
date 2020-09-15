#!/usr/bin/python3
import matplotlib.pyplot as plt
from voronoi import voronoi

generators = voronoi.Generators()
generators.randomize(num=50, weight_range=(5.0, 5.0))

voronoi_tessellation = voronoi.VoronoiTessellation()
voronoi_tessellation.compute_from_generators(generators)

voronoi_tessellation.plot()
plt.show()