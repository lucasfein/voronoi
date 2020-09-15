#!/usr/bin/python3
import matplotlib.pyplot as plt
import matplotlib
from voronoi import voronoi
import numpy as np

generators = voronoi.Generators()
generators.randomize(num=15, weight_range=(10.0, 10.0), min_distance=1.0)

voronoi_tessellation = voronoi.VoronoiTessellation()
voronoi_tessellation.compute_from_generators(generators)
voronoi_tessellation.plot(include_generators=False,
                          edge_color="gray",
                          edge_style="dashed")

generators.randomize_weights(weight_range=(10.0, 100.0), min_distance=0.25)
voronoi_tessellation.compute_from_generators(generators)
voronoi_tessellation.plot(generator_color="RdBu_r")

plt.show()
