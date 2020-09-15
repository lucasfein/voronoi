#!/usr/bin/python3
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from voronoi import voronoi

generators = voronoi.Generators()
generators.compute_from_image("data/1.jpg")

voronoi_tessellation = voronoi.VoronoiTessellation()
voronoi_tessellation.compute_from_generators(generators)

img = mpimg.imread("data/1.jpg")
plt.imshow(img)

voronoi_tessellation.plot(edge_color="white",
                          generator_color=(1.0, 1.0, 1.0, 0.5))

plt.tight_layout()
plt.show()
