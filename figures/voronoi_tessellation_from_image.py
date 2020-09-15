#!/usr/bin/python3
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from voronoi import voronoi

voronoi_tessellation = voronoi.VoronoiTessellation()
voronoi_tessellation.compute_from_image("data/2.jpg")

img = mpimg.imread("data/2.jpg")
plt.imshow(img)

voronoi_tessellation.plot(edge_color="red", edge_width=2.0)

plt.tight_layout()
plt.show()