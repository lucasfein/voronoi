#!/usr/bin/python3
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from voronoi import voronoi

num_generators = 25
generators = voronoi.Generators()
generators.randomize(num=num_generators, weight_range=(5.0, 5.0))

voronoi_tessellation = voronoi.VoronoiTessellation()
voronoi_tessellation.compute_from_generators(generators)

centroids = voronoi_tessellation.face_centroids(weight=5.0)

centroids_voronoi_tessellation = voronoi.VoronoiTessellation()
centroids_voronoi_tessellation.compute_from_generators(centroids)

voronoi_tessellation.plot()
centroids_voronoi_tessellation.plot(generator_color="red",
                                    edge_color="red",
                                    edge_style="dashed")

plt.show()

centroid_distribution = centroids.distance(generators)
centroid_distribution_median = np.median(centroid_distribution)

plt.hist(centroid_distribution, bins=25, label="centroids")
plt.legend()
plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
plt.xlabel("distance to nearest generator")
plt.ylabel("number of generators")
plt.title("n = {}, µ = {:.4f}".format(num_generators,
                                      centroid_distribution_median))
plt.axvline(centroid_distribution_median,
            color='k',
            linestyle='dashed',
            linewidth=1)

plt.show()

(approximate_generators, accepted_errors,
 errors) = voronoi_tessellation.approximate_generators()

voronoi_tessellation.plot()
approximate_generators.plot(color="red")
approximate_voronoi = voronoi.VoronoiTessellation()
approximate_voronoi.compute_from_generators(approximate_generators)
approximate_voronoi.plot(include_generators=False,
                         edge_color="red",
                         edge_style="dashed")
plt.show()

plt.plot(errors)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().set_ylim(errors[-1] * 0.9, errors[0] * 1.1)
plt.axhline(errors[0], color="k", linestyle="dashed", linewidth=1)
plt.axhline(accepted_errors[-1, :][1],
            color="k",
            linestyle="dashed",
            linewidth=1)

plt.xlabel("iteration")
plt.ylabel("minimal mean error of vertices")

plt.show()

approximation_distribution = approximate_generators.distance(generators)
approximation_distribution_median = np.median(approximation_distribution)

delta = approximation_distribution_median - centroid_distribution_median
plt.hist([centroid_distribution, approximation_distribution],
         label=["centroids", "approximation"],
         bins=25)
plt.legend()
plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
plt.xlabel("distance to nearest generator")
plt.ylabel("number of generators")
plt.title("n = {}, Δµ = {:.4f}".format(num_generators, delta))
plt.axvline(centroid_distribution_median,
            color='k',
            linestyle='dashed',
            linewidth=1)
plt.axvline(approximation_distribution_median,
            color='k',
            linestyle='dashed',
            linewidth=1)

arrow_y_coord = 0.33 * (
    max(np.histogram(centroid_distribution, bins=25)[0]) +
    max(np.histogram(approximation_distribution, bins=25)[0]))

plt.arrow(centroid_distribution_median,
          arrow_y_coord,
          delta,
          0,
          length_includes_head=True,
          overhang=0.0,
          color="k",
          head_width=0.025 * abs(delta),
          head_length=0.25 * abs(delta))

plt.show()
