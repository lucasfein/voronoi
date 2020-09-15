#!/usr/bin/python3
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from voronoi import voronoi
import progressbar

num_generators = 25
error_delta_dist = []
min_median_dist = []
sample_size = 100
generators = voronoi.Generators()
voronoi_tessellation = voronoi.VoronoiTessellation()
bar = progressbar.ProgressBar(max_value=sample_size)

for i in range(sample_size):
    generators.randomize(num=num_generators)
    voronoi_tessellation.compute_from_generators(generators)

    (approximate_generators, _,
     _) = voronoi_tessellation.approximate_generators()

    centroid_distribution = voronoi_tessellation.face_centroids().distance(
        generators)
    approximation_distribution = approximate_generators.distance(generators)

    centroid_distribution_median = np.median(centroid_distribution)
    approximation_distribution_median = np.median(approximation_distribution)

    delta = approximation_distribution_median - centroid_distribution_median
    error_delta_dist.append(delta)
    min_median_dist.append(
        min(centroid_distribution_median, approximation_distribution_median))

    bar.update(i + 1)

mean = np.mean(error_delta_dist)
std = np.std(error_delta_dist)

plt.hist(error_delta_dist, bins=25)
plt.title("n = {}, µ' = {:.4f}, σ = {:.4f}".format(sample_size, mean, std))
plt.ylabel("number of approximations")
plt.xlabel("Δµ")
plt.axvline(mean, color='k', linestyle='dashed', linewidth=1)
plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))

plt.show()

mean = np.mean(min_median_dist)
std = np.std(min_median_dist)

plt.hist(min_median_dist, bins=25)
plt.title("n = {}, µ' = {:.4f}, σ = {:.4f}".format(sample_size, mean, std))
plt.ylabel("number of approximations")
plt.xlabel("µ")
plt.axvline(mean, color='k', linestyle='dashed', linewidth=1)
plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))

plt.show()
