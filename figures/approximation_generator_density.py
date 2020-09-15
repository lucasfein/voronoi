#!/usr/bin/python3
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from voronoi import voronoi
import progressbar
import datetime

generators = voronoi.Generators()
voronoi_tessellation = voronoi.VoronoiTessellation()

sample_size = 25
intervals = (10, 101, 10)
num_intervals = len(range(*intervals))
bar = progressbar.ProgressBar(max_value=num_intervals * sample_size)
bar.update(0)
min_means = []
centroid_means = []
run_time = []

iteration = 0
for i in range(*intervals):
    temp_min_medians = []
    temp_centroid_means = []
    temp_run_time = []
    for j in range(sample_size):
        generators.randomize(num=i)
        voronoi_tessellation.compute_from_generators(generators)

        t_0 = datetime.datetime.now()
        (approximate_generators, _,
         _) = voronoi_tessellation.approximate_generators()
        t_1 = datetime.datetime.now()

        delta = t_1 - t_0

        temp_run_time.append(delta.total_seconds())

        centroid_distribution = voronoi_tessellation.face_centroids().distance(
            generators)
        approximation_distribution = approximate_generators.distance(
            generators)

        temp_min_medians.append(
            min(np.median(centroid_distribution),
                np.median(approximation_distribution)))
        temp_centroid_means.append(np.median(centroid_distribution))

        iteration += 1
        bar.update(iteration)

    min_means.append(np.mean(temp_min_medians))
    centroid_means.append(np.mean(temp_centroid_means))
    run_time.append(np.median(temp_run_time))

fig, ax1 = plt.subplots()
plt.title("n = {}".format(sample_size))
ax1.plot(np.arange(*intervals), min_means, label="minimal error")
ax1.plot(np.arange(*intervals),
         centroid_means,
         label="distance from centroids")
ax1.set_xticks(np.arange(*intervals))

ax1.set_ylabel("error")
ax1.set_xlabel("|G|")
ax1.legend(loc="upper left")

ax2 = ax1.twinx()
ax2.plot(np.arange(*intervals), run_time, color="gray", label="run time")
ax2.set_ylabel("run time [s]")
ax2.set_yscale("log")
ax2.legend(loc="upper right")

plt.show()
