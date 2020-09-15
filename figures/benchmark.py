#!/usr/bin/python3
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from voronoi import voronoi
import progressbar
import datetime

sample_size = 100
intervals = (10, 101, 10)
bar = progressbar.ProgressBar(max_value=len(range(*intervals)) * sample_size)
bar.update(0)
min_means = []
run_time_voronoi = []
run_time_delaunay = []

iteration = 0
for i in range(*intervals):
    temp_run_time_voronoi = []
    temp_run_time_delaunay = []
    for j in range(sample_size):
        generators = voronoi.Generators()
        generators.randomize(num=i)

        t_0 = datetime.datetime.now()
        voronoi_tessellation = voronoi.VoronoiTessellation()
        voronoi_tessellation.compute_from_generators(generators)
        t_1 = datetime.datetime.now()

        delta = t_1 - t_0
        temp_run_time_voronoi.append(delta.total_seconds() * 1000)

        t_0 = datetime.datetime.now()
        delaunay_triangulation = voronoi.DelaunayTriangulation()
        delaunay_triangulation.compute_from_generators(generators)
        t_1 = datetime.datetime.now()

        delta = t_1 - t_0
        temp_run_time_delaunay.append(delta.total_seconds() * 1000)

        iteration += 1
        bar.update(iteration)

    run_time_voronoi.append(np.median(temp_run_time_voronoi))
    run_time_delaunay.append(np.median(temp_run_time_delaunay))

plt.plot(np.arange(*intervals), run_time_voronoi, label="Voronoi Tessellation")
plt.plot(np.arange(*intervals),
         run_time_delaunay,
         label="Delaunay Triangulation")
plt.xticks(np.arange(*intervals))

plt.title("n = {}".format(sample_size))
plt.ylabel("run time [ms]")
plt.xlabel("|G|")

plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.legend()

plt.show()
