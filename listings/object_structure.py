#!/usr/bin/python3
import pprint
from voronoi import voronoi

generators = voronoi.Generators()
generators.randomize(num=5)

voronoi_tessellation = voronoi.VoronoiTessellation()
voronoi_tessellation.compute_from_generators(generators)

pp = pprint.PrettyPrinter(width=60)

print("voronoi_tessellation.nodes:")
pp.pprint(voronoi_tessellation.nodes)
print()

print("voronoi_tessellation.edges:")
pp.pprint(voronoi_tessellation.edges)
print()

print("voronoi_tessellation.faces:")
pp.pprint(voronoi_tessellation.faces)
