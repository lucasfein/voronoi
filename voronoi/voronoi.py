# -*- coding: utf-8 -*-
"""
voronoi
-------

Voronoi Tessellations of weighted generators in the Euclidean plane
"""
import os
import csv
import uuid
import subprocess
import warnings
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import scipy.ndimage
import skimage
import skimage.feature
import skimage.filters
import skimage.morphology
import skimage.segmentation
import skimage.measure
import skimage.exposure
from lxml import etree

warnings.filterwarnings("ignore", module="networkx")
warnings.filterwarnings("ignore", module="matplotlib")
warnings.filterwarnings("ignore", module="skimage")

DPI = 1200


class Generators():
    """Generators"""
    def __init__(self):
        """
        Generators
        """
        self.x_coord = np.array([])
        self.y_coord = np.array([])
        self.weight = np.array([])
        self.dim = (0, 0)

    def __len__(self):
        """len()
        
        Parameters
        ----------

        Returns
        -------
        int :
            number of generators
        """
        return len(self.weight)

    def randomize(self,
                  dim=(1000, 1000),
                  num=100,
                  weight_range=(1.0, 5.0),
                  min_distance=1.0,
                  padding=0.125):
        """Initialize a set of generators randomly

        Parameters
        ----------
        dim : (int, int)
            upper bounds for the position of any generator
            (Default value = (1000, 1000))

        num : int
            number of generators (Default value = 100)

        weight_range : (float, float)
            range of generator weights (Default value = (1.0, 5.0))

        min_distance : float
            minimum proportion of the mean weight of between any two generators
            (Default value = 1.0)

        padding : float
            mininmum proportion of dimension at borders not to place generators in
            (Default value = 0.125)

        Returns
        -------
        None

        """

        self.dim = dim

        temp_x = []
        temp_y = []
        temp_r = []

        while len(temp_x) < num:
            x_coord = np.random.uniform(0, self.dim[0])
            y_coord = np.random.uniform(0, self.dim[1])
            weight = np.random.uniform(*weight_range)

            if (x_coord < self.dim[0] * padding + weight
                    or self.dim[0] - self.dim[0] * padding - weight < x_coord
                    or y_coord < self.dim[1] * padding + weight
                    or self.dim[1] - self.dim[1] * padding - weight < y_coord):
                continue

            for temp_xi, temp_yi, temp_ri in zip(temp_x, temp_y, temp_r):
                if (abs(temp_xi - x_coord) < temp_ri + weight
                        or abs(temp_yi - y_coord) < temp_ri + weight
                        or np.sqrt((temp_xi - x_coord)**2 +
                                   (temp_yi - y_coord)**2) -
                    (weight + temp_ri) <
                    (weight + temp_ri) * 0.5 * min_distance):
                    break

            else:
                temp_x.append(x_coord)
                temp_y.append(y_coord)
                temp_r.append(weight)

        self.x_coord = np.array(temp_x)
        self.y_coord = np.array(temp_y)
        self.weight = np.array(temp_r)

    def randomize_weights(self, weight_range=(1.0, 5.0), min_distance=1.0):
        """Randomize weights avoiding overlaps

        Parameters
        ----------

        weight_range : (float, float)
            range of generator weights (Default value = (1.0, 5.0))

        min_distance : float
            minimum proportion of the mean weight of between any two generators
            (Default value = 1.0)

        Returns
        -------
        None

        """

        weights = []
        index = 0
        overlap = False
        while len(weights) < len(self):
            weight = np.random.uniform(*weight_range)

            if (self.x_coord[len(weights)] < weight
                    or self.dim[0] - weight < self.x_coord[len(weights)]
                    or self.y_coord[len(weights)] < weight
                    or self.dim[1] - weight < self.y_coord[len(weights)]):
                overlap = True
            else:
                overlap = False

            if not overlap:
                for i in range(len(self)):
                    if i != len(weights):
                        if (abs(self.x_coord[i] - self.x_coord[len(weights)]) <
                                self.weight[i] + weight
                                or abs(self.y_coord[i] -
                                       self.y_coord[len(weights)]) <
                                self.weight[i] + weight
                                or np.sqrt((self.x_coord[i] -
                                            self.x_coord[len(weights)])**2 +
                                           (self.y_coord[i] -
                                            self.y_coord[len(weights)])**2) -
                            (self.weight[i] + weight) <
                            (self.weight[i] + weight) * 0.5 * min_distance):
                            overlap = True
                            break
                    else:
                        overlap = False

            if not overlap:
                weights.append(weight)

        self.weight = np.array(weights)

    def load(self, file_name):
        """Import the set of generators from a .csv file

        Each is assumed to specify the x-coordinate, y-coordinate 
        and weight of one generator

        Parameters
        ----------
        file_name : str
            file location of a .csv file to load the generators from

        Returns
        -------
        None

        """
        temp_x = []
        temp_y = []
        temp_r = []

        if not file_name.endswith(".csv"):
            file_name = os.path.splitext(file_name)[0] + ".csv"

        with open(file_name, "r") as csv_file:
            csv_reader = csv.reader(csv_file)
            for (x_entry, y_entry, r_entry) in csv_reader:
                temp_x.append(float(x_entry))
                temp_y.append(float(y_entry))
                temp_r.append(float(r_entry))

        self.x_coord = np.array(temp_x)
        self.y_coord = np.array(temp_y)
        self.weight = np.array(temp_r)

        self.dim = (np.amax(self.x_coord) +
                    self.weight[np.argmax(self.x_coord)],
                    np.amax(self.y_coord) +
                    self.weight[np.argmax(self.y_coord)])

    def save(self, file_name, overwrite=True):
        """Export the generators to a .csv file

        Parameters
        ----------
        file_name : str
            file location to save the generators to

        overwrite : bool
            overwrite the specified file location should it exist if True,
            append suffix to file name otherwise (Default value = True)

        Returns
        -------
        str
            file location that generators were saved to

        """
        if not file_name.endswith(".csv"):
            file_name = os.path.splitext(file_name)[0] + ".csv"

        if not overwrite and os.path.isfile(file_name):
            num = 0
            while os.path.isfile(file_name):
                num += 1
                file_name = os.path.splitext(file_name)[0] + "_{}".format(
                    num) + ".csv"

        with open(file_name, "w") as file:
            writer = csv.writer(file)
            for i in range(len(self)):
                writer.writerow([
                    str(self.x_coord[i]),
                    str(self.y_coord[i]),
                    str(self.weight[i])
                ])
        return file_name

    def plot(self, color="black", ax=None):
        """Create a Matplotlib axis resembling the generators on a given figure

        Parameters
        ----------
        color : Matplotlib color
            generator color (Default value = "black")

        ax : Matplotlib axis
            Axis to plot on, falls back on current ax for None (Default value = None)

        Returns
        -------
        None

        """
        if not ax:
            plt_ax = plt.gca()
        else:
            plt_ax = ax

        plt_ax.set_aspect("equal")
        plt_ax.axis("off")

        circles = [
            plt.Circle((xi, yi), radius=ri)
            for xi, yi, ri in zip(self.x_coord, self.dim[1] -
                                  self.y_coord, self.weight)
        ]
        if color in plt.colormaps():
            collection = mpl.collections.PatchCollection(circles, cmap=color)
            collection.set_array(self.generators.weight)
        else:
            collection = mpl.collections.PatchCollection(circles, color=color)
        plt_ax.add_collection(collection)

    def compute_from_image(self,
                           file_name,
                           intensity_thresh=(83.0, 90.0),
                           size_thresh=65.0,
                           channel="B"):
        """Determine a set of generators from an image

        Parameters
        ----------
        file_name : str
            file location of image to determine the generators from

        intensity_thresh : (float, float)
            percentile range of rescaled channel intensity below which to
            assign pixel to background and above which to
            assign pixel to nuclei,
            correspond to sensitivity and selectivity, respectively
            (Default = (83.0, 93.0))

        size_thresh : float
            percentile of object size below which to discard a potential
            nucleus as noise (Default = 40.0)

        channel : str, optional
            color channel ("R", "G", "B") of image to use for segmentation,
            falls back to gray scale for other entries but defaults to blue
            for DAPI (Default = "B")

        Returns
        -------
        None

        """
        if channel in ["R", "G", "B"]:
            img = skimage.exposure.rescale_intensity(
                skimage.io.imread(file_name)
                [:, :, ["R", "G", "B"].index(channel)])
        else:
            img = skimage.exposure.rescale_intensity(
                skimage.io.imread(file_name, as_gray=True))
        self.dim = (img.shape[1], img.shape[0])

        markers = np.zeros_like(img)
        markers[img < np.percentile(img, intensity_thresh[0])] = 1
        markers[img > np.percentile(img, intensity_thresh[1])] = 2

        segmentation = skimage.morphology.watershed(skimage.filters.sobel(img),
                                                    markers)

        segmentation = scipy.ndimage.binary_fill_holes(segmentation - 1)
        cell_nucleii, num_labels = scipy.ndimage.label(segmentation)

        nucleii_sizes = np.bincount(cell_nucleii.ravel())
        cell_nucleii[nucleii_sizes[cell_nucleii] < np.percentile(
            nucleii_sizes, size_thresh)] = 0
        cell_nucleii, num_labels = scipy.ndimage.label(cell_nucleii)

        centers_of_mass = scipy.ndimage.measurements.center_of_mass(
            segmentation, cell_nucleii, np.arange(1, num_labels + 1))
        self.x_coord = np.array([col for (row, col) in centers_of_mass])
        self.y_coord = np.array(
            [cell_nucleii.shape[0] - row for (row, col) in centers_of_mass])

        self.weight = np.array([
            np.mean([
                np.amax(np.count_nonzero(cell_nucleii[cell_nucleus], axis=0)),
                np.amax(np.count_nonzero(cell_nucleii[cell_nucleus], axis=1))
            ]) * 0.5
            for cell_nucleus in scipy.ndimage.find_objects(cell_nucleii)
        ])

    def distance(self, generators):
        """
        Map the closest pairs from two sets of generators injectively
        and return their distances

        Parameters
        ----------
        generators : generators
            a set of generators to compare with

        Returns
        -------
        list
            smallest pairwise distances of generators from both sets

        """
        dismat = [[
            max(
                np.sqrt((generators.x_coord[i] - self.x_coord[j])**2 +
                        (generators.y_coord[i] - self.y_coord[j])**2) -
                (generators.weight[i] + self.weight[j]), 0)
            for j in range(len(self))
        ] for i in range(len(generators))]
        min_errors = []

        i_dim = len(generators)
        j_dim = len(self)

        while i_dim > 0 and j_dim > 0:
            minimum = float("inf")
            for i in range(len(generators)):
                for j in range(len(self)):
                    if dismat[i][j] < minimum:
                        minimum = dismat[i][j]
                        min_i = i
                        min_j = j

            min_errors.append(minimum)

            for i in range(len(generators)):
                dismat[i][min_j] = float("inf")

            for j in range(len(self)):
                dismat[min_i][j] = float("inf")

            i_dim -= 1
            j_dim -= 1
        return np.array(min_errors)


class DelaunayTriangulation():
    """Delaunay Trinagulation"""
    def __init__(self):
        """
        Delaunay Triangulation
        """
        self.generators = Generators()
        self.graph = nx.Graph()
        self.positions = {}
        self.graphml = None
        self.root = None
        self.nodes = []
        self.edges = []

    def compute_from_generators(self, generators):
        """Generate the Delaunay Triangulation from a set of generators

        Parameters
        ----------
        generators : generators
            generators to compute the Delaunay Triangulation from

        Returns
        -------
        None

        """
        file_name = "{}{}".format(str(uuid.uuid4()), ".csv")
        self.generators = generators
        self.generators.save(os.path.join(os.path.dirname(__file__),
                                          file_name))
        subprocess.call(
            "./voronoi.exe --delaunay-triangulation {} {} {}".format(
                file_name, self.generators.dim[0], self.generators.dim[1]),
            shell=True,
            cwd=os.path.dirname(__file__))
        file_name = os.path.join(os.path.dirname(__file__), file_name)
        os.remove(file_name)
        file_name = file_name.replace(".csv", ".graphml")
        self.load(file_name)
        os.remove(file_name)

    def load(self, file_name):
        """Import the Delaunay Triangulation from a .graphml file

        Parameters
        ----------
        file_name : str
            file location of a .graphml file to load the
            Delaunay Triangulation from

        Returns
        -------
        None

        """
        if not file_name.endswith(".graphml"):
            file_name = os.path.splitext(file_name)[0] + ".graphml"
        self.graph = nx.readwrite.graphml.read_graphml(file_name)
        self.graphml = etree.parse(file_name.replace(".csv", ".graphml"))
        self.root = self.graphml.getroot()
        self.dim = (int(self.root.find("graph", self.root.nsmap).get("dim_x")),
                    int(self.root.find("graph", self.root.nsmap).get("dim_y")))
        self.positions = {
            node: np.array([node_data["x"], node_data["y"]])
            for node, node_data in self.graph.node.items()
        }
        keys = {}
        for key in self.root.findall("key", self.root.nsmap):
            keys[key.get("id")] = {
                "attr.name": key.get("attr.name"),
                "for": key.get("for"),
                "attr.type": key.get("attr.type")
            }

        max_num_angles = max([
            int(keys[key]["attr.name"].replace("angle", "")) for key in keys
            if keys[key]["attr.name"].startswith("angle")
        ])

        self.nodes = [{
            "position":
            (self.graph.nodes[node]["x"], self.graph.nodes[node]["y"]),
            "angles":
            tuple(self.graph.nodes[node]["angle{}".format(j)]
                  for j in range(max_num_angles + 1)
                  if self.graph.nodes[node].get("angle{}".format(j)))
        } for node in self.graph.nodes]

        self.edges = [{
            "nodes": set([edge[0], edge[1]]),
            "length": self.graph.edges[edge]["length"]
        } for edge in self.graph.edges]

    def save(self, file_name, overwrite=True):
        """Export the Delaunay Triangulation to a .graphml file

        Parameters
        ----------
        file_name : str
            file location to save the Delaunay Triangulation to

        overwrite : bool
            overwrite the specified file location should it exist if True,
            append suffix to file name otherwise (Default value = True)

        Returns
        -------
        str
            file location that Delaunay Triangulation was saved to

        """
        if not file_name.endswith(".graphml"):
            file_name = os.path.splitext(file_name)[0] + ".graphml"

        if not overwrite and os.path.isfile(file_name):
            num = 0
            while os.path.isfile(file_name):
                num += 1
                file_name = os.path.splitext(file_name)[0] + "_{}".format(
                    num) + ".graphml"

        with open(file_name, "wb") as graphml_file:
            self.graphml.write(graphml_file,
                               xml_declaration=True,
                               encoding=self.graphml.docinfo.encoding)
        return file_name

    def plot(self,
             include_generators=True,
             edge_color="black",
             generator_color="black",
             edge_style="solid",
             edge_width=1.0,
             ax=None):
        """Create a Matplotlib axis resembling the Delaunay Triangulation

        Parameters
        ----------
        edge_color : Matplotlib color
            edge color (Default value = "black")

        generator_color :
            generator color (Default value = "black")

        include_generators : bool
            if True, include generators the Delaunay Triangulation
            was generated from in the plot (Default value = True)

        edge_style : Matplotlib line style
            line style for edges (Default value = "solid")

        edge_width : float
            line width for edges

        ax : Matplotlib axis
            Axis to plot on, falls back on current ax for None (Default value = None)

        Returns
        -------
        None

        """
        if not ax:
            plt_ax = plt.gca()
        else:
            plt_ax = ax

        plt_ax.set_aspect("equal")
        plt.box(False)

        nx.draw_networkx(self.graph, {
            node: (x, self.dim[1] - y)
            for node, (x, y) in self.positions.items()
        },
                         node_size=0,
                         with_labels=False,
                         edge_color=edge_color,
                         width=edge_width,
                         ax=plt_ax,
                         style=edge_style)
        if self.generators and include_generators:
            circles = [
                plt.Circle((xi, yi), radius=ri) for xi, yi, ri in zip(
                    self.generators.x_coord, self.dim[1] -
                    self.generators.y_coord, self.generators.weight)
            ]

            if generator_color in plt.colormaps():
                generator_collection = mpl.collections.PatchCollection(
                    circles, cmap=generator_color)
                generator_collection.set_array(self.generators.weight)
            else:
                generator_collection = mpl.collections.PatchCollection(
                    circles, color=generator_color)
            plt_ax.add_collection(generator_collection)


class VoronoiTessellation():
    """Voronoi Tessellation"""
    def __init__(self):
        """
        Voronoi Tessellation
        """
        self.generators = Generators()
        self.graph = nx.Graph()
        self.positions = {}
        self.graphml = None
        self.root = None
        self.nodes = []
        self.edges = []
        self.faces = []
        self.dim = (0, 0)

    def compute_from_generators(self, generators):
        """Generate the weighted Voronoi Tessellation from a set of generators

        Parameters
        ----------
        generators : generators
            generators to compute the Voronoi Tessellation from

        Returns
        -------
        None

        """
        file_name = "{}{}".format(str(uuid.uuid4()), ".csv")
        self.generators = generators
        self.generators.save(os.path.join(os.path.dirname(__file__),
                                          file_name))
        subprocess.call("./voronoi.exe --voronoi-tessellation {} {} {}".format(
            file_name, self.generators.dim[0], self.generators.dim[1]),
                        shell=True,
                        cwd=os.path.dirname(__file__))
        file_name = os.path.join(os.path.dirname(__file__), file_name)
        os.remove(file_name)
        file_name = file_name.replace(".csv", ".graphml")
        self.load(file_name)
        os.remove(file_name)

    def compute_from_image(self,
                           file_name,
                           intensity_thresh=75.0,
                           size_thresh=70.0,
                           node_proximity=5,
                           channel="G"):
        """Determine a Voronoi Tessellation from an image

        Parameters
        ----------
        file_name : str
            file location of image to determine the generators from

        intensity_thresh : float
            percentile of channel intensity above which to assign pixel
            to membrane (Default = 75.0)

        size_thresh : float
            percentile of intermembrane objects' size below which to discard
            the face (Default = 70.0)

        node_proximity : int
            minimum number of pixels in between two separate nodes,
            corresponds to resolution of graph structure given
            the connectivity based skeletonization of pixels 
            assigned to membrane
            (Default = 5)

        channel : str, optional
            color channel ("R", "G", "B") of image to use for segmentation,
            falls back to gray scale, but defaults to green for
            immunofluorescence staining (Default = "G")

        Returns
        -------
        float
            mean deviation between the determined Voronoi Tessellation's vertices
            and those of the Voronoi Tessellation resulting from
            the approximated generators

        """
        if channel in ["R", "G", "B"]:
            img = skimage.exposure.rescale_intensity(
                skimage.io.imread(file_name)
                [:, :, ["R", "G", "B"].index(channel)])
        else:
            img = skimage.exposure.rescale_intensity(
                skimage.io.imread(file_name, as_gray=True))
        self.dim = (img.shape[1], img.shape[0])

        cell_membrane = np.zeros_like(img)
        cell_membrane[img > np.percentile(img, intensity_thresh)] = 1
        cell_membrane, num_labels = scipy.ndimage.label(cell_membrane)
        object_sizes = np.bincount(cell_membrane.ravel())
        cell_membrane[cell_membrane != np.argmax(object_sizes[1:]) + 1] = 0
        intermembrane_space = 1 - cell_membrane

        intermembrane_space, num_labels = scipy.ndimage.label(
            intermembrane_space)
        object_sizes = np.bincount(intermembrane_space.ravel())
        intermembrane_space[object_sizes[intermembrane_space] < np.percentile(
            object_sizes, size_thresh)] = 0
        intermembrane_space[intermembrane_space > 0] = 1
        cell_membrane = 1 - intermembrane_space

        cell_membrane = skimage.morphology.skeletonize(cell_membrane).astype(
            int)
        kernel = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]])
        num_neighbor_pixels = scipy.ndimage.convolve(cell_membrane,
                                                     kernel,
                                                     mode="constant",
                                                     cval=0.0)
        num_neighbor_pixels[cell_membrane == 0] = 0

        branching_pixels = np.zeros_like(num_neighbor_pixels)
        branching_pixels[num_neighbor_pixels > 2] = 1

        branching_positions = np.nonzero(branching_pixels)
        branching_positions = [
            (x, y)
            for x, y in zip(branching_positions[0], branching_positions[1])
        ]
        adj = [[False for j in branching_positions]
               for i in branching_positions]

        def dfs(start, current, path):
            if current != start and branching_pixels[current[0], current[1]]:
                adj[branching_positions.index(start)][
                    branching_positions.index(current)] = True
            else:
                for i in (-1, 0, 1):
                    for j in (-1, 0, 1):
                        if ((i != 0 or j != 0)
                                and (0 <= current[0] + i <
                                     num_neighbor_pixels.shape[0])
                                and (0 <= current[1] + j <
                                     num_neighbor_pixels.shape[1])
                                and num_neighbor_pixels[current[0] +
                                                        i, current[1] + j] and
                            (current[0] + i, current[1] + j) not in path):
                            dfs(start, (current[0] + i, current[1] + j),
                                path + [current])

        for pos in branching_positions:
            dfs(pos, pos, [])

        updated_adj = True
        while updated_adj:
            updated_adj = False
            for i in range(len(branching_positions)):
                for j in range(len(branching_positions)):
                    if i != j:
                        if (abs(branching_positions[i][0] -
                                branching_positions[j][0]) < node_proximity
                                and abs(branching_positions[i][1] -
                                        branching_positions[j][1]) <
                                node_proximity):
                            for row in range(len(adj)):
                                adj[row][i] = adj[row][i] or adj[row][j]
                            for col in range(len(adj)):
                                adj[i][col] = adj[i][col] or adj[j][col]
                            adj[i][i] = False
                            for row in range(len(adj)):
                                del adj[row][j]
                            del adj[j]
                            del branching_positions[j]
                            updated_adj = True
                            break
                if updated_adj:
                    break

        self.graph.add_nodes_from([i for i in range(len(branching_positions))])

        for i in self.graph.nodes:
            for j in self.graph.nodes:
                if adj[i][j]:
                    self.graph.add_edge(i, j)

        for node in self.graph.nodes:
            self.positions[node] = (branching_positions[node][1],
                                    img.shape[0] -
                                    branching_positions[node][0])

        while [
                node for node in self.graph.nodes
                if self.graph.degree[node] < 2
        ]:
            nodes_to_remove = [
                node for node in self.graph.nodes
                if self.graph.degree[node] < 2
            ]
            self.graph.remove_nodes_from(nodes_to_remove)
            for node in nodes_to_remove:
                del self.positions[node]

        reenumeration = {
            node: str(list(self.graph.nodes).index(node))
            for node in self.graph.nodes
        }

        self.graph = nx.relabel_nodes(self.graph, reenumeration)
        self.positions = {
            reenumeration[node]: position
            for node, position in self.positions.items()
        }

        for node, (x, y) in self.positions.items():
            self.graph.node[node]['x'] = x
            self.graph.node[node]['y'] = y

        self.root = etree.fromstring(
            chr(10).join(nx.generate_graphml(self.graph)))
        self.graphml = etree.ElementTree(self.root)
        self.generators = Generators()

    def approximate_generators(self, epsilon=0.0):
        """Approximate a set of generators given the Voronoi Tessellation
        and compute face data

        Parameters
        ----------
        epsilon : float
            threshold for regression (Default 0.0). If the error of the vertices can
            consistently not be reduced by more than epsilon, terminate.
            Epsilon := 0.0 for best approximation, runtime decreases with increasing epsilon.

        Returns
        -------
        generators
            approximated generators

        ndarray
            mean errors resulting from accepted translations

        ndarray
            mean errors resulting from any translations


        """
        file_name_base = str(uuid.uuid4())
        file_name = "{}{}".format(file_name_base, ".graphml")
        file_name = os.path.join(os.path.dirname(__file__), file_name)
        self.save(os.path.join(file_name))
        subprocess.call("./voronoi.exe --generators {} {}".format(
            file_name, epsilon),
                        shell=True,
                        cwd=os.path.dirname(__file__))

        self.load(file_name)
        os.remove(file_name)

        file_name = file_name.replace(".graphml", ".csv")
        generators = Generators()
        generators.load(file_name)
        os.remove(file_name)
        generators.dim = self.dim

        with open(
                os.path.join(
                    os.path.dirname(__file__),
                    "{}{}".format(file_name_base,
                                  "_accepted_errors.csv"))) as acc_err_file:
            accepted_errors = np.array([
                (float(line.split(",")[0]), float(line.split(",")[1]))
                for line in acc_err_file.read().split()
            ])
        os.remove(
            os.path.join(os.path.dirname(__file__),
                         "{}{}".format(file_name_base,
                                       "_accepted_errors.csv")))

        with open(
                os.path.join(os.path.dirname(__file__),
                             "{}{}".format(file_name_base,
                                           "_errors.csv"))) as err_file:
            errors = np.array([float(err) for err in err_file.read().split()])
        os.remove(
            os.path.join(os.path.dirname(__file__),
                         "{}{}".format(file_name_base, "_errors.csv")))

        return generators, accepted_errors, errors

    def face_centroids(self, weight=1.0):
        """Return face centroids as a set of generators given the Voronoi Tessellation.
        As it is necessary to compute faces for this, the Voronoi Tessellation is
        updated should face data not be available. Otherwise the face data is
        simply exported to a set of generators.

        Parameters
        ----------
        weight : float
            uniform weight of centroids

        Returns
        -------
        generators
            centroid generators

        """
        if self.faces:
            centroids = Generators()
            centroids.x_coord = np.array(
                [face["centroid"][0] for face in self.faces])
            centroids.y_coord = np.array(
                [face["centroid"][1] for face in self.faces])
            centroids.weight = np.array([weight for face in self.faces])
            centroids.dim = self.dim
            return centroids

        file_name_base = str(uuid.uuid4())
        file_name = "{}{}".format(file_name_base, ".graphml")
        file_name = os.path.join(os.path.dirname(__file__), file_name)
        self.save(os.path.join(file_name))
        subprocess.call("./voronoi.exe --centroids {}".format(file_name),
                        shell=True,
                        cwd=os.path.dirname(__file__))
        self.load(file_name)
        os.remove(file_name)

        file_name = file_name.replace(".graphml", ".csv")
        centroids = Generators()
        centroids.load(file_name)
        os.remove(file_name)
        centroids.dim = self.dim
        centroids.weight = np.array([weight for i in range(len(generators))])

        return centroids

    def distance(self, voronoi_tessellation):
        """
        Map the closest pairs from two sets of vertices of a Voronoi Tessellation injectively
        and return their distances

        Parameters
        ----------
        voronoi_tessellation : Voronoi Tessellation
            a Voronoi Tessellation to compare with

        Returns
        -------
        list
            smallest pairwise distances of generators from both sets

        """
        dismat = [[
            np.sqrt((voronoi_tessellation_node["position"][0] -
                     node["position"][0])**2 +
                    (voronoi_tessellation_node["position"][1] -
                     node["position"][1])**2) for node in self.nodes
        ] for voronoi_tessellation_node in voronoi_tessellation.nodes]
        min_errors = []

        i_dim = len(voronoi_tessellation.nodes)
        j_dim = len(self.nodes)

        while i_dim > 0 and j_dim > 0:
            minimum = float("inf")
            for i in range(len(voronoi_tessellation.nodes)):
                for j in range(len(self.nodes)):
                    if dismat[i][j] < minimum:
                        minimum = dismat[i][j]
                        min_i = i
                        min_j = j

            min_errors.append(minimum)

            for i in range(len(voronoi_tessellation.nodes)):
                dismat[i][min_j] = float("inf")

            for j in range(len(self.nodes)):
                dismat[min_i][j] = float("inf")

            i_dim -= 1
            j_dim -= 1
        return np.array(min_errors)

    def load(self, file_name):
        """Import the Voronoi )Tessellation from a .graphml file

        Parameters
        ----------
        file_name : str
            file location of a .graphml file to load the Voronoi Tessellation from

        Returns
        -------
        None

        """
        if not file_name.endswith(".graphml"):
            file_name = os.path.splitext(file_name)[0] + ".graphml"
        self.graph = nx.readwrite.graphml.read_graphml(file_name)
        self.graphml = etree.parse(file_name.replace(".csv", ".graphml"))
        self.root = self.graphml.getroot()
        self.dim = (int(self.root.find("graph", self.root.nsmap).get("dim_x")),
                    int(self.root.find("graph", self.root.nsmap).get("dim_y")))
        self.positions = {
            node: np.array([node_data["x"], node_data["y"]])
            for node, node_data in self.graph.node.items()
        }
        keys = {}
        for key in self.root.findall("key", self.root.nsmap):
            keys[key.get("id")] = {
                "attr.name": key.get("attr.name"),
                "for": key.get("for"),
                "attr.type": key.get("attr.type")
            }

        self.faces = []
        for face in self.root.find("graph", self.root.nsmap).findall(
                "face", self.root.nsmap):
            face_nodes = []
            for node in face.findall("node", self.root.nsmap):
                face_nodes.append(node.get("node"))
            face_area = 0.0
            face_centroid = [0.0, 0.0]
            for data in face.findall("data", self.root.nsmap):
                if keys[data.get("key")]["for"] == "face":
                    if (keys[data.get("key")]["attr.name"] == "area" and
                            keys[data.get("key")]["attr.type"] == "double"):
                        face_area = float(data.text)
                    elif (keys[data.get("key")]["attr.name"] == "centroid_x"
                          and keys[data.get("key")]["attr.type"] == "double"):
                        face_centroid[0] = float(data.text)
                    elif (keys[data.get("key")]["attr.name"] == "centroid_y"
                          and keys[data.get("key")]["attr.type"] == "double"):
                        face_centroid[1] = float(data.text)

            self.faces.append({
                "nodes": tuple(face_nodes),
                "area": face_area,
                "centroid": tuple(face_centroid)
            })

        max_num_angles = max([
            int(keys[key]["attr.name"].replace("angle", "")) for key in keys
            if keys[key]["attr.name"].startswith("angle")
        ])

        self.nodes = [{
            "position":
            (self.graph.nodes[node]["x"], self.graph.nodes[node]["y"]),
            "angles":
            tuple(self.graph.nodes[node]["angle{}".format(j)]
                  for j in range(max_num_angles + 1)
                  if self.graph.nodes[node].get("angle{}".format(j)))
        } for node in self.graph.nodes]

        self.edges = [{
            "nodes": set([edge[0], edge[1]]),
            "length": self.graph.edges[edge]["length"]
        } for edge in self.graph.edges]

    def save(self, file_name, overwrite=True):
        """Export the Voronoi Tessellation to a .graphml file

        Parameters
        ----------
        file_name : str
            file location to save the Voronoi Tessellation to

        overwrite : bool
            overwrite the specified file location should it exist if True,
            append suffix to file name otherwise (Default value = True)

        Returns
        -------
        str
            file location that Voronoi Tessellation was saved to

        """
        if not file_name.endswith(".graphml"):
            file_name = os.path.splitext(file_name)[0] + ".graphml"

        if not overwrite and os.path.isfile(file_name):
            num = 0
            while os.path.isfile(file_name):
                num += 1
                file_name = os.path.splitext(file_name)[0] + "_{}".format(
                    num) + ".graphml"

        with open(file_name, "wb") as graphml_file:
            self.graphml.write(graphml_file,
                               xml_declaration=True,
                               encoding=self.graphml.docinfo.encoding)
        return file_name

    def plot(self,
             include_generators=True,
             edge_color="black",
             generator_color="black",
             edge_style="solid",
             edge_width=1.0,
             ax=None):
        """Create a Matplotlib axis resembling the Voronoi Tessellation

        Parameters
        ----------
        include_generators : bool
            if True, include generators the Voronoi Tessellation was generated
            from in the plot (Default value = True)

        edge_color : Matplotlib color
            edge color (Default value = "black")

        generator_color : Matplotlib color
            generator color (Default value = "black")

        edge_style : Matplotlib line style
            line style for edges (Default value = "solid")

        edge_width : float
            line width for edges

        ax : Matplotlib axis
            Axis to plot on, falls back on current ax for None (Default value = None)

        Returns
        -------
        None

        """
        if not ax:
            plt_ax = plt.gca()
        else:
            plt_ax = ax

        plt_ax.set_aspect("equal")
        plt.box(False)

        nx.draw_networkx(self.graph, {
            node: (x, self.dim[1] - y)
            for node, (x, y) in self.positions.items()
        },
                         node_size=0,
                         with_labels=False,
                         edge_color=edge_color,
                         width=edge_width,
                         ax=plt_ax,
                         style=edge_style)
        if self.generators and include_generators:
            circles = [
                plt.Circle((xi, yi), radius=ri) for xi, yi, ri in zip(
                    self.generators.x_coord, self.dim[1] -
                    self.generators.y_coord, self.generators.weight)
            ]

            if generator_color in plt.colormaps():
                generator_collection = mpl.collections.PatchCollection(
                    circles, cmap=generator_color)
                generator_collection.set_array(self.generators.weight)

            else:
                generator_collection = mpl.collections.PatchCollection(
                    circles, color=generator_color)
            plt_ax.add_collection(generator_collection)


class CircularVoronoiTessellation():
    """Generalized Voronoi Tessellation yielding a morphology 
    widely similar to that of keratinocytes in epithelial cell layers

    "Generalized Voronoi Tessellation as a Model of
    Two-dimensional Cell Tissue Dynamics", Bock et al. (2009)

    """
    def __init__(self):
        """
        Voronoi Tessellation
        """
        self.generators = Generators()
        self.graph = nx.Graph()
        self.positions = {}
        self.faces = {}
        self.graphml = None
        self.root = None
        self.dim = (0, 0)

    def compute_from_generators(self, generators):
        """Generate the weighted Voronoi Tessellation from a set of generators

        Parameters
        ----------
        generators : generators
            generators to compute the Voronoi Tessellation from

        Returns
        -------
        None

        """
        file_name = "{}{}".format(str(uuid.uuid4()), ".csv")
        self.generators = generators
        self.generators.save(os.path.join(os.path.dirname(__file__),
                                          file_name))
        subprocess.call(
            "./voronoi.exe --circular-voronoi-tessellation {} {} {}".format(
                file_name, self.generators.dim[0], self.generators.dim[1]),
            shell=True,
            cwd=os.path.dirname(__file__))
        file_name = os.path.join(os.path.dirname(__file__), file_name)
        os.remove(file_name)
        file_name = file_name.replace(".csv", ".graphml")
        self.load(file_name)
        os.remove(file_name)

    def load(self, file_name):
        """Import the Voronoi Tessellation from a .graphml file

        Parameters
        ----------
        file_name : str
            file location of a .graphml file to load the Voronoi Tessellation from

        Returns
        -------
        None

        """
        if not file_name.endswith(".graphml"):
            file_name = os.path.splitext(file_name)[0] + ".graphml"
        if os.path.isfile(file_name.replace(".graphml", ".csv")):
            self.generators.load(file_name.replace(".graphml", ".csv"))
        self.graph = nx.readwrite.graphml.read_graphml(file_name)
        self.graphml = etree.parse(file_name.replace(".csv", ".graphml"))
        self.root = self.graphml.getroot()
        self.positions = {
            node: np.array([node_data["x"], node_data["y"]])
            for node, node_data in self.graph.node.items()
        }
        self.dim = (int(self.root.find("graph", self.root.nsmap).get("dim_x")),
                    int(self.root.find("graph", self.root.nsmap).get("dim_y")))

    def save(self, file_name, overwrite=True):
        """Export the Voronoi Tessellation to a .graphml file

        Parameters
        ----------
        file_name : str
            file location to save the Voronoi Tessellation to

        overwrite : bool
            overwrite the specified file location should it exist if True,
            append suffix to file name otherwise (Default value = True)

        Returns
        -------
        str
            file location that Voronoi Tessellation was saved to

        """
        if not file_name.endswith(".graphml"):
            file_name = os.path.splitext(file_name)[0] + ".graphml"

        if not overwrite and os.path.isfile(file_name):
            num = 0
            while os.path.isfile(file_name):
                num += 1
                file_name = os.path.splitext(file_name)[0] + "_{}".format(
                    num) + ".graphml"

        with open(file_name, "wb") as graphml_file:
            self.graphml.write(graphml_file,
                               xml_declaration=True,
                               encoding=self.graphml.docinfo.encoding)
        return file_name

    def plot(self,
             include_generators=True,
             edge_color="black",
             generator_color="black",
             edge_style="solid",
             edge_width=1.0,
             ax=None):
        """Create a Matplotlib axis resembling the Voronoi Tessellation

        Parameters
        ----------
        include_generators : bool
            if True, include generators the Voronoi Tessellation was generated from
            in the plot (Default value = True)

        edge_color : Matplotlib color
            edge color (Default value = "black")

        generator_color : Matplotlib color
            generator color (Default value = "black")

        edge_style : Matplotlib line style
            line style for edges (Default value = "solid")

        edge_width : float
            line width for edges

        ax : Matplotlib axis
            Axis to plot on, falls back on current ax for None (Default value = None)

        Returns
        -------
        None

        """
        if not ax:
            plt_ax = plt.gca()
        else:
            plt_ax = ax

        plt_ax.set_aspect("equal")
        plt.box(False)

        nx.draw_networkx(self.graph, {
            node: (x, self.dim[1] - y)
            for node, (x, y) in self.positions.items()
        },
                         node_size=0,
                         with_labels=False,
                         edge_color=edge_color,
                         width=edge_width,
                         ax=plt_ax,
                         style=edge_style)
        if self.generators and include_generators:
            circles = [
                plt.Circle((xi, yi), radius=ri) for xi, yi, ri in zip(
                    self.generators.x_coord, self.dim[1] -
                    self.generators.y_coord, self.generators.weight)
            ]

            if generator_color in plt.colormaps():
                generator_collection = mpl.collections.PatchCollection(
                    circles, cmap=generator_color)
                generator_collection.set_array(self.generators.weight)
            else:
                generator_collection = mpl.collections.PatchCollection(
                    circles, color=generator_color)
            plt_ax.add_collection(generator_collection)
